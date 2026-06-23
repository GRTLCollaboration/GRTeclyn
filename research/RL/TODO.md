This is one of the most notorious "silent killers" in Numerical Relativity. It
is the exact reason why even the most advanced supercomputer simulations of
binary black holes (like those from the SXS Collaboration) eventually stop after
the post-merger ringdown.

Sommerfeld boundaries mathematically assume that all outgoing radiation is
perfectly spherical (f(r - ct) / r) and traveling exactly perpendicular to the
boundary face. But a rotating gravitational motor produces quadrupolar waves.
When those hit a cubic grid boundary, they shear. The boundary math fails, and
high-frequency "junk" reflects back inward at the speed of light.

If you punish your RL agent because the simulation crashed at t=500 due to
boundary reflections, you are punishing the agent for a C++ math limitation it
cannot control. The agent will become paranoid and stop pumping energy entirely
just to avoid generating waves that might bounce back.

Here is how we engineer around the Sandbox Problem as Senior Developers.

1. Calculate the "Causal Clean Window"

We don't guess the truncation time; we derive it from the speed of light (c=1 in
code units) and your grid dimensions.

In your params.txt, your grid is L_full = 64.0, and your center is
32.0, 32.0, 0.0.

  - The distance from the motor to the nearest boundary is 32 units.
  - It takes t=32 for the first gravitational wave to hit the boundary.
  - It takes another t=32 for the corrupted reflection to bounce back and hit
    the center.

This means your simulation is perfectly "clean" and causally isolated from
boundary errors until t=64. After t=64, the agent is swimming in its own echoes.

2. The RL Fix: The Truncation Horizon

We set the Gym environment's maximum lifespan to just after the first echo
returns.

In env.py:

# The agent crashed the physics itself
terminated = bool(min_lapse < 0.05 or L2_Ham > 0.05)

# The Sandbox filled with boundary entropy. 
# Stop the run before the noise kills the simulation.
truncated = bool(sim_time >= 80.0) 

Why truncated vs terminated is critical in PPO

In reinforcement learning, there is a massive mathematical difference between
these two flags:

  - terminated = True (Death): The episode ended because the agent drove off a
    cliff. The Value Function is forced to 0. The agent learns: "Whatever I just
    did was fatal. Never do that again."
  - truncated = True (Time Limit): The episode ended because the simulation ran
    out of time, but the agent was doing great. Stable-Baselines3 handles this
    using Time-Limit Bootstrapping. It looks at the Critic Network's prediction
    for t=81 and adds it to the reward. The agent learns: "I was doing an
    awesome job, and if the gods hadn't paused the universe, I would have kept
    surviving forever."

3. The Kreiss-Oliger (KO) Sponge (C++ Mitigation)

If you eventually want the agent to run to t=200 to study long-term motor
thermodynamics, you have to damp the echoes.

You already have this in your codebase! It is the sigma = 2.0 parameter in your
params.txt (Kreiss-Oliger dissipation). It acts as a low-pass filter, killing
high-frequency grid noise while letting the smooth physics survive.

If boundary entropy becomes a problem even before t=80, you can add a Sponge
Layer to the C++ RHS: artificially increase sigma near the boundaries of the
grid, so the waves get absorbed like acoustic foam before they can reflect.

The Verdict on the Continuous Run

By using truncated at the causal light-crossing time (t \approx 80), we achieve
your vision: The agent continuously updates its policy without being penalized
by t=0 resets, it learns an infinite-horizon motor strategy, and we
systematically shield it from being blamed for the limitations of numerical
boundary conditions.

This is how you build a rigorous, physics-aware ML pipeline.


Ah, I see the exact nuance now. That is a critical distinction.

In a passive Initial Value Problem (IVP), the peak physical performance of a
high-potential genome is almost always a transient intermediate state (for
example, eval 063 hitting its peak of 19.3\% FTL at t=16 before dispersing, or
eval 046 peaking before hitting a horizon at t \approx 21).

If you only look at the "final" genomes as they exist at the end of the campaign
(t=30), you are looking at the parameter sets that survived, which often
represent passive, low-energy, uninteresting dynamics.

The Correct Strategy: The "High-Peak" Genome Extractor

Instead of dynamic midway checkpointing or final-state seeding, the correct
approach is to scan the historical trajectory telemetry of your MAP-Elites /
CMA-ES runs, find the t=0 genomes that produced the highest intermediate peak
scores, and use those as your starting pool.

The RL agent’s mission is clear: It starts at t=0 with a genome that is proven
to dynamically open a massive warp throat at t \approx 16. The agent lets the
simulation evolve, and as the system approaches that critical, unstable peak
phase, the agent actively drives the pump and gauge to sustain that peak and
prevent the collapse that would normally occur by t=30.

Here is how we implement this precisely in the Python wrapper:

1. The Extraction Pipeline (seed.py)

Instead of manually copying parameters, seed.py will programmatically mine your
existing MAP-Elites/CMA-ES run directories.

# grteclyn_wrapper/rl/seed.py
import json
import pathlib
import numpy as np

def extract_high_peak_genomes(archive_path: str, top_n: int = 50) -> list[dict]:
    """
    Scans the historical runs in a MAP-Elites/CMA-ES archive.
    Finds the t=0 genomes that produced the highest intermediate peak 
    of ftl_geo_evolving at any point during their evolution.
    """
    archive_dir = pathlib.Path(archive_path)
    candidates = []

    # Iterate through all evaluation folders (eval_000001, eval_000002, etc.)
    for eval_dir in archive_dir.glob("eval_*"):
        timeseries_file = eval_dir / "output" / "small_data" / "central_timeseries.dat"
        overrides_file = eval_dir / "overrides.json" # Your t=0 parameters
        
        if not (timeseries_file.exists() and overrides_file.exists()):
            continue
            
        try:
            # Read the step-by-step telemetry
            # (Assuming columns: time, chi, lapse, K, weyl4, ftl_geo_evolving, ...)
            data = np.loadtxt(timeseries_file, skiprows=1)
            if data.ndim == 1:
                data = np.expand_dims(data, axis=0)
                
            # Find the peak value of ftl_geo_evolving (e.g., column 5) at any intermediate step
            ftl_column = data[:, 5]
            peak_ftl = float(np.max(ftl_column))
            
            with open(overrides_file, "r") as f:
                overrides = json.load(f)
                
            candidates.append({
                "eval_id": eval_dir.name,
                "peak_ftl": peak_ftl,
                "overrides": overrides
            })
        except Exception:
            continue

    # Sort descending by the peak intermediate score
    candidates.sort(key=lambda x: x["peak_ftl"], reverse=True)
    
    # Return the top N t=0 genomes
    return [c["overrides"] for c in candidates[:top_n]]

2. Environment Initialization (env.py)

During env.reset(), the environment selects one of these high-potential t=0
genomes, ensuring training always starts on a physically interesting, highly
dynamic chassis.

# grteclyn_wrapper/rl/env.py
from grteclyn_wrapper.rl.seed import extract_high_peak_genomes

class SpacetimeFtlEnv(gym.Env):
    def __init__(self, archive_path, ...):
        # Extract the high-peak seed pool once at startup
        self.seed_pool = extract_high_peak_genomes(archive_path, top_n=50)
        
    def reset(self, seed=None, options=None):
        # Sample a high-potential t=0 genome
        genome = self.np_random.choice(self.seed_pool)
        
        # Write the t=0 overrides to params_rl.txt
        self._write_params_file(genome)
        
        # Launch the C++ simulation from t=0
        self.sim_process = subprocess.Popen(...)
        
        return self._wait_for_first_obs(), {}

Why this is mathematically and computationally elegant

1.  Pure t=0 Initial Value Problem: You preserve perfect Einstein constraint
    satisfaction. The simulation starts clean, and the C++ engine doesn't have
    to deal with the numerical noise of midway checkpoint files.
2.  No Dynamic I/O Overhead: You write a flat text file (params_rl.txt) at
    reset. No multi-gigabyte HDF5 reads are needed during training.
3.  Active Stabilization Training: The agent gets to experience the natural,
    passive buildup of the warp throat from t=0 \to t=16. It learns the
    "momentum" of the system, so it is fully prepared to apply the correct
    counter-measures (the pump and gauge steering) exactly when the instability
    begins to manifest.
4.  No Exploration Bottleneck: You are using your historical MAP-Elites search
    data as a "filter" to discard 99% of the uninteresting phase space. PPO only
    trains on initial geometries that are guaranteed to produce dynamic,
    high-scoring intermediate states.

This completely resolves the final conceptual loop. It is a highly robust,
mathematically rigorous, and sample-efficient way to execute the campaign.
