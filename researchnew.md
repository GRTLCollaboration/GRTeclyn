\documentclass[twocolumn,aps,prd,nofootinbib,natbib]{revtex4-2}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{graphicx}
\usepackage{amssymb}
\usepackage{amsmath}
\usepackage{xcolor}
\usepackage{hyperref}
\usepackage{textcase}

\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    filecolor=magenta,
    urlcolor=blue,
    citecolor=blue,
}
\urlstyle{same}

\begin{document}

\title{\MakeTextUppercase{Neural Spacetime Search with \texttt{GRTeclyn}: A Practical Baseline for Dynamical Metric Discovery}}

\author{\MakeTextUppercase{Nikita M. Shirokov}}
\email{shirokov.nm@phystech.edu}
\affiliation{\textit{Independent Researcher, Moscow, Russia}}
\date{\today}

\begin{abstract}
We outline a practical baseline for turning \texttt{GRTeclyn} into a closed-loop search engine for exotic spacetime metrics on a GPU cluster. Building on prior 3D numerical-relativity evolutions of the Ellis--Bronnikov wormhole, we present an automated framework that goes beyond static metric inspection. A Python scheduler proposes parameter vectors representing spatial and extrinsic curvature profiles. We outline two primary pathways to handle the matter sector: directly reconstructing the required matter distribution from the Einstein constraint equations at $t=0$, or employing an elliptic solver to force the starting geometry to be mathematically consistent with standard, positive-energy matter. The system is designed to run in batches on a GPU cluster, using early-stopping filters to optimize node hours. Promoted candidates are subjected to operational faster-than-light (FTL) tests, a resolution convergence ladder, and physical stability checks to distinguish physical FTL dynamics from coordinate artifacts.
\end{abstract}

\maketitle

\section{Introduction}

Most proposed interstellar metrics are designed by hand. One writes a static, highly symmetric metric ansatz, calculates the resulting stress-energy tensor using the Einstein field equations, and inspects whether it possesses FTL or shortcut properties. This traditional method is biased toward highly simplified, static geometries that humans know how to write down in closed form.

By treating metric discovery as a closed-loop optimization problem, we can search a much wider, non-spherical, and dynamical space of initial data. The previous \texttt{GRTeclyn} wormhole study demonstrated that 3D numerical relativity is highly suited for tracking the nonlinear evolution, instability, and collapse of exotic geometries. 

Rather than inspecting static equations, the proposed framework couples a Python-based optimization engine with the high-performance C++ GPU-accelerated AMR engine of \texttt{GRTeclyn}. This document outlines the physical, mathematical, and computational framework required to run this search engine across a GPU cluster, with a strict emphasis on deriving the matter sector and avoiding coordinate artifacts.

Throughout this draft, ``FTL'' means effective faster-than-light propagation through a curved spacetime geometry compared with flat-spacetime light travel. It does not imply local superluminal motion relative to the local light cone.

\section{Goal}

The objective is to construct an automated numerical-relativity laboratory for exotic spacetime metrics. The pipeline searches for geometries that display four classes of physical behavior:

\begin{enumerate}
    \item \textbf{FTL transport:} a compact passenger region moves effectively faster than a light signal in flat spacetime, while maintaining bounded, survivable tidal forces.
    \item \textbf{FTL communication:} a physical signal (null ray or scalar pulse) arrives at a detector earlier than it would in a flat-spacetime control run.
    \item \textbf{Portals:} two distinct mouth-like structures create a stable shortcut.
    \item \textbf{Wormholes:} a throat-like spatial bridge that survives for multiple light-crossing times without immediately collapsing into a horizon.
\end{enumerate}

To ensure physical viability, the ultimate success criterion is to find stable, non-collapsing FTL geometries that minimize—or entirely eliminate—violations of the standard energy conditions.

\subsection{Optimization as the Engine of Discovery}

A foundational question is why a project aimed at the \emph{discovery} of novel, unclassified spacetimes relies so heavily on numerical \emph{optimization}. The answer lies in the dimensionality and physical constraints of the search space.

The configuration space of all possible 3D initial geometries—defined by the spatial metric $\gamma_{ij}$ and the extrinsic curvature $K_{ij}$—is infinitely dimensional. If one conducts a blind, unguided search (i.e., random sampling), the probability of finding a viable candidate is vanishingly small. Virtually all randomly generated initial metrics fail immediately: they either violate the Einstein constraints, collapse into immediate singularity horizons, or disperse into trivial flat space. 

Optimization provides the steering mechanism necessary to navigate this infinite space. By defining a strict objective scoring function—such as maximizing FTL signal arrival time while minimizing tidal forces and penalizing negative energy density—we establish a physics-guided fitness landscape. 

The optimizer does not merely "tweak" a pre-existing geometry; instead, it acts as an evolutionary sculptor. Using a flexible basis of radial and angular functions (e.g., spherical harmonics) as its raw material, the optimizer dynamically combines these modes. Through this process, the system can discover highly non-trivial, asymmetric, or rotating geometries—such as multi-mouth portals or stabilized warp-like corridors—that have never been written down in closed-form by hand. In this framework, optimization is not merely parameter tuning; it is the mathematical mechanism of topological and geometric discovery.

\section{Initial Data and Matter Derivation}

In general relativity, we cannot prescribe an arbitrary initial geometry without satisfying the Einstein constraint equations at $t=0$. The initial spatial metric $\gamma_{ij}$ and the extrinsic curvature $K_{ij}$ are mathematically coupled to the matter's local energy density $\rho$ and momentum density $j_i$.

We propose two distinct pathways to handle the matter sector within the automated search loop:

\subsection{Path A: Algebraic Matter Reconstruction}

In this pathway, the optimizer is permitted to generate any arbitrary spatial geometry. We do not assume a specific matter model. Instead, we use the Einstein equations to algebraically derive the physical "cost" of the proposed geometry. 

Before the time evolution starts, \texttt{GRTeclyn} evaluates the Hamiltonian and Momentum constraints at every point on the 3D grid. We algebraically reconstruct the required energy density $\rho$ and momentum density $j_i$ at $t=0$:

\begin{equation}
    \rho = \frac{1}{16\pi} \left( R + K^2 - K_{ij}K^{ij} \right),
\end{equation}
\begin{equation}
    j_i = \frac{1}{8\pi} D_j\left(K^j_i - \delta^j_i K\right).
\end{equation}

Here, $R$ is the 3D Ricci scalar and $D_j$ is the covariant derivative. The Python script reads this derived $\rho$ grid from the diagnostic files. If $\rho < 0$ in any region, the script flags that this geometry requires exotic matter and calculates the integrated volume of negative energy.

\subsection{Path B: Elliptic Constraint Solving}

If we wish to strictly forbid exotic matter from the start, we cannot simply guess a geometry and reconstruct $\rho$. Instead, we fix the matter sector to a standard, positive-energy model (such as a standard scalar field with $\rho \ge 0$). 

We then employ an elliptic initial-data constraint solver (such as \texttt{GRTresna}, designed to work with the AMReX grid). The optimizer proposes a "target" geometry, and the elliptic solver slightly deforms this target until it mathematically satisfies the constraints while coupled only to our positive-energy matter field. The resulting constraint-satisfying initial data is then passed to \texttt{GRTeclyn} for evolution.

\section{Candidate Parameterization}

To avoid setting up million-point grids directly from the optimizer, the Python script generates a compact parameter vector $\theta$. A decoder turns $\theta$ into smooth initial fields on the AMReX grid.

\subsection{Level 1: Radial Recipes}

We define the initial fields at $t=0$ using spherically symmetric radial profiles:
\begin{equation}
    q(r) = q_\infty + \sum_n a_n B_n(r),
\end{equation}
where $q$ represents the CCZ4 fields $\chi$, $\alpha$, $\beta^i$, or $K$. The coefficients $a_n$ are the parameters sampled by the optimizer, and $B_n(r)$ are smooth radial basis functions (such as Gaussians).

\subsection{Level 2: Angular Recipes}

To allow for non-spherical, asymmetric, or rotating geometries, we expand the parameterization using spherical harmonics:
\begin{equation}
    q(r,\theta,\varphi) = q_\infty + \sum_{n\ell m} a_{n\ell m}B_n(r)Y_{\ell m}(\theta,\varphi).
\end{equation}
This enables the search engine to discover rotating warp fields, quadrupolar configurations, and asymmetric wormhole mouths.

\subsection{Level 3: Neural Field Decoder}

In the advanced phase, we can parameterize the initial fields using an implicit coordinate-based neural network (such as a SIREN network):
\begin{equation}
    (x,y,z) \xrightarrow{\text{NN}(\theta)} \{\chi,\alpha,\beta^i,K,\tilde A_{ij},\phi,\Pi\}.
\end{equation}
This is the most flexible option but requires a stable training loop and will be implemented only after the Level 1 and Level 2 pipelines are thoroughly validated.

\subsection{Level 4: Hybrid Meta-Optimization via LLM Agents}

While Level 1 and Level 2 parameter searches are managed by derivative-free numerical optimizers (such as CMA-ES), we introduce a high-level LLM-based research agent to manage the meta-optimization loop. The LLM agent does not propose raw floating-point coefficients. Instead, it operates on a feedback loop utilizing the diagnostic output of completed batches:

\begin{enumerate}
    \item \textbf{Automated Diagnostics and Debugging:} If a cluster of candidates terminates early with solver failures or constraint explosions, the LLM agent parses the raw compilation and execution logs. It diagnoses whether the failure is physical (e.g., a physical singularity) or numerical (e.g., coordinate shear, boundary reflection, or grid resolution exhaustion) and modifies the CCZ4 gauge parameters or grid settings accordingly.
    \item \textbf{Dynamic Hypothesis and Basis Selection:} If the numerical optimizer plateaus in a local minimum, the LLM agent generates new physical hypotheses. It instructs the Python scheduler to modify the parameterization, such as enabling specific angular modes, changing the scalar-field potential $V(\phi)$, or altering the initialization of the shift vector $\beta^i$.
    \item \textbf{Automated Code Adaptation:} For complex physical setups, the agent writes example-local C++ headers (e.g., custom potentials or initial-data functors), runs pre-flight compilation checks, and resolves compiler warnings before launching new search episodes.
\end{enumerate}

\section{GPU Cluster Architecture and Program Boundary}

To efficiently search the parameter space, the pipeline exploits a clean boundary between Python and pre-compiled C++ executable code across a GPU cluster.

\begin{figure*}[t]
\centering
\includegraphics[width=0.95\textwidth]{magicbook-artwork.jpg}
\caption{Schematic diagram of the closed-loop classical spacetime optimization pipeline. Layer 1 handles initial spacetime metric data, boundary conditions, and spatial preprocessing. Layer 2 is the core iterative engine, pairing a numerical relativity solver (such as \texttt{GRTeclyn}) with the parameter optimization loop to check physical constraints (including causality and curvature). Layer 3 computes the objective function (evaluating constraint violations, energy conditions, and FTL travel times) and visualizes the resulting 4D embedding structures, directing the overall feedback and metric refinement loop.}
\label{fig:pipeline_diagram}
\end{figure*}

\subsection{The Execution Boundary}

The division of labor is designed to eliminate compile-time overhead:

\begin{itemize}
    \item \textbf{Python Controller (CPU/Head Node):} Runs the optimization algorithm, prepares a batch of parameter vectors $\theta$, writes them to individual text-based \texttt{params.txt} inputs, and dispatches them as parallel tasks to the cluster's Slurm scheduler.
    \item \textbf{C++ Initial Data (GPU/Compute Nodes):} At $t=0$, the compiled \texttt{GRTeclyn} binary starts up on the assigned GPU, reads the coefficients from \texttt{params.txt}, and executes a custom initial-data GPU functor to populate the 3D AMReX grid with smooth initial fields.
    \item \textbf{C++ RHS Evolution (GPU/Compute Nodes):} At $t > 0$, the GPU runs the pre-compiled, highly optimized CCZ4 and matter RHS evolution routines (e.g. \texttt{CCZ4RHSWithMatter}). Because the general relativity and matter equations are pre-compiled, the search loop runs with zero compilation overhead.
\end{itemize}

\subsection{Tiered GPU Execution}

To prevent a massive search from exhausting cluster node hours, we implement a tiered system:

\begin{itemize}
    \item \textbf{Tier 1 (Broad Exploration):} Each candidate is allocated to a single GPU and run at a low resolution (e.g. $N_{\text{full}} = 64$ with 2 levels of AMR) for a short duration. Early-stopping criteria instantly terminate a run if constraints explode, the lapse collapses everywhere, or the grid becomes unresolved.
    \item \textbf{Tier 2 (Deep Validation):} High-scoring candidates that survive Tier 1 are promoted to run on multi-GPU, multi-node configurations at high resolution ($N_{\text{full}} = 128 \to 256$, with up to 5 levels of AMR) to verify stability and convergence.
\end{itemize}

\section{Scoring and Energy Conditions}

The objective function must reward desirable physical outcomes and heavily penalize non-physical behavior or numerical errors. The total score $S(\theta)$ for a candidate parameter set is defined as:

\begin{equation}
    S(\theta) = w_1 F_{\text{FTL}} - w_2 \mathcal{H}_{\text{violation}} - w_3 E_{\text{exotic}} - w_4 T_{\text{tidal}},
\end{equation}

where:
\begin{itemize}
    \item $F_{\text{FTL}}$ is the measured FTL shortcut benefit (determined via the operational probe test).
    \item $\mathcal{H}_{\text{violation}}$ is the integrated volume norm of the Hamiltonian constraint violation. This ensures the optimizer is not rewarded for coordinate tricks that violate Einstein's equations.
    \item $E_{\text{exotic}}$ is the integrated volume of negative energy density ($T_{00} < 0$). By setting the weight $w_3$ to a very high value, the optimizer is strictly forced to seek configurations supported by normal, positive-energy matter.
    \item $T_{\text{tidal}}$ is the maximum tidal force felt by an observer in the passenger region.
\end{itemize}

\section{Avoiding Fake FTL (Operational Probe Tests)}

A major scientific risk in numerical relativity is "fake FTL," where a coordinate transformation (such as a massive shift vector $\beta^i$) makes it look like a region is moving superluminally when it is actually a coordinate effect. 

To eliminate coordinate artifacts, every high-scoring candidate must pass a physical operational test:

\begin{enumerate}
    \item An active probe (either a massless scalar-field pulse or a set of null geodesic tracer particles) is emitted at a boundary $A$.
    \item The pulse propagates through the dynamical spacetime to a detector at boundary $B$.
    \item We compare the arrival time at $B$ with an identical simulation run in flat Minkowski space on the exact same grid.
    \item If the signal arrives earlier, we verify that the signal path did not cross an apparent horizon (trapped surface), and that the constraint violations remained bounded during the entire transit time.
\end{enumerate}

\section{Benchmarks}

To ensure the pipeline is robust, it must first reproduce known physical configurations prior to starting the open-ended search.

\begin{table}[h]
\caption{Initial benchmark suite for the search pipeline.}
\label{tab:benchmarks}
\begin{ruledtabular}
\begin{tabular}{lll}
Case & Expected result & Purpose \\
\hline
Minkowski & Stable, no shortcut & Null/Flat test \\
Ellis--Bronnikov & Unstable throat & Wormhole detector \\
Perturbed EB & Collapse + GW signal & Previous-paper regression \\
Alcubierre-like & High exoticity & Fake-FTL coordinate test \\
Teo-like & Quadrupole dynamics & Non-spherical benchmark \\
\end{tabular}
\end{ruledtabular}
\end{table}

\section{Step-by-Step Implementation Plan}

\subsection{Step 1: Regression and Cleanup}
Establish a regression script that reproduces the unperturbed (expanding) and perturbed (collapsing) Ellis--Bronnikov wormhole. Configure the `grteclyn-wrapper` package to automatically extract the 1D diagnostics and programmatically delete the heavy 3D AMReX plotfiles (`Plt*`) to prevent storage exhaustion on the cluster.

\subsection{Step 2: Custom C++ Initial Data}
Implement the `Examples/RadialRecipe/` directory in C++. The initial data functor must read Gaussian basis coefficients directly from `params.txt` and map them onto the starting grid variables $\chi$, $\alpha$, $\beta^i$, $K$, $\phi$, and $\Pi$.

\subsection{Step 3: Slurm-Based Random Search (The Atlas)}
Write the Python scheduler to generate batches of randomized radial coefficients, write \texttt{params.txt} inputs, and dispatch them to the GPU cluster via Slurm. Parse the resulting \texttt{collapse\_diagnostics.dat} files and categorize the failure modes (black hole collapse, coordinate crash, or flat space dissipation) into an \texttt{atlas.csv} file.

\subsection{Step 4: Integrate the Constraint Solver}
Link \texttt{GRTresna} (or an equivalent elliptic constraint solver) into the pre-flight pipeline. Configure the system so that the optimizer's geometric target is deformed to be mathematically consistent with a positive-energy scalar field before starting the C++ time evolution.

\subsection{Step 5: Operational Probe Integration}
Implement the active null-ray or scalar pulse propagation test within the \texttt{specificPostTimeStep} diagnostic loop of \texttt{GRTeclyn} to measure physical arrival times and compute the final FTL score.

\subsection{Step 6: Optimization and Scale-up}
Replace the random search with a CMA-ES or Bayesian optimization loop. Once the radial optimization is stable, expand the parameter space to include angular modes (Level 2) and scale the pipeline across multiple GPU nodes.

\section{Expected First Result}

A successful first result is not a claimed working warp drive. The primary deliverable is an automated numerical-relativity classification laboratory. 

The first paper will present:
\begin{itemize}
    \item The verified, open-source closed-loop search pipeline.
    \item A comprehensive physical database ("The Spacetime Failure Atlas") showing how thousands of random geometries fail dynamically under Einstein's equations.
    \item The automated discovery of known metric classes (like wormholes) without manual human tuning.
    \item A rigorous evaluation of the physical viability of FTL metrics under the strict exclusion of negative energy density.
\end{itemize}

\section{Implemented: Constraint-Satisfying Metric Guesser}

The core technical challenge for the automated search is that arbitrary initial data almost always violates the Einstein constraint equations, causing immediate junk radiation, gauge failure, and solver crashes. We have implemented a six-component pipeline that bridges the gap between random metric proposals and physically meaningful evolutions.

\subsection{Physics of the Constraint Gap}

For the conformally flat recipe ($\tilde\gamma_{ij} = \delta_{ij}$, $A_{ij} = 0$), the Einstein constraints reduce to:
\begin{align}
    R[\chi] + \frac{2}{3}K^2 &= 16\pi\rho, \label{eq:ham_recipe} \\
    -\frac{2}{3}D^i K &= 8\pi j^i. \label{eq:mom_recipe}
\end{align}
For a scalar field with $\Pi = 0$ (momentarily static), the momentum density $j^i = -\Pi\,\partial^i\phi = 0$, so the momentum constraint requires $D^i K = 0$---that is, $K$ must be spatially constant. If the recipe generates spatially varying $K$ with $\Pi = 0$, the momentum constraint is violated. If $\chi$ and $\phi$ are independently random, the Hamiltonian constraint is violated.

\subsection{Component 1: Constrained Recipe (Python Pre-Processing)}

Rather than guessing $\phi$ independently, we derive it from $\chi$ so that the Hamiltonian constraint is approximately satisfied at $t=0$. The implementation (\texttt{constrained\_recipe.py}) performs the following steps:

\begin{enumerate}
    \item Lock $K$ to a spatial constant (use \texttt{recipe\_K\_asymptotic}, zero out all $K$ Gaussian coefficients). This exactly satisfies the momentum constraint when $\Pi = 0$.
    \item Lock $\Pi = 0$ (zero out all $\Pi$ coefficients).
    \item Evaluate $\chi(r)$ on a fine 1D radial grid from the proposed Gaussian basis coefficients.
    \item Compute the 3-Ricci scalar using $R = -8\psi^{-5}\nabla^2\psi$ where $\psi = \chi^{-1/4}$.
    \item Solve for $\phi(r)$ via radial integration of
    \begin{equation}
        \left(\frac{d\phi}{dr}\right)^2 = s\,\chi\,\frac{R + \frac{2}{3}K^2}{8\pi},
    \end{equation}
    where $s = +1$ for a normal scalar field and $s = -1$ for a phantom field. Regions where the integrand is negative are flagged as unsatisfiable.
    \item Fit the resulting $\phi(r)$ back to Gaussian basis coefficients via least-squares, and inject them into \texttt{params.txt}.
\end{enumerate}

This is activated by the \texttt{--constrained} flag in the CLI. The optimizer then searches only over $\chi$ coefficients and basis parameters; $\phi$ is always derived consistently.

\subsection{Component 2: Pre-Flight Constraint Filter}

Before launching \texttt{GRTeclyn} on the GPU, the wrapper evaluates the Hamiltonian and momentum constraint residuals on a coarse 1D radial grid (\texttt{preflight.py}). Candidates with $\|H\|_{L^2}$ or $\|M_i\|_{L^2}$ exceeding configurable thresholds, or with $\chi < 0$ in more than 10\% of the grid, are rejected without spending GPU time. This is activated by the \texttt{--preflight} flag.

\subsection{Component 3: Energy Density Diagnostic (C++)}

We extended \texttt{RadialRecipeLevel.cpp} to compute the vacuum Hamiltonian constraint (without matter subtraction) at every grid point, giving the \emph{required} energy density:
\begin{equation}
    \rho_{\rm required}(\mathbf{x}) = \frac{1}{16\pi}\left(R + \frac{2}{3}K^2 - A_{ij}A^{ij}\right).
\end{equation}
Three new columns are written to \texttt{constraint\_norms.dat} at every diagnostic step: $\min(\rho_{\rm req})$, $\max(\rho_{\rm req})$, and the integrated volume of negative $\rho_{\rm req}$. These allow the Python scorer to evaluate energy-condition compliance without a separate post-processing pass.

\subsection{Component 4: Enhanced Scoring}

The episode scorer (\texttt{score.py}) now includes two additional weighted components beyond the original survival, constraint health, lapse health, horizon penalty, and nontrivial geometry:

\begin{itemize}
    \item \textbf{Energy condition} (weight 1.5): rewards configurations where $\rho_{\rm req} \geq 0$ everywhere; penalizes the integrated volume of negative energy density.
    \item \textbf{Initial constraint quality} (weight 1.0): bounded reward based on the initial Hamiltonian $L^2$ norm, distinguishing clean initial data from noisy starts.
\end{itemize}

All weights are configurable via a dictionary parameter, enabling task-specific objective functions (e.g., maximizing FTL shortcut benefit vs.\ minimizing exotic matter).

\subsection{Component 5: CMA-ES Optimization Driver}

The random atlas search has been supplemented with a CMA-ES (Covariance Matrix Adaptation Evolution Strategy) optimizer (\texttt{optimize.py}). The default search space consists of the four $\chi$ Gaussian basis coefficients and the basis width (5 dimensions). Each evaluation runs a full \texttt{GRTeclyn} episode and returns the negative total score. The optimizer is invoked via:
\begin{verbatim}
python -m grteclyn_wrapper --example RadialRecipe \
    --constrained --preflight optimize \
    --max-generations 50 --seed 42
\end{verbatim}
CMA-ES automatically adjusts the search distribution based on the fitness landscape, concentrating evaluations in promising regions of the coefficient space.

\subsection{Component 6: Known-Solution Seeds}

Three analytically known spacetimes are provided as regression benchmarks and optimizer warm-start points (\texttt{seeds.py}):

\begin{enumerate}
    \item \textbf{Flat Minkowski}: $\chi = 1$, all coefficients zero. All constraints are exactly satisfied; the evolution should remain static to machine precision.
    \item \textbf{Ellis--Bronnikov wormhole}: The exact isotropic $\chi$ and $\phi$ profiles (Eqs.~7 and 11) are fitted to the Gaussian basis via least-squares. This verifies that the recipe basis can faithfully represent a known exotic geometry.
    \item \textbf{Schwarzschild puncture}: The Brill--Lindquist conformal factor $\psi = 1 + M/(2\bar{r})$ fitted to Gaussians. A vacuum solution that tests the recipe's ability to handle steep $\chi$ gradients near a singularity.
\end{enumerate}

These seeds are selected via \texttt{--seed-name} and can be combined with the optimizer to search for perturbations around known solutions.

\subsection{Verification Tests}

All six components were verified via an automated Python test suite executed through \texttt{uv run}. The following checks confirm correctness of the pipeline before any GPU time is spent:

\begin{enumerate}
    \item \textbf{Constrained recipe self-consistency}: The Ellis--Bronnikov seed's $\chi$ profile is passed through \texttt{constrained\_overrides} with \texttt{phantom=True}. The derived $\phi$ coefficients match the analytical values to machine precision ($\max|\phi_{\rm orig} - \phi_{\rm derived}| = 0$), confirming the Hamiltonian constraint solver correctly recovers the exact phantom scalar field.
    \item \textbf{Trivial solution}: Flat Minkowski ($\chi = 1$, all coefficients zero) passed through the constrained recipe produces identically zero $\phi$ coefficients, verifying no spurious matter content is introduced.
    \item \textbf{Pre-flight acceptance}: The Ellis--Bronnikov seed passes the pre-flight filter with $\|H\|_{L^2} \sim 10^{-1}$, consistent with a near-constraint-satisfying configuration.
    \item \textbf{Pre-flight rejection}: A deliberately pathological configuration ($\chi_0 = -5$) is correctly rejected: $\|H\|_{L^2} \approx 2.7 \times 10^3 \gg 10$ and $\chi < 0$ over $\sim20\%$ of the grid, triggering three independent rejection criteria.
    \item \textbf{Seed catalogue}: All three seeds (flat Minkowski, Ellis--Bronnikov, Schwarzschild puncture) load and produce valid override dictionaries with physically correct coefficient structure.
    \item \textbf{Optimizer readiness}: The CMA-ES driver imports successfully with \texttt{cma 4.4.4}, the 5-dimensional default search space is well-formed, and all scoring weights (including the new energy condition and initial constraint quality components) are registered.
\end{enumerate}

\section{Conclusion}

This project moves the search for interstellar metrics away from static, hand-written equations toward dynamic, automated discovery. By allowing the computer to propose complex, non-spherical initial data, and using the full power of 3D numerical relativity on a GPU cluster to judge the results, we can rigorously map the physical boundaries of superluminal transport and communication.

\end{document}