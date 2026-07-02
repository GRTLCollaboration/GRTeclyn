### 3. Bicomplex Magnetic Shielding (The Matter/Anti-Matter Core)
**The Observation:** The lumps dispersed to 23% confinement because centripetal forces ripped them apart. But we have the `bicomplex_scalar` model (Normal Matter $\Phi^+$ and Exotic Matter $\Phi^-$).
**The New Idea:** Stop making the exotic matter try to hold itself together. Use the normal matter as a physical containment shell.
*   **How it works:** In `spaces.py`, force the AI to build concentric pairs. Place a ring of normal matter at radius $R=5.0$, and a ring of exotic matter at $R=4.8$. 
*   **The Physics:** Because normal matter and exotic matter have opposite gravitational signs, they will push/pull against each other. If the AI spins them at the right speeds, the outward centripetal force of the exotic matter will be perfectly blocked by the inward gravitational pressure of the normal matter. You essentially build an invisible "magnetic bottle" out of gravity, allowing you to spin the exotic matter at near-light speed without it dispersing.


### 4. The "Spark Plug" (Multi-Pulse Injection)
**The Observation:** The FTL channel reliably lasts about 16 to 20 code units before the matter dissolves. 
**The New Idea:** Car engines don't ignite the fuel once; they pulse the spark plugs thousands of times a minute. We need to do the same to the spacetime grid.
*   **How it works:** Modify the `TrajectoryEvaluator` pump to turn on and off. Instead of just placing the Q-balls and watching them die, add a `pulse_interval` dimension to your search space.
*   **The Physics:** The engine injects Q-ball matter at $t=0$. The FTL wave peaks at $t=15$. At exactly $t=15$, the engine injects a *second* pulse of matter into the moving wake, re-energizing the metric. 
*   **Implementation:** If the AI can optimize the exact timing of these injections, you turn a transient 20-second FTL fade into a permanent, steady-state FTL corridor. 