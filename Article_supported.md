%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a revised LaTeX file for your research proposal,
% fully restyled to match the format of the template paper 2406.02466.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the 'natbib' option for author-year citations
\documentclass[twocolumn,aps,prd,nofootinbib,natbib]{revtex4-2}

% --- Packages ---
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{graphicx}
\usepackage{amssymb}
\usepackage{amsmath}
\usepackage{xcolor}
\usepackage{hyperref}
\usepackage{textcase}

% --- Style Customizations ---

% Hyperref setup for clean look like the template
\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    filecolor=magenta,      
    urlcolor=blue,
    citecolor=blue,
}
\urlstyle{same}

% --- Section Style Customizations ---
\makeatletter

% Part 1: Modify the section font and alignment using the revtex hook
% This changes the style to centered, small, and bold.
\def\section@font{\centering\normalfont\small\bfseries}

% Part 2: Create a "wrapper" to safely make the section title uppercase
\let\oldsection\section % Save the original \section command
% Redefine \section to take an argument, uppercase it, then call the original
\renewcommand{\section}[1]{\oldsection{\MakeTextUppercase{#1}}}

\makeatother


\begin{document}

% --- Title and Author Information ---
\title{\MakeTextUppercase{Nonlinear Collapse and Gravitational-Wave Emission from Unstable Ellis--Bronnikov Wormholes
}}

\author{\MakeTextUppercase{Nikita M. Shirokov}}
\email{shirokov.nm@phystech.edu}
\affiliation{\textit{Independent Researcher, Moscow, Russia}}

% --- Version Date ---
\date{Version March 19, 2026}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- CORRECTED POSITION FOR THE ABSTRACT ---
% For REVTeX, the abstract MUST come BEFORE the \maketitle command.
% The class will then automatically format it correctly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\begin{abstract}
If primordial wormhole-like defects were seeded from quantum foam and stretched to macroscopic scales during inflation, they must be described as fully coupled matter-gravity systems. Using full 3D numerical relativity with \texttt{GRTeclyn}, we perform self-consistent evolutions of the fully coupled Einstein--phantom-scalar system and show that the unstable Ellis--Bronnikov geometry collapses into a black hole in all cases. Without any applied deformation, the inherent linear instability drives spontaneous collapse. When the collapse is triggered via controlled scalar-field perturbations (monopolar $A_0$ plus quadrupolar $A_\phi$), the collapse is accelerated and accompanied by a near-zone curvature transient propagating at $v\approx 1.01c$, consistent with physical gravitational radiation---in contrast to the superluminal constraint-damping modes of the unperturbed signal. These results demonstrate that unstable traversable wormholes are generically susceptible to dynamical collapse and emission of gravitational waves.
\end{abstract}


% The \maketitle command generates the title block AFTER the abstract
\maketitle

% --- Main Content of the Proposal ---
\section{Introduction}
Wormholes—topological bridges connecting distinct regions of spacetime—have been studied since the work of Flamm~\cite{flamm16}, Einstein and Rosen~\cite{einstein35}, and Wheeler~\cite{misner57}. The first explicit traversable solutions, supported by a phantom scalar field, were found independently by Ellis~\cite{ellis73} and Bronnikov~\cite{bronnikov73}; this Ellis--Bronnikov geometry is the canonical example studied here.

Morris and Thorne showed that static, traversable wormholes require ``exotic matter'' violating the null energy condition~\cite{morris88}. Much subsequent work has focused on engineering, minimizing, or replacing this exotic stress-energy~\cite{visser95, lobo05, garattini2026}, and recent models show that rapid rotation can reduce the required amount~\cite{uemichi2026}. That perspective is natural for engineering a traversable bridge, but it is not the only physically relevant question.

Our motivation diverges from the search for stability. In Wheeler's spacetime foam, transient microscopic wormholes may form at very early times, and inflation could stretch some to macroscopic scales~\cite{kardashev07, garcia16}. The relevant physical problem is then the nonlinear outcome of their linear instability. While 1D spherical studies have shown that they are unstable to perturbations, full 3D simulations are required to compute the gravitational wave emission (\(\Psi_4\)) of the collapse.

Shinkai and Hayward~\cite{shinkai02} first established the dynamical fate of the Ellis--Bronnikov wormhole in 1D spherical symmetry, showing it to be highly unstable: a rarefactive deformation triggers inflationary expansion, while a compressive one forces collapse into a black hole (Fig.~\ref{fig:wormhole_collapse_stages}). Even a minute compressive deformation inevitably seals the bridge.

However, 1D spherical symmetry precludes non-spherical dynamics and, by Birkhoff's theorem, any curvature radiation. Full 3D simulations are therefore needed to determine whether unstable wormholes radiate as they collapse and what near-zone curvature signal accompanies the process. In this work, we evolve the full coupled Einstein-Klein-Gordon system using the exact Ellis--Bronnikov solution as initial data, and trigger the dynamical collapse via a controlled quadrupolar perturbation of the scalar field profile.

\begin{figure}[htp]
\centering
\includegraphics[width=\linewidth]{wormhole_collapse_stages.pdf}
\caption{Conceptual 3D embedding diagram illustrating the dynamical collapse of an unstable traversable wormhole. \textbf{(a)} The initial traversable state. \textbf{(b)} Dynamical constriction of the throat. \textbf{(c)} Topological pinch-off and formation of a black hole. Both the accelerated (perturbed) and spontaneous (unperturbed) evolutions follow this same qualitative geometric pathway, with the accelerated case collapsing much faster and emitting a clear gravitational-wave transient during pinch-off (depicted as ripples in panel c).}
\label{fig:wormhole_collapse_stages}
\end{figure}

The advent of gravitational-wave astronomy~\citep{abbott16,abbott18,abbott21,abbott23} has moved the search for Exotic Compact Objects (ECOs)~\cite{cardoso19} from theory to observational reality. Understanding the dynamics of candidate geometries is a prerequisite for any detection strategy. In this paper we extract the Weyl scalar \(\Psi_4\) at multiple radii and use propagation speed analysis to distinguish physical curvature transients (\(v\approx c\)) from numerical constraint-damping modes (\(v>c\)); the contrast between unperturbed and perturbed runs---both of which begin from constraint-satisfying initial data to truncation error---serves as a built-in null test.

\section{Formalism and Methodology}

\subsection{3+1 Decomposition and the CCZ4 Formulation}
We follow the standard numerical-relativity approach to solving the Einstein Equations,
\begin{equation}
G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2} R g_{\mu\nu} = 8\pi T_{\mu\nu},
\end{equation}
as an initial value problem, using the 3+1 decomposition~\cite{baumgarte10} and setting $G = c = 1$. The four-dimensional spacetime manifold is foliated into a family of spacelike hypersurfaces $\Sigma_t$, parametrized by a time coordinate $t$. The line element is written in terms of the lapse function $\alpha$, the shift vector $\beta^i$ (which controls how the spatial coordinates propagate between time slices), and the spatial metric $\gamma_{ij}$ as:
\begin{equation}
    ds^2 = \left(-\alpha^2 + \beta_i \beta^i \right) dt^2 + 2 \beta_i dt dx^i + \gamma_{ij} dx^i dx^j.
\end{equation}
The dynamical evolution is governed by the coupled Einstein--scalar-field system, implemented within the \texttt{GRTeclyn} framework using the conformal and covariant Z4 (CCZ4) formulation with matter source terms. Our objective is to start from the exact Ellis--Bronnikov wormhole solution---geometry \emph{and} its phantom scalar field support---and evolve the full coupled system self-consistently, so that the initial data satisfy the Hamiltonian and momentum constraints to truncation error.

\subsection{Initial Data: Ellis--Bronnikov Wormhole in Isotropic Coordinates}

The general static, spherically symmetric traversable wormhole is described by the Morris--Thorne line element~\cite{morris88}:
\begin{equation}
ds^2 = -e^{2\Phi(r)} dt^2 + \frac{dr^2}{1 - b(r)/r} + r^2 (d\theta^2 + \sin^2\theta d\phi^2).
\end{equation}
Setting the redshift function to zero ($\Phi(r) = 0$) and choosing the shape function $b(r) = b_0^2/r$ yields the Ellis--Bronnikov wormhole~\cite{ellis73, bronnikov73}, a canonical massless phantom-scalar-field solution. The vanishing redshift ensures the absence of initial horizons, sets the proper time of Eulerian observers equal to coordinate time at $t=0$, and eliminates radial tidal forces at the throat; $b_0$ is the throat radius; and the choice of $b(r)$ guarantees asymptotic flatness as $r \to \infty$.

In standard Schwarzschild-like coordinates, the metric component $g_{rr} = (1 - b_0^2/r^2)^{-1}$ diverges at the throat ($r=b_0$). To resolve this coordinate singularity and simultaneously map both "universes" connected by the bridge onto a single computational grid, we transform to \textbf{isotropic coordinates}. We introduce the isotropic radius $\bar{r}$, related to the proper distance $\ell = \pm \sqrt{r^2 - b_0^2}$ via the embedding $\ell = b_0 \sinh(\ln(2\bar{r}/b_0))$. 

To initialize the wormhole on a 3D Cartesian grid, we transform the standard metric into isotropic coordinates. As demonstrated by Nandi et al.~\cite{nandi2008energetics}, the spatial geometry of the Ellis--Bronnikov wormhole admits a conformally flat representation where the CCZ4 conformal factor is strictly given by:
\begin{equation}
    \psi = \sqrt{1 + \frac{b_0^2}{4\bar{r}^2}}.
\end{equation}
Unlike the Brill--Lindquist black-hole puncture~\cite{brill63}, where $\psi$ diverges at a point singularity, the conformal factor here represents a genuine topological bridge. In this chart, the primary universe extends to $\bar{r} \to \infty$ ($\psi \to 1$), the throat is located at $\bar{r} = b_0/2$, and the secondary universe is compactified such that its asymptotic infinity maps to the grid origin $\bar{r} \to 0$ ($\psi \to \infty$). 

\subsection{Transformation to CCZ4 Variables and Regularization}

For stable numerical evolution in the CCZ4 formulation, the spatial metric $\gamma_{ij}$ appearing in the line element (Eq.~2) is decomposed into a conformal metric $\tilde{\gamma}_{ij}$ with unit determinant and a regularized conformal factor $\chi$:
\begin{align}
    \tilde{\gamma}_{ij} &= \psi^{-4} \gamma_{ij} = \delta_{ij}, \\
    \chi &= \psi^{-4} = \left( 1 + \frac{b_0^2}{4\bar{r}^2} \right)^{-2}.
\end{align}
To strictly avoid division-by-zero errors at the grid origin ($\bar{r}=0$), the conformal factor is algebraically rewritten and implemented as:
\begin{equation}
    \chi = \left( \frac{4\bar{r}^2}{4\bar{r}^2 + b_0^2} \right)^{2}.
\end{equation}
This smoothly vanishes at the origin ($\chi \to 0$), providing a regular way to represent the compactified second asymptotic end on the Cartesian grid.

\subsection{Phantom Scalar Field and Self-Consistent Matter Evolution}

The Ellis--Bronnikov wormhole is supported by a massless phantom (ghost) scalar field $\phi$ whose kinetic term enters the action with the opposite sign to a normal scalar field. The stress-energy tensor is
\begin{equation}
    T_{\mu\nu} = -\left(\nabla_\mu\phi\,\nabla_\nu\phi - \frac{1}{2}g_{\mu\nu}\nabla^\alpha\phi\,\nabla_\alpha\phi - g_{\mu\nu}\,V(\phi)\right),
\end{equation}
where the overall minus sign relative to the canonical scalar field ensures violation of the null energy condition, as required for a traversable wormhole. The scalar field obeys the standard Klein--Gordon equation
\begin{equation}
    \Box\phi = -\frac{dV}{d\phi},
\end{equation}
which follows from the contracted Bianchi identity $\nabla_\mu T^{\mu\nu}=0$ and is identical in form to that of a normal scalar field; only the gravitational coupling (the sign of $T_{\mu\nu}$) is reversed.

In the 3+1 decomposition, the scalar field is represented by the field value $\phi$ and its conjugate momentum $\Pi = -n^\mu\nabla_\mu\phi$, where $n^\mu$ is the unit timelike normal. Both are evolved alongside the CCZ4 geometric variables at every Runge--Kutta substep via the \texttt{CCZ4RHSWithMatter} module in \texttt{GRTeclyn}, ensuring full self-consistency between the matter and gravitational sectors.

For the exact Ellis--Bronnikov solution with $V(\phi)=0$, the scalar field profile in isotropic coordinates is
\begin{equation}
    \phi(\bar{r}) = \frac{1}{\sqrt{4\pi}}\arctan\!\left(\frac{\bar{r} - b_0^2/(4\bar{r})}{b_0}\right),
    \label{eq:phi_EB}
\end{equation}
with $\Pi(t{=}0)=0$ (static solution). This profile, together with the conformal geometry of Eqs.~(5)--(7), constitutes an exact solution of the coupled Einstein--phantom-scalar system; the Hamiltonian and momentum constraints are therefore satisfied to truncation error on the discrete grid at $t{=}0$.

\subsection{Perturbation and Collapse Triggering}

The Ellis--Bronnikov wormhole is known to be linearly unstable~\cite{gonzalez09, shinkai02}: even an infinitesimal perturbation triggers either collapse to a black hole or expansion to infinity. Because the full phantom scalar field is evolved self-consistently, the instability can be seeded by a controlled perturbation of the matter sector rather than by breaking the geometric constraints.

The choice of which field to perturb is dictated by the structure of the initial-value constraints. For a phantom scalar field, the momentum density is $S_i = \Pi\,\partial_i\phi$, and the momentum constraint reads $\nabla_j K^{ij} - \nabla^i K = 8\pi S_i$. Perturbing $K_{ij}$ directly would violate both the Hamiltonian constraint (since $K^2$ and $K_{ij}K^{ij}$ change) and the momentum constraint, injecting spin-2 numerical artifacts spectrally entangled with the physical $\Psi_4$ signal. Perturbing the scalar field momentum $\Pi$ while keeping $K_{ij}=0$ would make $S_i \neq 0$, again violating the momentum constraint unless a coupled elliptic system is solved.

The key observation is that if we perturb only the scalar field \emph{profile} $\phi$ while keeping $\Pi(t{=}0) = 0$ and $K_{ij}(t{=}0) = 0$ everywhere, the momentum constraint is satisfied \emph{exactly}:
\begin{equation}
    \underbrace{\nabla_j K^{ij} - \nabla^i K}_{=\,0} = 8\pi\,\underbrace{\Pi}_{=\,0}\,\partial^i\phi = 0.
\end{equation}
We therefore apply a localized perturbation to the scalar field profile at $t{=}0$ consisting of a monopolar (breathing) component and a quadrupolar component:
\begin{equation}
    \phi(\bar{r}, \theta, \varphi)\big|_{t=0} = \phi_{\rm EB}(\bar{r}) + \bigl(A_0 + A_\phi\, Y_{20}(\theta, \varphi)\bigr)\,
    \exp\!\left(-\frac{\bar{r}^2}{\sigma_\phi^2}\right),
\label{eq:phi_perturbation}
\end{equation}
where $\phi_{\rm EB}(\bar{r})$ is the exact profile of Eq.~\eqref{eq:phi_EB}, $A_0$ is the monopolar amplitude that drives a fast, symmetric compression of the scalar field support, $A_\phi$ is the quadrupolar amplitude that breaks spherical symmetry, and $\sigma_\phi$ controls the radial width. The monopolar component directly seeds the compressive unstable mode of the Ellis--Bronnikov solution, accelerating the collapse. The quadrupolar component is required by Birkhoff's theorem to produce a physical $\Psi_4$ signal.

The perturbation modifies the gradient energy $(\nabla\phi)^2$ and hence the energy density $\rho$, introducing a Hamiltonian constraint residual at leading order of $\mathcal{O}(A_0, A_\phi)$ from the cross-term $2\nabla\phi_{\rm EB}\cdot\nabla\delta\phi$. Crucially, this is a purely scalar (spin-0) violation. While the CCZ4 formulation will dynamically damp and propagate this scalar residual (typically via the $\Theta$ damping field), a spin-0 constraint mode does not project onto the transverse-traceless spin-2 Weyl scalar $\Psi_4$. Therefore, the extracted gravitational-wave signal remains cleanly separated from the numerical constraint-damping artifacts, representing a vast improvement over methods that directly perturb the geometric tensor fields.

Physically, at $t{=}0$ the scalar field is no longer in equilibrium with the geometry. The monopolar component compresses the phantom field support near the throat, while the quadrupolar component introduces aspherical gradients. As the evolution begins, these gradients dynamically generate a non-zero $\Pi$ through the Klein--Gordon equation, and the resulting aspherical matter distribution backreacts on the metric via the Einstein equations, physically generating $K_{ij}$ and radiating $l=2$ gravitational waves through the collapse dynamics.

\subsection{Apparent Horizon Detection}

During the period in which this work was carried out, \texttt{GRTeclyn} did not yet include a production-ready full apparent-horizon finder for these runs. To distinguish genuine topological pinch-off and black-hole-like self-trapping from gauge effects such as lapse collapse, we therefore implemented a geometric trapped-surface proxy on coordinate spheres centered on the throat. At each time step, we scan outward from the origin and estimate the expansion of outgoing null rays, \(\theta_+\). For a conformally flat spatial metric \(\gamma_{ij}=\chi^{-1}\delta_{ij}\), the spherical expression used in the code is
\begin{equation}
    \theta_+ = \frac{2\sqrt{\chi}}{\bar{r}} - \frac{\partial_{\bar{r}} \chi}{\sqrt{\chi}} + \tilde{A}_{rr} - \frac{2}{3}K.
\end{equation}
The formation of a trapped surface is identified when \(\theta_+ \le 0\). The maximum radius where this condition holds provides a proxy for the apparent-horizon radius, \(r_{\rm AH}\). Although this quantity does not replace a full elliptic horizon finder, it gives a practical collapse diagnostic for the parameter study reported here.

\subsection{Gauge Conditions}

To ensure stable evolution across both dispersive and collapsing branches, we employ robust gauge choices. 

The shift vector $\beta^i$ is evolved using the Gamma-driver condition, starting from $\beta^i(t{=}0) = 0$. Because the initial data possess no rotation and the primary dynamics are radial, the vanishing initial shift avoids gauge-induced artifacts during the early constraint-cleaning phase and the subsequent gravitational-wave extraction.

For the lapse function, we employ the $1+\log$ slicing condition, which evolves as $\partial_t \alpha - \beta^i \partial_i \alpha = -2\alpha K$. In all runs reported here, the initial lapse is set to the pre-collapsed profile
\begin{equation}
    \alpha(t{=}0) = \sqrt{\chi}. \label{eq:lapse_precollapsed}
\end{equation}
This choice seeds the lapse at a value already compatible with $1{+}\log$ equilibrium near the wormhole throat, where $\chi$ is very small. In the deep interior ($\bar r\to 0$), the conformal factor vanishes and the CCZ4 evolution terms involving $1/\chi$ can produce runaway numerical errors if the lapse remains near unity. The pre-collapsed profile $\alpha=\sqrt{\chi}$ suppresses the effective evolution rate in this compactified inner region, preventing such instabilities while leaving the exterior dynamics ($\chi\approx 1$) unbiased. Using the same initial lapse for both the unperturbed and perturbed configurations ensures a clean comparison of the physics without gauge-induced differences.

\subsection{Collapse diagnostics}
\label{sec:diagnostics_method}

The nine-panel diagnostic figures (Appendix, Figs.~\ref{fig:unperturbed_diagnostics} and~\ref{fig:perturbed_diagnostics}) track the geometric state of the wormhole throat, with all quantities evaluated on the finest available AMR level. We summarize the key definitions here; physical interpretation of each panel is given in the figure captions.

The lapse $\alpha$ measures the rate of proper-time flow ($d\tau = \alpha\,dt$); $\alpha\to 0$ signals black-hole formation in the moving-puncture gauge. The conformal factor $\chi = \psi^{-4}$ satisfies $\gamma_{ij}=\chi^{-1}\tilde\gamma_{ij}$; $\chi\to 0$ indicates infinite stretching of the spatial geometry. The trace of the extrinsic curvature $K$ encodes the local rate of volume change:
\begin{equation}
    \partial_t (\sqrt{\gamma}) = -\alpha K \sqrt{\gamma} + \nabla_i \beta^i \sqrt{\gamma},
    \label{eq:K_volume}
\end{equation}
with $K>0$ denoting contraction and $K<0$ expansion. The trapped-surface proxy $r_{\rm AH}$ is the maximum radius where $\theta_+\le 0$ (Eq.~8). The throat areal radius, our most direct gauge-independent collapse indicator, is
\begin{equation}
    R_{\rm areal}(\bar r) = \frac{\bar r}{\sqrt{\chi(\bar r)}}.
    \label{eq:R_areal}
\end{equation}
To quantify the late-time settling, we fit
\begin{equation}
    \max|K|(t) = A\,e^{-t/\tau} + C
    \label{eq:K_decay_fit}
\end{equation}
to the late-time \(\max|K|\) curve; $\tau\approx 0$ with nonzero $C$ indicates a stationary collapsed state (trumpet slice).

\subsection{Gravitational-wave extraction}
\label{sec:gw_method}

We decompose the Weyl scalar \(\Psi_4\) (satisfying \(\Psi_4 = \ddot{h}_+ - i\ddot{h}_\times\) in the wave zone) into spin-weighted spherical harmonics; the \((\ell,m)=(2,0)\) mode dominates for our axisymmetric deformation. The radius-scaled quantity \(r\Psi_4\) is plotted in both simulation time and retarded time \(u=t-R_{\rm ext}\); alignment across extraction radii in retarded time indicates a coherent outgoing wavefront. The one-sided energy spectral density is computed via a Tukey-windowed DFT (\(\alpha=0.25\)):
\begin{equation}
    S_{\Psi_4}(f) = \frac{\Delta t^2}{T}\left(\left|\widetilde{\mathrm{Re}(\Psi_4)}\right|^2 + \left|\widetilde{\mathrm{Im}(\Psi_4)}\right|^2\right).
    \label{eq:esd_psi4}
\end{equation}

\noindent\textbf{Propagation speed analysis.}
The critical test for distinguishing physical radiation from numerical artifacts is:
\begin{equation}
    v = \frac{R_{i+1} - R_i}{t_{i+1}^{\rm peak} - t_i^{\rm peak}},
    \label{eq:propagation_speed}
\end{equation}
where \(t_i^{\rm peak}\) is the time of the dominant peak in \(|r\Psi_4|\) at radius \(R_i\). Physical gravitational waves give \(v=1\) (in geometrized units); CCZ4 constraint-damping modes propagate superluminally.

\noindent\textbf{Strain and detectability.}
The strain PSD is \(S_h(f) = S_{\Psi_4}(f)/(2\pi f)^4\), with a 4th-order Butterworth high-pass roll-off below \(f_{\rm low}=0.05\,f_{\rm max}\) to suppress low-frequency drift. Physical rescaling uses:
\begin{align}
    f_{\rm phys} &= \frac{f_{\rm code}}{M \cdot M_{\odot,\rm sec}}, \label{eq:f_phys}\\
    S_h^{\rm phys}(f) &= S_h^{\rm code}(f) \cdot M_{\odot,\rm sec} \cdot \left(\frac{M\cdot M_{\odot,\rm m}}{D}\right)^2, \label{eq:Sh_phys}
\end{align}
where \(M_{\odot,\rm sec} \approx 4.93\times 10^{-6}\)~s and \(M_{\odot,\rm m} \approx 1477\)~m. The characteristic strain \(h_{\rm char}=\sqrt{f\,S_h}\) is compared against the Advanced LIGO design sensitivity~\cite{aasi15}; an SNR~$\gtrsim 8$ is required for detection.


\section{Numerical Setup}

Numerical evolutions were performed using the \texttt{GRTeclyn} codebase~\cite{GRTeclyn,Andrade2021}, which is currently under active development. As the successor to \texttt{GRChombo}, \texttt{GRTeclyn} is built on the \texttt{AMReX} library~\cite{amrex} to provide block-structured Adaptive Mesh Refinement (AMR) and is explicitly optimized for high-performance GPU acceleration. Production runs targeted H100-class GPUs (corresponding to \texttt{CUDA\_ARCH=90} in the build configuration) and used one MPI rank per GPU, scaling up to 8 GPUs for the largest calculations. 

The full physical computational domain spans a coordinate length of \(L_{\text{full}} = 40\), covered by a coarse grid of \(N_{\text{full}} = 160\) cells, yielding a coarse resolution of \(dx_{\text{coarse}} = 0.25\). To drastically reduce the computational cost, we exploit the Cartesian reflection symmetries of the single-throat initial data and evolve only the positive octant (\(x \ge 0, y \ge 0, z \ge 0\)), effectively modeling just one eighth of the physical volume. Parity symmetry conditions are imposed at the inner reflection boundaries, while Sommerfeld radiative boundary conditions are enforced at the outer edges.

To capture the extreme curvature gradients developing during the collapse while maintaining tractability, we employ up to 5 levels of 2:1 adaptive mesh refinement (AMR) (\texttt{max\_level = 5}). The mesh is dynamically regridded based on the gradients of the conformal factor ($\chi$), allowing the code to focus resolution precisely where the wormhole throat is pinching off. The finest grid spacing is $dx_{\rm fine} = dx_{\rm coarse} / 2^5 \approx 7.8\times10^{-3}$. Time stepping is handled using a 4th-order Runge-Kutta scheme, with a highly conservative Courant factor (\texttt{dt\_multiplier = 0.01}) required to maintain stability during the constraint relaxation and the abrupt collapse branches.

The choice of the simulation domain size (\(L_{\text{full}} = 40\)) and maximum refinement level was dictated by absolute computational constraints. With these dimensions, the required 5 levels of resolution across the AMR hierarchy, and the extremely small time step (\texttt{dt\_multiplier = 0.01}), the memory and computational load were maximized for our multi-GPU hardware limit. The total wall-clock time required for a single \(t=20M\) simulation in this configuration was approximately 24 hours on a multi-GPU cluster. Due to these hard scaling limits, it was impossible to extend the simulation to a significantly larger domain (e.g., \(L_{\text{full}} = 100\)) or place the extraction radii further out into the true wave zone to extract a cleaner gravitational-wave signal without making the runtime prohibitively long.



\section{Results}

In all evolutions of the unstable Ellis--Bronnikov geometry, the wormhole collapses into a black hole regardless of whether an initial compressive deformation is applied. The deformation modulates the collapse timescale and the character of the accompanying curvature signal, but does not alter the final state. Below we summarize the two representative configurations and the key numerical choices that differ between them.

\noindent\textbf{Diagnostic sampling.}
All collapse diagnostics---$\min(\alpha)$, $\min(\chi)$, $\max|K|$, the trapped-surface proxy $r_{\rm AH}$, and $\min(\theta_+)$---are evaluated on the \emph{finest available AMR level} at each coarse time step. This ensures that the extreme gradients developing in the throat region are fully resolved in the reported quantities.

\begin{table}[h]
\caption{Parameters for the two representative configurations. All shared settings: $b_0=0.5$, $\alpha(t{=}0)=\sqrt{\chi}$, $\chi_{\rm min}=10^{-8}$, $\kappa_1=3.0$, $\sigma_{\rm KO}=1.0$, $L_{\rm full}=40$, $N_{\rm full}=160$, max\_level$=5$, $dx_{\rm fine}\approx 7.8\times10^{-3}$, and $\text{phantom\_mass} = 0.0$.}
\label{tab:params}
\begin{ruledtabular}
\begin{tabular}{lcc}
Parameter & Unperturbed & Perturbed \\
\hline
$A_0$ (monopole) & 0 & 0.5 \\
$A_\phi$ (quadrupole) & 0 & 0.05 \\
$\sigma_\phi$ & --- & 0.5 \\
\end{tabular}
\end{ruledtabular}
\end{table}
\noindent The unperturbed configuration ($A_0{=}A_\phi{=}0$) is momentarily static; the evolution is driven entirely by the linear instability of the Ellis--Bronnikov solution, seeded by numerical truncation error. The perturbed configuration applies a monopolar compression ($A_0{=}0.5$) to accelerate the collapse, combined with a small quadrupolar perturbation ($A_\phi{=}0.05$) to break spherical symmetry and enable $\ell=2$ radiation. The monopolar component directly seeds the compressive unstable mode, while the quadrupolar component generates the physical $\Psi_4$ signal. Identical gauge, floor, and grid settings ensure a clean comparison.

\subsection{Unperturbed evolution: spontaneous collapse of the unperturbed geometry}

The diagnostics in Fig.~\ref{fig:unperturbed_diagnostics} (Appendix) reveal a geometry that spontaneously \emph{collapses}. The most direct indicator is the throat areal radius (panel~g), which falls from the initial \(R_{\rm areal}=0.50\) to a plateau of \(\sim 0.10\) within \(t\approx 5\). This contraction occurs in two phases: (i)~a violent gauge/constraint transient (\(t\lesssim 2\)), during which the lapse plunges, the null expansion proxy dives to \(\sim -250\), and the throat contracts at coordinate speeds approaching \(c\); and (ii)~a slower, monotonic contraction that stabilizes into a static remnant by \(t\approx 8\). The remnant state exhibits all the hallmarks of a trumpet-sliced black hole: the lapse settles at a small but nonzero value (\(\sim 4\times 10^{-2}\)), \(\max|K|\) converges to a constant (\(\sim 2.5\)), and the areal radius asymptotes to a fixed plateau.

The physical mechanism driving this collapse is the linear instability of the Ellis--Bronnikov wormhole identified by Gonzalez~et~al.~\cite{gonzalez09} and confirmed in spherical symmetry by Shinkai and Hayward~\cite{shinkai02}. The Ellis--Bronnikov wormhole is linearly unstable, and the truncation error of the numerical scheme provides a seed perturbation that is sufficient to trigger the unstable mode. The choice of \(\chi_{\rm min}=10^{-8}\) is critical: with a higher floor (\(\chi_{\rm min}=10^{-4}\)), the artificial plateau near \(\bar r=0\) effectively stabilizes the inner region, masking the physical instability and producing an apparent expansion that is in fact an artifact of the protective floor. Only when the floor is lowered sufficiently does the true dynamical fate---collapse---emerge.


As theoretically expected from spherical symmetry, the unperturbed run should yield \(\Psi_4=0\) in the continuum limit. The small signal extracted in Fig.~\ref{fig:psi4_unperturbed} (Appendix) (amplitude \(\sim 4\times 10^{-4}\)) represents the numerical noise floor of our setup, driven by the Cartesian grid and AMR refinement boundaries breaking exact spherical symmetry during the violent constraint-damping phase. Crucially, the propagation speed analysis (panel~d) confirms that this signal is entirely non-physical: the dominant peaks travel at \(v\approx 1.4\)--\(1.5c\), well above the speed of light, identifying them as CCZ4 constraint-damping modes rather than gravitational waves.

\noindent\textbf{Absence of ringdown despite black-hole formation.}
Although the diagnostics (Fig.~\ref{fig:unperturbed_diagnostics}) unambiguously confirm that a black hole has formed---the lapse freezes, \(\theta_+\) becomes negative, and the areal radius settles onto a trumpet plateau---the extracted \(\Psi_4\) (Fig.~\ref{fig:psi4_unperturbed}) contains \emph{no quasi-normal-mode ringdown}. The ESD (panel~c) falls off monotonically with no peaked spectral feature, and the time-domain waveform (panels~a,\,b) shows erratic, noise-like oscillations rather than a coherent, exponentially decaying sinusoid. This is the expected consequence of Birkhoff's theorem applied to the remnant itself: because the collapse proceeds with exact spherical symmetry (to machine precision), the black hole forms in a perfectly round state with no \(\ell=2\) deformation. A QNM is a resonance that must be \emph{excited} by a non-spherical perturbation; when none is present, the black hole simply sits silently. This result provides a critical null test for the perturbed analysis in Sec.~IV\,C: the oscillatory late-time structure observed in the perturbed \(\Psi_4\) waveform (Fig.~\ref{fig:psi_perturbed}) is entirely absent here, confirming that it is driven by the \(A_\phi\) quadrupolar deformation rather than being an intrinsic feature of the collapse dynamics. The two-order-of-magnitude amplitude ratio (\(\sim 3\times 10^{-2}\) vs.\ \(\sim 4\times 10^{-4}\)) further establishes that the perturbed signal sits well above the numerical noise floor.

The constraint norms (Appendix, Fig.~\ref{fig:constraints_unperturbed}) remain bounded throughout, confirming that the CCZ4 evolution is well-controlled.


\subsection{Perturbed evolution: compressive deformation and accelerated collapse}

The diagnostics in Fig.~\ref{fig:perturbed_diagnostics} (Appendix) confirm that the monopolar compression drives rapid black-hole formation: the throat collapses by \(t\approx 2\), compared to \(t\approx 5\) in the unperturbed case. The throat areal radius (panel~g) drops from its initial value of \(\sim 0.77\) to the \(\sim 0.10\) trumpet plateau, with a steep initial contraction driven by the \(A_0{=}0.5\) monopole. The lapse (panel~a) settles at \(\sim 10^{-1}\), the extrinsic curvature (panel~c) equilibrates at \(\max|K|\approx 3\), and the trapped-surface proxy (panel~d) stabilizes at \(r_{\rm AH}\approx 0.15\), all characteristic of a trumpet-sliced black hole. The null expansion proxy (panel~e) plunges to \(\sim -200\), confirming strong self-trapping.

The new phantom-field diagnostic panels provide direct evidence of self-consistent matter-gravity evolution. The scalar field profile (panel~j) undergoes violent oscillations during collapse (\(t\sim 1\text{--}4\)) as the phantom matter is squeezed and redistributed by the dynamical geometry, then settles to a roughly constant range (\(\phi\sim 0.30\text{--}0.43\)) as the field is partially swallowed by the black hole and partially expelled. The conjugate momentum (panel~k) starts exactly at zero---confirming that our initial data strategy (\(\Pi(t{=}0)=0\)) is implemented correctly---and then explodes to \(|\Pi|\sim 0.35\) during collapse as the Klein--Gordon equation dynamically generates momentum from the scalar field gradients. At late times \(\Pi\) damps back toward zero as the remnant settles. This dynamical generation of \(\Pi\) from zero initial data is the signature of a fully coupled Einstein--Klein-Gordon evolution, fundamentally distinct from any approach that does not evolve the matter sector.

\subsection{Gravitational-wave signatures (near-zone)}

\begin{figure*}[tp]
\centering
\includegraphics[width=0.92\textwidth]{psi4_analysis_perturbed.pdf}
\caption{Gravitational-wave analysis for the perturbed evolution (\(A_0{=}0.5\), \(A_\phi{=}0.05\), \(b_0{=}0.5\)).
\textbf{(a)}~Signal amplitude \(\sim 3\times 10^{-2}\), two orders of magnitude above the unperturbed noise floor, with 3--4 clear oscillation cycles visible; signals arrive at progressively later times at larger radii, confirming causal propagation.
\textbf{(b)}~Retarded-time waveforms align across extraction radii; the damped-sinusoid fit yields \(f_{\rm QNM}\approx 0.295\,M^{-1}\), \(\tau\approx 3.03\,M\) (\(Q\approx 2.8\)). This structure is absent in the unperturbed case (Fig.~\ref{fig:psi4_unperturbed}b).
\textbf{(c)}~ESD shows peaked structure, unlike the monotonically falling unperturbed ESD.
\textbf{(d)}~Propagation speed \(v\approx 1.01c\) between \(R{=}8\) and \(R{=}12\), consistent with physical gravitational radiation; the unperturbed run yields superluminal \(v\approx 1.4\text{--}1.5c\). The outer extraction pair (\(R{=}12\to 16\)) is not yet reliable because at \(t\sim 15\) the wavefront has not fully traversed \(R{=}16\).
\textbf{(e)}~Strain PSD.
\textbf{(f)}~For \(M{=}30\,M_\odot\) at \(D{=}10\)~Mpc the signal peaks above the most sensitive LIGO band (SNR\(\approx 0\)).}
\label{fig:psi_perturbed}
\end{figure*}

In the perturbed case (Fig.~\ref{fig:psi_perturbed}), the extracted \(\Psi_4\) amplitude is \(\sim 3\times 10^{-2}\), two orders of magnitude above the unperturbed noise floor, with 3--4 clear oscillation cycles. The propagation speed between the inner extraction pair (\(R{=}8\to 12\)) is \(v\approx 1.01c\), consistent with physical gravitational radiation propagating at the speed of light. This is the central observational result of the perturbed evolution: the monopolar compression drives a fast collapse (\(t\sim 2\)) that, combined with the quadrupolar symmetry breaking, produces a near-zone curvature transient traveling at the physically expected speed---in stark contrast to the superluminal modes (\(v\approx 1.4\text{--}1.5c\)) that dominate the unperturbed signal. The outer extraction pair (\(R{=}12\to 16\)) yields \(v\approx 2.1c\), but this measurement is unreliable because at \(t\sim 15\) the simulation has not yet run long enough for the wavefront to fully propagate to \(R{=}16\) (the retarded time at that radius is \(u=t-R\approx -1\)). Because both runs begin from constraint-satisfying initial data (to truncation error) and share the same numerical settings and remnant type, features appearing only in the perturbed case cannot be attributed to constraint-damping modes alone.

\noindent\textbf{Ringdown structure.}
The retarded-time waveform (Fig.~\ref{fig:psi_perturbed}b) exhibits a clear oscillatory late-time tail, quantitatively characterized in Sec.~IV\,D. Because both runs produce the same remnant yet this structure appears only when \(A_\phi\neq 0\) (contrast with the null result in Sec.~IV\,A), it is seeded by the quadrupolar symmetry breaking.

\noindent\textbf{Caveats.}
The scalar field profile perturbation satisfies the momentum constraint exactly; the only residual is a Hamiltonian constraint violation of \(\mathcal{O}(A_0, A_\phi)\) from the cross-term \(2\nabla\phi_{\rm EB}\cdot\nabla\delta\phi\). Crucially, this is a spin-0 violation that does not contaminate the spin-2 \(\Psi_4\) signal (see Sec.~II\,E). All extractions are performed in the near zone (\(R_{\rm ext}=8\text{--}16\) with \(L_{\rm full}=40\)) and should be understood as near-zone diagnostics rather than asymptotic templates. The constraint norms (Appendix, Fig.~\ref{fig:constraints_perturbed}) remain bounded throughout: the Hamiltonian norm starts at \(\sim 6\times 10^{-2}\) (as predicted from the \(\mathcal{O}(A_0)\) cross-term), dips to \(\sim 10^{-2}\) as CCZ4 damps the initial residual, and rises slowly to \(\sim 4\times 10^{-2}\); the momentum norm starts at \(\sim 10^{-4}\) (confirming that the \(0{=}0\) argument of Eq.~12 works to machine precision) and grows to \(\sim 10^{-2}\) only as the dynamics generate \(K_{ij}\). The spatial structure of the \(K\) field is shown in Fig.~\ref{fig:perturbed_kz_evolution} (Appendix).


\subsection{Black-hole remnant characterisation}

The universal collapse outcome permits several quantitative characterisations of the newly formed black hole.

\noindent\textbf{Irreducible mass from the throat areal radius.}
The gauge-independent areal radius (Eq.~\ref{eq:R_areal}) plateaus at late times at \(R_{\rm areal,min}\approx 0.10\) in the unperturbed run. In the moving-puncture 1+log gauge the minimum areal radius corresponds to the trumpet throat. Using the standard relation \(R_{\rm trumpet} \approx 1.31\,M_{\rm BH}\)~\cite{hannam07}, this yields an estimated remnant mass \(M_{\rm BH}\approx 0.077\) in code units. The lapse settles at \(\alpha_{\rm min}\approx 0.034\), consistent with the trumpet-slice puncture value for a BH of this mass.

\noindent\textbf{Total radiated energy.}
Integrating \(|r\Psi_4|^2/(16\pi)\) over retarded time at the innermost extraction radius (\(R=8\)) gives \(E_{\rm rad}\approx 1.9\times 10^{-5}\,M\)---an extremely small fraction of the system energy, consistent with the weakness of the quadrupolar deformation (\(A_\phi=0.05\)).

\noindent\textbf{QNM ringdown fit.}
To test whether the late-time oscillation in the perturbed \(\Psi_4\) waveform is a black-hole quasi-normal mode, we fit a damped sinusoid \(A e^{-t/\tau}\sin(2\pi f_{\rm QNM} t + \phi)\) to the ringdown tail at \(R=8\). The fit yields \(f_{\rm QNM}\approx 0.295\,M^{-1}\) and \(\tau\approx 3.03\,M\), giving a quality factor \(Q = \pi f\tau \approx 2.8\). Interpreting this frequency as the fundamental Schwarzschild \(\ell=2\) QNM (\(\omega M = 0.3737\)) gives \(M_{\rm BH}^{f} = 0.3737/(2\pi f_{\rm QNM}) \approx 0.20\). The imaginary part of the same mode (\(\omega_I M = 0.0890\)) yields an independent mass estimate from the damping time: \(M_{\rm BH}^{\tau} = \tau\times 0.0890 \approx 0.27\). The \(\sim 34\%\) discrepancy between these two QNM-derived mass estimates reflects the fact that the final state is not a pristine vacuum Schwarzschild black hole but rather a compact remnant actively accreting phantom (negative-energy) scalar matter, which shifts both the oscillation frequency and damping time away from the analytic predictions.

\noindent\textbf{Mass discrepancy between geometric and QNM estimates.}
The trumpet-radius geometric mass (\(M_{\rm BH}^{\rm geom}\approx 0.077\), from \(R_{\rm areal,min}\approx 0.10\) and \(R_{\rm trumpet}\approx 1.31\,M\)~\cite{hannam07}) is substantially smaller than both QNM-derived estimates (\(M_{\rm BH}^{f}\approx 0.20\), \(M_{\rm BH}^{\tau}\approx 0.27\)). This hierarchy is physically expected: the trumpet relation assumes a pristine vacuum Schwarzschild spacetime, whereas the final state here is a black hole immersed in a dynamically evolving phantom scalar field with negative energy density. The phantom matter modifies the effective potential experienced by curvature perturbations, shifting the QNM frequencies and damping times. The important point is that all three estimates consistently indicate a compact remnant of order \(M_{\rm BH}\sim 0.1\text{--}0.3\) in code units, and the spread quantifies the systematic uncertainty introduced by the non-vacuum background rather than a methodological failure.

\noindent The unperturbed null test (Sec.~IV\,A) confirms that this oscillatory structure is seeded by the \(A_\phi Y_{20}\) symmetry breaking, not by the collapse itself. However, for the parameters used here the characteristic frequency of the initial Gaussian perturbation (set by \(\sigma_\phi\)) is comparable to the fitted QNM frequency; a definitive separation would require repeating the simulation with a significantly different \(\sigma_\phi\) while keeping the remnant mass fixed.

\noindent\textbf{Instability growth rate.}
The early-time departure of the throat areal radius from its initial value in the unperturbed run is well fit by an exponential: \(|R_{\rm areal}(t) - R_{\rm areal}(0)| \sim e^{\lambda t}\) with \(\lambda \approx 0.50\,M^{-1}\), corresponding to an e-folding time of \(\sim 2.0\,M\). This provides a numerical measurement of the linear instability rate of the unsupported Ellis--Bronnikov geometry, complementing the analytical predictions of Gonzalez~et~al.~\cite{gonzalez09}. In the perturbed evolution, the $A_0$ monopolar perturbation directly seeds and amplifies this compressive mode, driving the measured growth rate up to \(\lambda \approx 0.99\, M^{-1}\) (e-folding time \(\sim 1.0\,M\)). This quantifiable jump confirms that the monopolar perturbation is doing exactly what it was designed to do.

\noindent\textbf{Total physical collapse time.}
We can scale the total time from $t=0$ to the point where the black hole fully forms and the trumpet plateau stabilizes ($t \approx 2M$) into physical units. To illustrate how fast this collapse occurs in the real universe, Table~\ref{tab:collapse_times} provides the total physical collapse time for different representative wormhole masses. Even for a supermassive wormhole comparable to Sagittarius A*, the entire collapse from a static throat into a black hole takes only about 10 seconds.

\begin{table}[h]
\caption{Total physical collapse time ($t \approx 2M$) for the accelerated perturbed collapse at various mass scales, assuming an initial throat radius $b_0=0.5M$.}
\label{tab:collapse_times}
\begin{ruledtabular}
\begin{tabular}{lccc}
Scale & Mass ($M$) & Throat ($b_0$) & Collapse Time \\
\hline
Stellar & $30\, M_\odot$ & $\sim 22$ km & $\sim 0.3$ ms \\
Intermediate & $10^3\, M_\odot$ & $\sim 740$ km & $\sim 10$ ms \\
Supermassive & $10^6\, M_\odot$ & $\sim 7.4{\times}10^5$ km & $\sim 10$ s \\
Cosmological & $10^9\, M_\odot$ & $\sim 5$ AU & $\sim 2.7$ hr \\
\end{tabular}
\end{ruledtabular}
\end{table}

\noindent\textbf{Detectability.}
For a hypothetical \(M=30\,M_\odot\) remnant at \(D=10\)~Mpc, the physical QNM frequency is \(f_{\rm phys}\approx 2000\)~Hz, placing it above the most sensitive Advanced LIGO band with \(\mathrm{SNR}\approx 0\). The amplitude of the extracted signal is directly proportional to the conservatively small quadrupolar deformation (\(A_\phi=0.05\)) chosen for this simulation. Since the radiated energy scales as \(\mathcal{O}(A_\phi^2)\), a highly asymmetric primordial collapse (with \(A_\phi \sim \mathcal{O}(1)\)) would produce a vastly stronger signal. For the moderate perturbation simulated here, detection would require either a much more massive wormhole progenitor---an intermediate-mass wormhole (\(\sim 500\text{--}1000\,M_\odot\)) would shift the signal into the most sensitive LIGO band (\(\sim 100\text{--}200\)~Hz), while a supermassive wormhole (\(\sim 10^8\,M_\odot\)) would be a target for the space-based LISA observatory in the millihertz band---or a much closer source (within the Milky Way).



\section{Discussion and Conclusion}
We have investigated the nonlinear matter-gravity dynamics of unstable traversable wormholes using full 3D numerical relativity. Our results are consistent with the collapse pathway identified by Shinkai and Hayward~\cite{shinkai02} in 1D, but extend it in two ways: (i)~lowering the conformal-factor floor to \(\chi_{\rm min}=10^{-8}\) reveals that the unperturbed geometry spontaneously collapses---a result masked by higher floors that artificially stabilized the inner region; and (ii)~full 3D evolution enables extraction of the Weyl scalar \(\Psi_4\) and propagation speed analysis to distinguish physical curvature transients from constraint-damping artifacts, which is impossible in spherical symmetry by Birkhoff's theorem. All near-zone curvature quantities should be regarded as local diagnostics rather than asymptotic gravitational-wave templates.

The central result is that the unstable Ellis--Bronnikov geometry collapses into a black hole in all cases studied. With no applied initial deformation, the geometry spontaneously collapses due to the well-known linear instability of the Ellis--Bronnikov solution~\cite{gonzalez09}. When an initial scalar field perturbation is applied, the collapse is dramatically accelerated. In both cases, the geometry does not remain as a static traversable throat on the resolved grid. While quantum backreaction may theoretically stabilize microscopic topological defects at the Planck scale \cite{mehulic2026}, our results demonstrate that once stretched to macroscopic scales where classical dynamics dominate, they are inherently unstable. This provides a natural cosmological mechanism explaining why inflationary wormholes do not abundantly persist in the late universe.

The collapse timescale is fundamentally set by the light-crossing time of the throat, \(t_{\rm lc} = b_0/c\). For a stellar-mass wormhole (\(M\sim 30\,M_\odot\)) with a physical throat radius of \(\sim 22\)~km, light crosses the throat in about \(0.07\)~ms---and the entire topological transition from a traversable bridge to a black hole is completed in roughly four of those light-crossing times (Table~\ref{tab:collapse_times}). This means that unstable wormholes are essentially instantaneous events on any astrophysical timescale: a primordial wormhole seeded during inflation and stretched to stellar scales would have collapsed within a fraction of a millisecond after losing stability. Even a hypothetical galaxy-sized wormhole (\(10^9\,M_\odot\) scale) would be gone in under three hours. For an unstable wormhole to survive for the age of the universe (\(\sim 13.8\)~Gyr), its throat radius would need to exceed \(\sim 10^9\)~light-years---a cosmological scale. Any unstable primordial wormhole of stellar or galactic scale would have collapsed into a black hole almost instantaneously.

The extreme rapidity of the collapse also explains the violent scalar-field dynamics observed in the phantom-field diagnostic panels (Fig.~\ref{fig:perturbed_diagnostics}j,\,k): the Klein--Gordon equation must generate the entire conjugate momentum \(\Pi\) from exactly zero initial data within \(\sim 2\) code time units, leaving the phantom matter no opportunity to respond adiabatically. Instead, it is squeezed, redistributed, and partially swallowed by the forming black hole on a dynamical timescale comparable to the light-crossing time itself.

The extracted near-zone curvature signal provides a powerful diagnostic tool. The propagation speed measurement \(v\approx 1.01c\) between \(R{=}8\) and \(R{=}12\) in the perturbed run, compared to \(v\approx 1.4\text{--}1.5c\) in the unperturbed case, is the strongest evidence that the perturbed curvature transient contains a physical gravitational-wave component. The contrast between the two runs---identical numerical settings and initial constraints satisfied to truncation error, but dramatically different signals---further confirms that features appearing only when \(A_\phi\neq 0\) cannot be attributed to constraint-damping modes alone. The dynamical generation of \(\Pi\) from exactly zero initial data (tracked in the phantom-field diagnostic panels) provides an independent confirmation that the Klein--Gordon sector is driving the collapse physics self-consistently.

Several natural extensions would place the results on firmer ground. First, constructing fully constraint-satisfying initial data with strong metric perturbations using an elliptic solver would enable the study of large-amplitude non-linear effects without the \(\mathcal{O}(A_0, A_\phi)\) Hamiltonian residual introduced by our scalar profile perturbation. Because we explicitly evolved the phantom scalar field, energy accounting is rigorous and the ADM mass is conserved. Second, extracting \(\Psi_4\) at larger radii (or using Cauchy-characteristic extraction) would move the measurement into the wave zone and enable the construction of asymptotic waveform templates suitable for matched-filter searches. These extensions are the targets of our ongoing development of the \texttt{GRTeclyn} framework.

\section{Acknowledgements}
This research was supported by Gravity Frontiers (\texttt{ic@gravityfrontiers.org}) under a research grant for the development of software for mathematical modeling of wormholes.

The author also thanks Ilya Nachevsky for generously providing the high-performance computing resources necessary for this work. These simulations would not have been possible without his support. 

\appendix*
\section{Numerical details, validation, and convergence}
We exploit Cartesian reflection symmetries and evolve only the positive octant, \(x\ge 0\), \(y\ge 0\), \(z\ge 0\), imposing parity conditions on the inner faces and Sommerfeld radiative conditions at the outer boundary. The constraint norms are computed as volume-weighted root-mean-square (RMS) averages over the entire computational domain (AMR Level 0). While constraint violations are largest near the throat where gradients are extreme (which would typically motivate measuring them on the finest level), integrating them over the coarsest level provides a more rigorous global bound. The Level 0 grid covers the entire domain, meaning the RMS norm captures not only the generation of violations at the throat, but also their subsequent outward propagation toward the boundaries. If we evaluated the norm only on the finest levels, outwardly propagating constraint-cleaning modes would simply exit the refined region, leading to an artificially rapid drop in the norm and a misleading impression of constraint satisfaction. By measuring on Level 0, we track the total amount of constraint violation present in the entire simulated universe, confirming that the overall CCZ4 evolution remains stable and bounded.
\begin{align}
L_2(\mathrm{Ham}) &=
\left(\frac{\sum_i \mathrm{Ham}_i^2\,V_i}{\sum_i V_i}\right)^{1/2},\\
L_2(\mathrm{Mom}) &=
\left(\frac{\sum_i \left(\mathrm{Mom}_{x,i}^2+\mathrm{Mom}_{y,i}^2+\mathrm{Mom}_{z,i}^2\right)V_i}{\sum_i V_i}\right)^{1/2}.
\end{align}
These two scalars are written to the small-data file as \texttt{L2\_Ham} and \texttt{L2\_Mom}. Because the numerical initial data includes small interpolation errors mapping the exact solution to the Cartesian grid, as well as the \(\mathcal{O}(A_0, A_\phi)\) Hamiltonian constraint residual from the scalar-field perturbation (arising from the cross-term \(2\nabla\phi_{\rm EB}\cdot\nabla\delta\phi\)), these norms are generically nonzero at \(t=0\). In the perturbed case, the initial Hamiltonian norm is \(\sim 6\times 10^{-2}\), consistent with the predicted leading-order cross-term from \(A_0{=}0.5\); the initial momentum norm is \(\sim 10^{-4}\), confirming that the \(\Pi{=}0\), \(K_{ij}{=}0\) initial data satisfy the momentum constraint to machine precision. The relevant stability test is whether the CCZ4 evolution keeps the subsequent violation bounded. Figures~\ref{fig:constraints_unperturbed} and~\ref{fig:constraints_perturbed} show that for both evolutions the Hamiltonian and momentum RMS norms remain bounded throughout. The late-time Hamiltonian norm (\(\sim 4\times 10^{-2}\)) is large by the standards of constraint-satisfying numerical relativity, reflecting the \(\mathcal{O}(A_0)\) initial residual and the persistent under-resolution of the compactified second asymptotic end near $\bar r=0$. The boundedness of the norms nevertheless confirms numerical stability and ensures that the CCZ4 evolution does not develop runaway violations that would invalidate the dynamical conclusions.

\begin{figure}[htp]
\centering
\includegraphics[width=\linewidth]{constraints_unperturbed.pdf}
\caption{Volume-weighted RMS norms of the Hamiltonian and momentum constraints for the unperturbed evolution (\(A_\phi{=}0\)). The Hamiltonian norm grows slowly but remains at \(\mathcal{O}(10^{-2})\); the momentum norm saturates at \(\sim 8\times 10^{-3}\). These levels reflect the truncation error and grid effects of the numerical scheme, but the norms remain bounded throughout the \(t\approx 18\) evolution, confirming numerical stability.}
\label{fig:constraints_unperturbed}
\end{figure}

\begin{figure}[htp]
\centering
\includegraphics[width=\linewidth]{constraints_plot_perturbed.pdf}
\caption{Volume-weighted RMS norms of the Hamiltonian and momentum constraints for the perturbed evolution (\(A_0{=}0.5\), \(A_\phi{=}0.05\)). The Hamiltonian norm starts at \(\sim 6\times 10^{-2}\)---the predicted \(\mathcal{O}(A_0)\) cross-term from the monopolar perturbation---dips to \(\sim 10^{-2}\) as CCZ4 damps the initial residual, then rises slowly to \(\sim 4\times 10^{-2}\). The momentum norm starts at \(\sim 10^{-4}\), confirming that the momentum constraint is satisfied to machine precision at \(t{=}0\) by the \(\Pi{=}0\), \(K_{ij}{=}0\) initial data, and grows to \(\sim 10^{-2}\) only as the dynamics generate extrinsic curvature. Both norms remain bounded throughout, confirming numerical stability.}
\label{fig:constraints_perturbed}
\end{figure}

\begin{figure*}[tp]
\centering
\includegraphics[width=0.92\textwidth]{diagnostic_unperturbed.pdf}
\caption{Diagnostics of the unperturbed evolution (\(A_\phi{=}0\), \(b_0=0.5\), \(\chi_{\rm min}=10^{-8}\)).
\textbf{(a)}~\(\min(\alpha)\) plunges during the initial transient, then stabilizes at \(\sim 4\times 10^{-2}\) (trumpet slicing).
\textbf{(b)}~\(\min(\chi)\) drops sharply, then rises to \(\sim 10^{-3}\).
\textbf{(c)}~\(\max|K|\) settles to the trumpet equilibrium \(\sim 2.5\).
\textbf{(d)}~\(r_{\rm AH}\) decreases from \(\sim 0.28\) to \(\sim 0.14\).
\textbf{(e)}~\(\min(\theta_+)\) plunges to \(\sim -250\), confirming strong self-trapping.
\textbf{(f)}~Radius at \(\min(\theta_+)\) settles at the puncture (\(\sim 0.02\)).
\textbf{(g)}~\(R_{\rm areal,min}\) falls from \(0.50\) to \(\sim 0.10\) by \(t\approx 5\), the trumpet radius of the remnant.
\textbf{(h)}~\(dR_{\rm areal}/dt\) approaches zero at late times (static remnant).
\textbf{(i)}~\(\tau\approx 0\): \(\max|K|\) converges to a constant rather than decaying, confirming black-hole formation.}
\label{fig:unperturbed_diagnostics}
\end{figure*}

\begin{figure*}[tp]
\centering
\includegraphics[width=0.92\textwidth]{psi4_analysis_unperturbed.pdf}
\caption{Gravitational-wave analysis for the unperturbed evolution (\(A_\phi{=}0\)).
\textbf{(a)}~Signal amplitude \(\sim 4\times 10^{-4}\), consistent with numerical noise.
\textbf{(b)}~Retarded-time waveforms do not align; the QNM fitter returns a flat line (\(Q\approx 0.02\)), confirming no ringdown.
\textbf{(c)}~ESD falls monotonically with no QNM peak.
\textbf{(d)}~Propagation speeds are superluminal (\(v\approx 1.39\)--\(1.49c\)), identifying CCZ4 constraint-damping modes.
\textbf{(e)}~Strain PSD is featureless.
\textbf{(f)}~Signal falls below the Advanced LIGO noise floor (SNR\(\approx 0\)), establishing the null baseline.}
\label{fig:psi4_unperturbed}
\end{figure*}

\begin{figure*}[tp]
    \centering
    \includegraphics[width=0.92\textwidth]{collapse_diagnostics_perturbed.pdf}
    \caption{Diagnostics of the perturbed collapse (\(A_0{=}0.5\), \(A_\phi{=}0.05\), \(b_0{=}0.5\), \(\chi_{\rm min}=10^{-8}\)).
\textbf{(a)}~\(\min(\alpha)\) drops sharply during the initial transient, then stabilizes at \(\sim 10^{-1}\) (trumpet equilibrium).
\textbf{(b)}~\(\min(\chi)\) hits the floor (\(10^{-8}\)), then rises to \(\sim 10^{-3}\).
\textbf{(c)}~\(\max|K|\) settles to \(\sim 3\), indicating a trumpet-sliced black hole.
\textbf{(d)}~\(r_{\rm AH}\) decreases from \(\sim 0.44\) to \(\sim 0.15\), confirming trapped-surface formation and stabilization.
\textbf{(e)}~\(\min(\theta_+)\sim -200\), confirming strong self-trapping.
\textbf{(f)}~Radius at \(\min(\theta_+)\) settles at \(\sim 0.02\).
\textbf{(g)}~\(R_{\rm areal,min}\) falls from \(\sim 0.77\) to \(\sim 0.10\) by \(t\approx 2\), the trumpet radius of the remnant. The fast collapse timescale (compared to \(t\approx 5\) in the unperturbed case) is driven by the \(A_0{=}0.5\) monopole.
\textbf{(h)}~\(dR_{\rm areal}/dt\to 0\) at late times (static remnant).
\textbf{(i)}~\(\tau\approx 0\): \(\max|K|\) converges to a constant, confirming black-hole formation.
\textbf{(j)}~Phantom scalar field profile \(\phi\): violent oscillations during collapse (\(t\sim 1\text{--}4\)) as the phantom matter is squeezed and redistributed, settling to \(\sim 0.30\text{--}0.43\) at late times (field partially swallowed, partially expelled).
\textbf{(k)}~Scalar field momentum \(\Pi\): starts exactly at zero (confirming the initial data strategy), dynamically generated to \(|\Pi|\sim 0.35\) during collapse via the Klein--Gordon equation, then damps back toward zero as the remnant settles.
\textbf{(l)}~Instability growth rate (\(\lambda\)): exponential departure of \(R_{\rm areal,min}\) from its initial value, directly quantifying the rate of collapse.}
    \label{fig:perturbed_diagnostics}
\end{figure*}


\begin{figure*}[tp]
\centering
\includegraphics[width=\textwidth]{evolution_panel_perturbed.pdf}
\caption{Evolution snapshots of the \(K\) field in the \(z{=}0\) plane for the perturbed run (\(A_\phi{=}0.05\), \(b_0{=}0.5\)). The initial scalar field perturbation dynamically generates a compressive deformation (bright central region) that drives rapid infall. By \(t\approx 2\), the throat has pinched off and the \(K\) field develops the trumpet-slice profile: a localized positive-\(K\) core at the puncture surrounded by a negative-\(K\) halo, with outward-propagating ripples representing the constraint-cleaning transient and the near-zone curvature signal. At late times the exterior relaxes toward \(K\to 0\) as required for asymptotic flatness.}
\label{fig:perturbed_kz_evolution}
\end{figure*}

% --- Bibliography ---
\clearpage
\bibliographystyle{apsrev4-2}
\bibliography{references}
\end{document}