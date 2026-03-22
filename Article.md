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
\title{\MakeTextUppercase{Nonlinear Evolution of Primordial Wormhole-Like Topological Defects
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
Wormholes are topological bridges permitted by General Relativity, but the best-known traversable solutions require exotic matter that may not survive beyond the earliest stages of the universe. If primordial wormhole-like defects were seeded from quantum foam and subsequently stretched to macroscopic scales during inflation, their later evolution would be governed not by a sustained exotic source but by the unsupported relaxation of the geometry itself. In this work, we investigate the nonlinear dynamics of such unsupported wormhole-like topological defects under localized extrinsic-curvature perturbations using full 3D numerical relativity with \texttt{GRTeclyn}. We find three distinct dynamical regimes: (i)~unperturbed evolution ($A_0{=}A_2{=}0$), in which the bare topology undergoes vacuum expansion driven solely by its inherent constraint defect; (ii)~sub-critical relaxation under a moderate kick, where the geometry disperses without detected trapping; and (iii)~super-critical self-trapping, in which a strong implosive kick drives topological pinch-off while the extracted near-zone curvature is dominated by transient junk radiation. These results clarify the vacuum dynamical fate of unsupported exotic geometries and provide a concrete framework for studying relic topological defects without assuming persistent exotic matter.
\end{abstract}


% The \maketitle command generates the title block AFTER the abstract
\maketitle

% --- Main Content of the Proposal ---
\section{Introduction}
Wormholes—topological bridges connecting distinct regions of spacetime or separate universes—stand among the most fascinating predictions of General Relativity. The geometric possibility of such structures was first identified by Ludwig Flamm in 1916 \cite{flamm16}, shortly after the discovery of the Schwarzschild solution. The concept was formally introduced into physics by Einstein and Rosen in 1935 as an attempt to model elementary particles as non-singular field configurations \cite{einstein35}, and in 1957 John Archibald Wheeler coined the term ``wormhole'' while investigating the topological foam of spacetime at the quantum scale \cite{misner57}. The first explicit traversable wormhole solutions supported by a phantom scalar field were found independently by Ellis~\cite{ellis73} and Bronnikov~\cite{bronnikov73} in 1973; this Ellis--Bronnikov geometry remains the canonical example studied in the present work.

For decades, these solutions were considered topological curiosities, primarily because the standard Einstein-Rosen bridge is non-traversable; the throat pinches off faster than even a photon can cross it. The modern era of wormhole physics accelerated when Morris and Thorne demonstrated that static, traversable wormholes are valid solutions to the Einstein equations, provided the throat is threaded by ``exotic matter'' that violates the null energy condition \cite{morris88}. Consequently, much of the literature has focused on how such stress-energy might be engineered, minimized, replaced in modified gravity \cite{visser95, lobo05}, or sourced by astrophysical dark matter profiles \cite{garattini2026}. Additionally, recent stationary models demonstrate that rapid rotation can arbitrarily reduce the amount of exotic matter required to keep the throat open \textbf{\cite{uemichi2026}}, pointing toward spin as a potential stabilizing mechanism. That perspective is natural if one is attempting to engineer or stabilize a wormhole for transport, but it is not the only physically relevant question.

The motivation for the present work, however, diverges entirely from the search for stability. In Wheeler's picture of spacetime foam, transient microscopic wormholes may be produced at very early times, and inflation may in principle stretch some topological defects to much larger scales \cite{kardashev07, garcia16}. If such relic objects ever existed, the exotic component that originally sustained them need not remain present at late times; it may have been tied to an inflationary sector, a transient field configuration, or some other early-universe mechanism that subsequently decayed away. In that case the relevant physical problem is not how to maintain a supported traversable wormhole indefinitely, but how an unsupported wormhole-like geometry relaxes, disperses, or collapses once left to evolve under ordinary vacuum gravity. Indeed, recent theoretical work by \citet{dimaschko25} suggests that any bound state of a traversable wormhole is non-stationary. Rather than trying to keep the bridge open, this work investigates the subsequent dynamics of unsupported topological defects under external perturbation.

The dynamical fate of such traversable wormholes was first established in a pioneering numerical study by Shinkai and Hayward \cite{shinkai02}. Evolving the classic Ellis--Bronnikov metric in 1D spherical symmetry, they demonstrated that these structures are highly unstable against perturbations. They observed that the wormhole throat suffers a bifurcation of horizons, leading either to an inflationary expansion (when the perturbation is rarefactive, $K<0$) or a rapid collapse into a black hole (when the perturbation is compressive, $K>0$), as conceptually illustrated in Fig.~\ref{fig:wormhole_collapse_stages}. Crucially, even a minute compressive perturbation---representing, for instance, a traveler attempting to traverse the bridge---inevitably forces the wormhole to collapse into a black hole, sealing off causal connection between the two regions. Their work elegantly proved that exotic topological structures can be rigorously studied using the same local, dynamical methods applied to standard black holes, effectively unifying the two phenomena. 

While Shinkai and Hayward established the ultimate instability and evolutionary paths (contraction vs. expansion) of the wormhole topology, their strict assumption of 1D spherical symmetry inherently precluded the study of non-spherical dynamics and any extracted curvature signal. By Birkhoff's theorem, a purely spherically symmetric vacuum evolution cannot radiate. To determine how unsupported wormhole-like defects behave once symmetry is broken, one must therefore move beyond 1D restrictions and perform full 3D, non-linear dynamical simulations.

Simulating an unsupported exotic geometry in full 3D is computationally demanding: the extreme curvature gradients at the wormhole throat require deep adaptive mesh refinement, and the long evolution times needed to distinguish dispersal from collapse push memory and throughput requirements well beyond what CPU-only codes can deliver in a practical timeframe. A useful methodological precedent is the recent 3D study by \citet{clough24}, who evolved the Alcubierre warp-bubble geometry as an initial-value problem in vacuum, demonstrating that numerical relativity can follow the nonlinear fate of an exotic spacetime once the sustaining source is absent. Following a similar philosophy, we study unsupported wormhole dynamics with the \texttt{GRTeclyn} framework—a GPU-native numerical-relativity code built upon the \texttt{AMReX} infrastructure \cite{amrex} and descended from \texttt{GRChombo} \cite{Andrade2021}—which implements block-structured adaptive mesh refinement explicitly optimized for high-performance GPU acceleration on H100-class hardware. This capability is essential: the parameter surveys reported here require dozens of full 3D evolutions, each utilizing up to 8 GPUs, a workload that would be prohibitive without dedicated GPU acceleration. Our aim is not to assume away the missing exotic matter, but to understand the nonlinear fate of the residual topological geometry once that support is gone.

Moreover, any primordial wormhole that survived inflationary stretching would not remain perfectly spherically symmetric: tidal interactions with surrounding density perturbations, the stochastic nature of reheating, and subsequent structure formation all introduce asphericity. The extrinsic-curvature perturbation we impose at $t{=}0$ therefore models this physically expected post-inflationary deformation, rather than being a purely artificial device.

\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{wormhole_collapse_stages.pdf}
\caption{Conceptual 3D embedding diagram illustrating the dynamical evolutionary paths of an unstable traversable wormhole. \textbf{Top Row (a--c):} The super-critical branch. A strong compressive perturbation ($K>0$) drives the throat to constrict (b) until its topology dynamically pinches off, resulting in the formation of a black hole (c) that seals the two universes. \textbf{Bottom Row (d--f):} The sub-critical and unperturbed branches. With no perturbation or with a moderate kick, the resolved exterior geometry expands and weakens (e) as the negative spatial Ricci curvature, no longer balanced by exotic matter, drives volumetric expansion ($K<0$). We claim relaxation toward a flatter weak-curvature state on the resolved grid; the compactified second asymptotic end at $\bar r=0$ remains a coordinate artifact on the finite domain, precluding claims about the global topology of the full two-ended manifold (f).}
\label{fig:wormhole_collapse_stages}
\end{figure}

The ultimate goal of this research is to investigate the nonlinear dynamics of unsupported wormhole-like defects and to determine which perturbations lead to dispersal, which lead to trapping, and what gravitational-wave signal accompanies each regime. The advent of GW astronomy by the LIGO-Virgo-KAGRA collaborations \citep{abbott16,abbott18,abbott21,abbott23} has moved the search for Exotic Compact Objects (ECOs) \cite{cardoso19} from pure theory to observational reality, but the first theoretical task is to understand the dynamics of candidate geometries themselves. We extract the Weyl scalar \(\Psi_4\) at multiple coordinate radii in the near zone; retarded-time alignment of the waveforms across these radii reveals an outgoing curvature transient in the sub-critical regime whose amplitude far exceeds the numerical noise floor, while the super-critical regime is dominated by junk radiation and constraint-cleaning transients. Because the initial extrinsic-curvature kick does not satisfy the momentum constraint, the extracted signals inevitably contain a contribution from constraint-cleaning modes; a definitive separation of physical radiation from this contamination awaits constraint-satisfying initial data. Nevertheless, these near-zone diagnostics clearly distinguish the dynamical regimes from unperturbed vacuum expansion through sub-critical dispersal to super-critical pinch-off.

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
For the experiments studied here, the dynamical evolution is governed by the pure vacuum Einstein equations, implemented within the \texttt{GRTeclyn} framework using the conformal and covariant Z4 (CCZ4) formulation. Our objective is to start from the geometry of an Ellis--Bronnikov wormhole at $t=0$, subject to an initial extrinsic-curvature perturbation. The evolution is then driven by the nonlinear adjustment of this unsupported geometry via the CCZ4 constraints.

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
This smoothly vanishes at the origin ($\chi \to 0$), providing a regular way to represent the compactified second asymptotic end on the Cartesian grid. For the unsupported runs analyzed in this manuscript, the initial data passed to the solver consists only of the CCZ4 geometric variables with \(\tilde{\gamma}_{ij} = \delta_{ij}\); the explicitly supported scalar-field evolution is implemented separately in the codebase and is not used here.

\subsection{Dynamical Triggering via Extrinsic Curvature Perturbation}

The classic Ellis--Bronnikov traversable wormhole is an exact solution to the Einstein equations sourced by a massless phantom scalar field that violates the Null Energy Condition and thereby supports the throat. In this work, we investigate the dynamical fate of this topology in the limit where the phantom source is suddenly removed. Following the same ``sudden approximation'' employed by \citet{clough24} for the Alcubierre bubble, we initialize the exact spatial geometry of the Ellis--Bronnikov wormhole at $t=0$ but deliberately set the stress-energy tensor to zero ($T_{\mu\nu}=0$), rather than defining and dynamically depleting an active matter field.

By starting with this "empty" spacetime, we are directly simulating the bare residual topology the exact moment after the supporting matter has been removed. Consequently, the initial slice purposefully carries a Hamiltonian constraint defect that corresponds exactly to the missing exotic matter. This allows the CCZ4 formulation to naturally evolve the vacuum relaxation of the unsupported geometry.

In the ADM convention used throughout this paper, the trace of the extrinsic curvature $K$ encodes the local rate of change of spatial volume: $K>0$ corresponds to volumetric contraction (compression) and $K<0$ to volumetric expansion.

When evolved with no perturbation ($A_0{=}A_2{=}0$), the unsupported geometry \emph{expands}: the negative spatial Ricci curvature ($R^{(3)}<0$) at the throat, no longer balanced by the exotic source, drives the local volume outward. This establishes that the default vacuum fate of the bare topology is dispersive, not static. To map the full phase space---in particular, the threshold between this natural expansion and genuine gravitational collapse---we introduce a localized, tunable compressive perturbation ($K>0$) at $t{=}0$:
\begin{equation}
    K(\bar{r}, \theta, \phi) = \left[A_0 + A_2 Y_{20}(\theta, \phi)\right]
    \exp\left(-\frac{\bar{r}^2}{\sigma^2}\right).
\label{eq:kick_profile}
\end{equation}
The monopolar amplitude $A_0$ controls the overall compressive strength. All kicks in this work use $A_0\ge 0$ (compressive). For moderate amplitudes the natural expansion wins and the geometry disperses (sub-critical); for sufficiently large $A_0$ the compressive momentum overcomes the expansion and drives the throat through dynamical pinch-off into black-hole formation (super-critical).

The quadrupolar amplitude $A_2$ explicitly breaks the exact spherical symmetry of the initial data. Because Birkhoff's theorem forbids gravitational radiation from a spherically symmetric spacetime, this $Y_{20}$ deformation is required to produce a physical $\Psi_4$ signal. Physically, perfect spherical symmetry is unrealistic for any post-inflationary relic, so the quadrupolar component models the asphericity that tidal interactions and density perturbations would naturally impart.

\subsection{Apparent Horizon Detection}

During the period in which this work was carried out, \texttt{GRTeclyn} did not yet include a production-ready full apparent-horizon finder for these runs. To distinguish genuine topological pinch-off and black-hole-like self-trapping from gauge effects such as lapse collapse, we therefore implemented a geometric trapped-surface proxy on coordinate spheres centered on the throat. At each time step, we scan outward from the origin and estimate the expansion of outgoing null rays, \(\theta_+\). For a conformally flat spatial metric \(\gamma_{ij}=\chi^{-1}\delta_{ij}\), the spherical expression used in the code is
\begin{equation}
    \theta_+ = \frac{2\sqrt{\chi}}{r} - \frac{\partial_r \chi}{\sqrt{\chi}} + \tilde{A}_{rr} - \frac{2}{3}K.
\end{equation}
The formation of a trapped surface is identified when \(\theta_+ \le 0\). The maximum radius where this condition holds provides a proxy for the apparent-horizon radius, \(r_{\rm AH}\). Although this quantity does not replace a full elliptic horizon finder, it gives a practical collapse diagnostic for the parameter study reported here.

\subsection{Gauge Conditions}

To ensure stable evolution across both dispersive and collapsing branches, we employ robust gauge choices. 

The shift vector $\beta^i$ is evolved using the Gamma-driver condition, starting from $\beta^i(t{=}0) = 0$. Because the initial data possess no rotation and the primary dynamics are radial, the vanishing initial shift avoids gauge-induced artifacts during the early constraint-cleaning phase and the subsequent gravitational-wave extraction.

For the lapse function, we employ the $1+\log$ slicing condition, which evolves as $\partial_t \alpha - \beta^i \partial_i \alpha = -2\alpha K$. Two initialization strategies are used depending on the dynamical regime:
\begin{align}
    \text{Super-critical:}\quad & \alpha(t{=}0) = 1, \label{eq:lapse_unit}\\
    \text{Sub-critical \& unperturbed:}\quad   & \alpha(t{=}0) = \sqrt{\chi}. \label{eq:lapse_precollapsed}
\end{align}
In the super-critical regime, a strong initial $K$ kick drives the $1{+}\log$ condition to rapidly collapse the lapse toward zero near the throat, providing a natural singularity-avoidance mechanism; a unit initial lapse is therefore sufficient and avoids introducing gauge bias.
In the sub-critical regime (moderate $K$ kick) and in the unperturbed regime ($A_0{=}A_2{=}0$), the $K$-driven lapse collapse is either too slow or entirely absent, and the lapse remains near unity in regions where $\chi$ is very small. The CCZ4 evolution terms involving $1/\chi$ can then produce runaway numerical errors before the gauge has time to respond. The pre-collapsed profile $\alpha=\sqrt{\chi}$ seeds the lapse at a value already compatible with $1{+}\log$ equilibrium and suppresses the effective evolution rate in the compactified inner region of the wormhole, allowing the subsequent evolution to proceed stably without biasing the exterior dynamics where $\chi\approx 1$.





\section{Numerical Setup}

Numerical evolutions were performed using the \texttt{GRTeclyn} codebase, which is currently under active development. Built upon the \texttt{AMReX} framework~\cite{amrex}, \texttt{GRTeclyn} implements block-structured Adaptive Mesh Refinement (AMR) and is explicitly optimized for high-performance GPU acceleration. Production runs targeted H100-class GPUs (corresponding to \texttt{CUDA\_ARCH=90} in the build configuration) and used one MPI rank per GPU, scaling up to 8 GPUs for the largest calculations. In practice, the simulation workflow follows the same broad philosophy as other recent vacuum relaxations of exotic geometries, including the unsupported Alcubierre-bubble study of \citet{clough24}: evolve the bare geometry directly and diagnose whether it disperses, self-traps, or converts its initial defect into propagating transients.

The full physical computational domain spans a coordinate length of \(L_{\text{full}} = 40\), covered by a coarse grid of \(N_{\text{full}} = 160\) cells, yielding a coarse resolution of \(dx_{\text{coarse}} = 0.25\). To drastically reduce the computational cost, we exploit the Cartesian reflection symmetries of the single-throat initial data and evolve only the positive octant (\(x \ge 0, y \ge 0, z \ge 0\)), effectively modeling just one eighth of the physical volume. Parity symmetry conditions are imposed at the inner reflection boundaries, while Sommerfeld radiative boundary conditions are enforced at the outer edges.

To capture the extreme curvature gradients developing during the collapse while maintaining tractability, we employ up to 6 levels of 2:1 adaptive mesh refinement (AMR). The mesh is dynamically regridded based on the gradients of the conformal factor ($\chi$), allowing the code to focus resolution precisely where the wormhole throat is pinching off. The finest grid spacing ranges from $dx_{\rm fine}\approx 7.8\times10^{-3}$ (5 levels) to $dx_{\rm fine}\approx 3.9\times10^{-3}$ (6 levels), depending on the run. Time stepping is handled using a 4th-order Runge-Kutta scheme, with a conservative Courant factor (dt\_multiplier) reduced to 0.05 to maintain stability during the constraint relaxation and the abrupt collapse branches.

We focus on three representative evolutions---unperturbed, sub-critical, and super-critical---whose full parameter sets are given at the beginning of Sec.~IV.


\section{Results}

Before presenting the dynamics we summarize the three representative configurations and the key numerical choices that differ between them.

\noindent\textbf{Diagnostic sampling.}
All collapse diagnostics---$\min(\alpha)$, $\min(\chi)$, $\max|K|$, the trapped-surface proxy $r_{\rm AH}$, and $\min(\theta_+)$---are evaluated on the \emph{finest available AMR level} at each coarse time step. This ensures that the extreme gradients developing in the throat region are fully resolved in the reported quantities, rather than being smoothed by the coarse base grid. Note that $\max|K|$ reports the absolute value; the signed spatial structure of $K$ (distinguishing compressive $K>0$ from expansive $K<0$ regions) is shown separately in the evolution snapshots.

\noindent\textbf{Super-critical configuration.}
For the collapsing branch we use throat radius $b_0=0.2$, a strong monopolar kick $A_0=10.0$ with quadrupolar deformation $A_2=0.05$ and width $\sigma_K=0.2$, yielding rapid self-trapping. The initial lapse is set to $\alpha(t{=}0)=1$ (type~0); the strong $K$ kick triggers immediate $1{+}\log$ lapse collapse, which freezes the singular region and stabilizes the evolution. We employ $\kappa_1=0.1$, $\sigma_{\rm KO}=0.5$, $L_{\rm full}=40$, $N_{\rm full}=160$, and 6 levels of AMR ($\text{max\_level}=6$, $dx_{\rm fine}\approx 3.9\times10^{-3}$). The conformal-factor floor is set to $\chi_{\rm min}=10^{-10}$.

\noindent\textbf{Sub-critical configuration.}
For the dispersive branch with a moderate perturbation we use throat radius $b_0=0.5$ with no monopolar kick ($A_0=0.0$), a quadrupolar deformation $A_2=0.5$, and width $\sigma_K=1.0$. This kick is strong enough to perturb the geometry but insufficient to trigger collapse. The initial lapse is set to the pre-collapsed profile $\alpha(t{=}0)=\sqrt{\chi}$ (type~1), and the conformal-factor floor is raised to $\chi_{\rm min}=10^{-4}$ (see below). We use constraint damping $\kappa_1=3.0$ and Kreiss--Oliger dissipation $\sigma_{\rm KO}=1.0$, with $L_{\rm full}=40$, $N_{\rm full}=160$, and 5 levels of AMR ($\text{max\_level}=5$, $dx_{\rm fine}\approx 7.8\times10^{-3}$).

\noindent\textbf{Unperturbed configuration.}
To isolate the natural vacuum fate of the unsupported geometry we also evolve a run with $A_0=A_2=0$ (identically vanishing extrinsic curvature at $t{=}0$) and throat radius $b_0=0.5$. The evolution is therefore driven entirely by the Hamiltonian constraint defect inherent in the unsupported wormhole. Because the absence of any $K$ kick removes the natural gauge-protection mechanism (no $1{+}\log$-driven lapse collapse), two additional measures are required to maintain numerical stability: (i)~the initial lapse is set to the pre-collapsed profile $\alpha(t{=}0)=\sqrt{\chi}$ (type~1), which suppresses the effective evolution speed where the conformal factor is small, and (ii)~the conformal-factor floor is raised to $\chi_{\rm min}=10^{-4}$, creating a protective plateau over the innermost cells near $\bar r=0$ that represent the compactified second asymptotic end rather than physical strong-field geometry. We use $\kappa_1=3.0$, $\sigma_{\rm KO}=1.0$, $L_{\rm full}=40$, $N_{\rm full}=160$, and 5 levels of AMR ($\text{max\_level}=5$, $dx_{\rm fine}\approx 7.8\times10^{-3}$).

\subsection{Vacuum evolution of an unsupported topology}
In classical General Relativity, a traversable wormhole inherently requires exotic matter violating the Null Energy Condition to remain static. Recent no-go theorems confirm that this requirement is a deeply fundamental consequence of the geometric flaring-out condition, persisting even in modified frameworks such as Unimodular Gravity \cite{cataldo2026}. In particular, in the Ellis--Bronnikov geometry the spatial Ricci scalar on the initial hypersurface is negative in the strong-field throat region, \(R^{(3)} < 0\). Consequently, when this geometry is initialized and evolved in pure vacuum (\(T_{\mu\nu}=0\)) without exotic support, the initial slice contains a substantial Hamiltonian constraint defect, and any nontrivial localized \(K\) kick introduces additional momentum-constraint violation. We therefore interpret the unsupported evolution as an out-of-equilibrium vacuum-relaxation problem rather than as a stationary vacuum solution.

Within the CCZ4 formulation, the constraint-damping mechanism rapidly converts the initial defect into propagating transient modes. After this early-time relaxation, the subsequent dynamics exhibit a clear nonlinear threshold separating three regimes: \emph{unperturbed vacuum expansion} (\(A_0{=}A_2{=}0\)), \emph{sub-critical dissipation} (weak kick), and \emph{super-critical dynamical pinch-off} (strong compressive kick).

\subsection{Unperturbed evolution: vacuum expansion of the bare topology}

The most fundamental question about an unsupported wormhole is what happens to the bare geometry when no external perturbation is applied at all. To answer this, we evolve the Ellis--Bronnikov initial data with \(A_0=A_2=0\) (identically vanishing extrinsic curvature at \(t=0\)) and throat radius \(b_0=0.5\). The evolution is therefore driven entirely by the Hamiltonian constraint defect inherent in the unsupported geometry.

\noindent\textbf{Run parameters.} We use the pre-collapsed initial lapse \(\alpha(t{=}0)=\sqrt{\chi}\) (Eq.~\ref{eq:lapse_precollapsed}) and conformal-factor floor \(\chi_{\rm min}=10^{-4}\), with constraint damping \(\kappa_1=3.0\) and Kreiss--Oliger dissipation \(\sigma_{\rm KO}=1.0\). The grid uses \(L_{\rm full}=40\), \(N_{\rm full}=160\), and 5 levels of AMR. Diagnostics are sampled on the finest available level.

\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{diagnostic_unperturbed.pdf}
\caption{Collapse diagnostics for the unperturbed evolution (\(A_0{=}A_2{=}0\), \(b_0=0.5\)). After an early constraint-cleaning transient, \(\min(\alpha)\) recovers monotonically from \(\sim 7\times 10^{-4}\) toward \(\sim 3.5\times 10^{-2}\), and \(\min(\chi)\) rises from the floor \(10^{-4}\) to \(\sim 1.7\times 10^{-3}\). The trapped-surface proxy \(r_{\rm AH}\) starts at \(\sim 0.28\) (reflecting the initial wormhole throat topology) but decreases steadily to \(\sim 0.13\), demonstrating that the residual trapping weakens over time. \(\min(\theta_+)\) remains negative but becomes less so, consistent with a geometry dispersing away from self-trapping. \(\max|K|\) (absolute value) shows a brief spike during constraint cleaning and then settles to moderate values (\(\sim 2.5\)); the spatial snapshots in Fig.~\ref{fig:unperturbed_kz_evolution} confirm that the signed \(K\) is predominantly negative (expansive) in the throat region.}
\label{fig:unperturbed_diagnostics}
\end{figure}

The diagnostics in Fig.~\ref{fig:unperturbed_diagnostics} reveal a geometry that is \emph{expanding}, not collapsing. Despite the initial trapped-surface proxy being nonzero---a consequence of the wormhole throat topology itself---\(r_{\rm AH}\) monotonically decreases throughout the evolution, indicating that the trapping inherited from the initial data is weakening rather than growing. Concurrently, the minimum lapse recovers upward (from \(\sim 7\times 10^{-4}\) to \(\sim 3.5\times 10^{-2}\) by \(t\approx 18\)), confirming the absence of singularity formation.

The spatial structure of the trace of the extrinsic curvature provides the clearest picture of the dynamics. By \(t\approx 15\), the throat region exhibits a persistent \(K<0\) (negative) halo extending to \(r\lesssim 2.5\), while a tiny positive-\(K\) residual remains at the grid origin within the \(\chi\)-floor plateau.

\begin{figure*}[t]
\centering
\includegraphics[width=\textwidth]{evolution_K_z_panel_unperturbed.pdf}
\caption{Evolution snapshots of the \(K\) field in the \(z=0\) plane for the unperturbed run (\(A_0{=}A_2{=}0\)). At early times, the constraint-cleaning transient generates outward-propagating ripples. At late times, the throat region is dominated by \(K<0\) (dark regions), indicating volumetric expansion. The dark spot at the origin is a numerical residual within the \(\chi\)-floor plateau and carries no physical significance.}
\label{fig:unperturbed_kz_evolution}
\end{figure*}

Recall that $K<0$ corresponds to local volumetric expansion of the spatial slices. Physically, this expansion is the geometry's natural response to the removal of the exotic matter that previously balanced the negative spatial Ricci curvature ($R^{(3)}<0$) at the throat. Without the phantom field supplying inward-directed pressure, the bare throat has no mechanism to resist the outward push of its own negative curvature, and the geometry inflates. This is precisely the ``inflationary expansion'' pathway identified by Shinkai and Hayward~\cite{shinkai02}, except here no explicit rarefactive perturbation is required---the unsupported vacuum geometry \emph{naturally} follows this pathway when left unperturbed.

\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{psi_unperturbed.pdf}
\caption{Extracted radius-scaled curvature signal \(r\,\mathrm{Re}(\Psi_4^{2,0})\) for the unperturbed evolution (\(A_0{=}A_2{=}0\)). Because Birkhoff's theorem forbids gravitational radiation from a spherically symmetric spacetime, the unperturbed run should yield \(\Psi_4=0\). The small signal extracted here (amplitude \(\sim 10^{-3}\)) represents the numerical noise floor of our setup, driven by the Cartesian grid and AMR boundaries breaking exact spherical symmetry during the violent constraint-damping phase. This provides a baseline to confirm that the signals in our perturbed runs are physical.}
\label{fig:psi4_unperturbed}
\end{figure}

Because Birkhoff's theorem forbids gravitational radiation from a spherically symmetric spacetime, the unperturbed run should yield \(\Psi_4=0\). The small signal extracted in Fig.~\ref{fig:psi4_unperturbed} (amplitude \(\sim 10^{-3}\)) represents the numerical noise floor of our setup, driven by the Cartesian grid and AMR boundaries breaking exact spherical symmetry during the violent constraint-damping phase. This provides a baseline to confirm that the signals in our perturbed runs are physical. The constraint norms (Appendix, Fig.~\ref{fig:constraints_unperturbed}) remain bounded throughout, confirming that the CCZ4 evolution is well-controlled.

\noindent\textbf{Summary.} An unsupported Ellis--Bronnikov wormhole with no applied perturbation does not collapse. Instead, the bare geometry undergoes volumetric expansion driven by its own negative spatial curvature. The throat opens, the trapped-surface proxy weakens, and the lapse recovers---all signatures of a dispersive, non-collapsing evolution. The missing exotic matter manifests as a Hamiltonian constraint defect that is pushed outward by the CCZ4 constraint-damping mechanism, leaving behind a geometry that slowly relaxes toward a flatter configuration. We emphasize that this specific evolutionary track is mediated by the CCZ4 damping field \(\Theta\): the negative \(R^{(3)}\) acts as a source for \(\Theta\), which in turn drives \(K\) negative. The qualitative outcome---expansion rather than collapse---is therefore an interplay between the geometric properties of the unsupported initial data and the formulation-dependent mechanism by which constraint violations are resolved. Only a sufficiently strong compressive kick can overcome this expansion tendency and force the throat to pinch off into a black hole.

\subsection{Sub-critical branch: dissipation without trapping}
\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{collapse_diagnostics_subcritical.pdf} % Replace with your sub-critical diagnostic plot
\caption{Diagnostics of a sub-critical evolution. The minimum lapse \(\min(\alpha)\) and minimum conformal factor \(\min(\chi)\) do not collapse toward zero, and no trapped surface is detected (\(r_{\rm AH}=0\), equivalently \(\theta_+>0\) everywhere for the spherical proxy). This indicates relaxation/expansion rather than prompt self-trapping.}
\label{fig:subcritical_diagnostics_plot}
\end{figure}
In the absence of a sufficiently strong implosive perturbation, the unsupported geometry follows a dispersive branch. After an initial gauge adjustment, both \(\alpha\) and \(\chi\) relax toward their asymptotic values and the trapped-surface proxy remains zero (\(r_{\rm AH}=0\)), indicating that no spherical trapped region forms on the resolved exterior grid in this branch.

\noindent\textbf{Run parameters.} For the sub-critical evolution shown in Fig.~\ref{fig:subcritical_diagnostics_plot}, we used \(b_0=0.5\) with a modest quadrupolar kick of amplitude \(A_2=0.5\) and width \(\sigma_K=1.0\) (no monopole component). To damp junk radiation and early constraint-cleaning transients, we used comparatively strong damping/dissipation parameters \(\kappa_1=3.0\) and \(\sigma_{\rm KO}=1.0\).

\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{psi_sub_critical.pdf}
\caption{Extracted radius-scaled curvature signal \(r\,\mathrm{Re}(\Psi_4^{2,0})\) for the sub-critical evolution. The retarded-time alignment of the waveforms across the two extraction radii is consistent with an outgoing curvature transient, with an amplitude ($\sim 0.10$) two orders of magnitude above the numerical noise floor established by the unperturbed run (Fig.~\ref{fig:psi4_unperturbed}). Because the initial \(K\) kick does not satisfy the momentum constraint, this signal inevitably contains a contribution from constraint-cleaning modes propagated outward by the CCZ4 evolution; the retarded-time coherence and the large amplitude ratio relative to the noise floor are nevertheless suggestive of a physical component.}
\label{fig:psi_sub_critical}
\end{figure}


\begin{figure*}[t]
\centering
\includegraphics[width=\textwidth]{evolution_panel_sub_critical.pdf}
\caption{Sub-critical evolution snapshots of the \(K\) field in the \(z\)-plane, assembled from selected simulation frames. The grayscale palette highlights the relaxation/dispersion dynamics without the formation of a trapped region.}
\label{fig:subcritical_kz_evolution_panel}
\end{figure*}

\subsection{Super-critical branch: trapping and black hole formation}
In the super-critical regime, the trapped-surface proxy introduced in Sec.~II becomes essential for separating genuine self-trapping from pure gauge behavior. For sufficiently strong inward forcing, \(\theta_+\) crosses zero and the proxy radius \(r_{\rm AH}\) becomes positive, indicating the onset of gravitational trapping and topological pinch-off.

\noindent\textbf{Run parameters.} For the successful collapse shown in Fig.~\ref{fig:supercritical_diagnostics}, we used \(b_0=0.2\) with a strong monopolar compressive kick \(A_0=10.0\), a small quadrupolar deformation \(A_2=0.05\), and width \(\sigma_K=0.2\). In contrast to the sub-critical case, we did not use large damping/dissipation; we set \(\kappa_1=0.1\) and \(\sigma_{\rm KO}=0.5\) to avoid over-damping the nonlinear collapse while still controlling numerical noise.

\begin{figure*}[t]
    \centering
    \includegraphics[width=0.92\textwidth]{collapse_diagnostics_supercritical.pdf}
    \caption{Diagnostics of a super-critical collapse. The minimum lapse (\(\alpha\)) and conformal factor (\(\chi\)) plunge toward zero, demonstrating the time-freezing and infinite spatial stretching characteristic of the moving-puncture gauge. Concurrently, the minimum null expansion proxy \(\min(\theta_+)\) crosses zero and becomes negative, causing the trapped-surface radius proxy \(r_{\rm AH}\) to jump from zero to a finite value. The subsequent step-wise growth in \(r_{\rm AH}\) is an expected artifact of the discrete Cartesian AMR grid tracking the outward coordinate stretching of the trapped region. Together, these features are consistent with dynamical self-trapping and pinch-off in the unsupported experiment.}
    \label{fig:supercritical_diagnostics}
\end{figure*}
When the compressive amplitude exceeds a critical threshold, the implosion overcomes the topological repulsion of the throat. In this super-critical branch, \(\min(\theta_+)\) crosses zero, \(r_{\rm AH}\) becomes positive, and the evolution becomes consistent with black-hole-like self-trapping and throat pinch-off.

\subsection{Gravitational-wave signatures (near-zone)}
The extracted Weyl signal is strongly regime-dependent. In the unperturbed case (\(A_0{=}A_2{=}0\)), as dictated by Birkhoff's theorem, the spherically symmetric setup does not radiate; the small \(\Psi_4\) signal (Fig.~\ref{fig:psi4_unperturbed}) represents the numerical noise floor. In the sub-critical case with a quadrupolar kick, the retarded-time waveform in Fig.~\ref{fig:psi_sub_critical} exhibits an outgoing curvature transient whose amplitude (\(\sim 0.10\)) is two orders of magnitude above the noise floor. In the super-critical case, by contrast, the extracted \(\Psi_4\) is dominated by early-time junk radiation and constraint-cleaning transients. We therefore interpret the super-critical waveform as a near-zone diagnostic of the violent vacuum relaxation rather than as evidence for a clean outgoing gravitational-wave burst.

An important caveat applies to all perturbed runs. The initial extrinsic-curvature kick (Eq.~\ref{eq:kick_profile}) modifies \(K\) without simultaneously solving for a consistent traceless part \(\tilde A_{ij}\), thereby introducing a momentum-constraint violation in addition to the Hamiltonian defect. The CCZ4 evolution damps these violations by propagating them outward at approximately the speed of light, and these constraint-cleaning modes inevitably couple to the Weyl scalar. Consequently, the extracted \(\Psi_4\) contains a superposition of any physical curvature radiation and constraint-driven transients. The retarded-time coherence across extraction radii and the large amplitude ratio relative to the spherically-symmetric noise floor are suggestive of a physical component in the sub-critical signal, but a definitive separation of physical radiation from constraint contamination would require constraint-satisfying initial data---an important target for future work.

Furthermore, all extractions are performed in the near zone (\(R_{\rm ext}=8\text{--}16\) with \(L_{\rm full}=40\)); the signals reported here should therefore be understood as near-zone curvature diagnostics rather than asymptotic gravitational-wave templates.

Note that the extracted waveform amplitude in the sub-critical case is substantially larger than in the super-critical case; this is an expected consequence of seeding the sub-critical run with a much larger initial quadrupolar deformation (\(A_2=0.5\) vs \(0.05\)).
\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{psi_super_critical.pdf}
\caption{Radius-scaled curvature waveform \(r\,\mathrm{Re}(\Psi_4^{2,0})\) for the super-critical run, shown in simulation time \(t\) (top) and retarded time \(u=t-R_{\rm ext}\) (bottom). In this regime the extracted signal is dominated by junk radiation and early transient contamination, so the plot is best interpreted as a near-zone diagnostic of the violent vacuum relaxation rather than as a clean outgoing gravitational-wave template.}
\label{fig:psi_super_critical}
\end{figure}

\section{Discussion: 3D Dynamics versus 1D Instability Proofs}

Our 3D numerical results are consistent with the qualitative bifurcation identified by Shinkai and Hayward \cite{shinkai02}, but the emphasis here is different. Rather than revisiting the 1D instability proof itself, we use full 3D numerical relativity to study how an unsupported wormhole-like geometry relaxes once its sustaining source is absent. Shinkai and Hayward showed that the Ellis--Bronnikov wormhole lies near a bifurcation between contraction and expansion; in their language, a compressive ($K>0$) perturbation drives the throat to collapse into a black hole, while a rarefactive perturbation leads to inflationary expansion. Our calculations revisit that vacuum relaxation problem in a fully 3D setting where non-spherical dynamics and extracted curvature can also be monitored:
\begin{itemize}
    \item \textbf{Unperturbed vacuum expansion (\(A_0{=}A_2{=}0\)):} When no perturbation is applied, the unsupported geometry \emph{naturally expands}. The throat region develops persistent \(K<0\) (volumetric expansion) driven by the negative spatial Ricci curvature that is no longer balanced by exotic matter. The trapped-surface proxy weakens monotonically and the lapse recovers, with no indication of collapse. This establishes that expansion---not stasis---is the default vacuum fate of an unsupported wormhole.
    \item \textbf{Super-critical trapping:} In the super-critical branch, a sufficiently large compressive kick drives the unsupported geometry through a trapping transition in the spherical proxy. The throat pinches off, \(\theta_+\) becomes negative, and the trapped-surface proxy turns on. This is the regime in which the vacuum defect is strong enough to self-trap rather than disperse.
    \item \textbf{Sub-critical relaxation:} In the sub-critical branch, a weak or moderate kick is insufficient to overcome the natural expansion tendency. The negative-curvature region relaxes, the lapse and conformal factor recover, and the resolved exterior spacetime tends toward a flatter weak-curvature configuration. Without a continuously maintained exotic source, the unsupported throat does not remain stationary.
\end{itemize}

Crucially, our 3D approach reveals behavior that 1D codes inherently cannot capture. Because Shinkai and Hayward utilized a 1D spherically symmetric simulation, Birkhoff's theorem prevented any nontrivial curvature signal from being extracted. By moving the problem to full 3D and explicitly breaking exact symmetry with a quadrupolar perturbation, we can inspect \(\Psi_4\) directly. The key lesson of the present simulations is not that every pinch-off yields a clean gravitational-wave burst, but that the near-zone curvature signal depends strongly on the dynamical regime: the sub-critical branch produces the clearest outgoing structure, while the super-critical branch is dominated by junk radiation. We stress that because the initial \(K\) kick does not satisfy the momentum constraint, CCZ4 constraint-cleaning modes contribute to the extracted Weyl scalar; disentangling physical radiation from this contamination is an important direction for future work with constraint-satisfying initial data.

\section{Conclusion}
We have investigated the nonlinear vacuum dynamics of unsupported wormhole-like topological defects using full 3D numerical relativity. The central question in this work is not how to sustain a traversable wormhole with exotic matter, but what becomes of the residual geometry once that support is absent and the object is perturbed.

The evolution exhibits a clear three-regime structure. With no applied perturbation, the bare geometry naturally expands: the throat's negative spatial curvature, no longer balanced by exotic matter, drives volumetric expansion (\(K<0\)) and the trapped-surface proxy weakens monotonically. Sub-critical perturbations produce a similar dispersive outcome, while sufficiently strong compressive kicks produce a super-critical branch in which the throat pinches off and the trapped-surface proxy becomes positive. In this sense, unsupported wormholes in vacuum CCZ4 behave as dynamical topological defects whose default fate is expansion, with a threshold to self-trapping that can only be crossed by a strong compressive perturbation. In particular, the unsupported branch studied here does not remain as a static traversable throat on the resolved grid. While quantum backreaction may theoretically stabilize microscopic topological defects at the Planck scale \cite{mehulic2026}, our results demonstrate that once stretched to macroscopic scales where classical vacuum dynamics dominate, they lack the sustained exotic support required to survive. This provides a natural cosmological mechanism explaining why inflationary wormholes do not abundantly persist in the late universe.

The extracted near-zone curvature signal is likewise regime-dependent. In the sub-critical regime, retarded-time alignment of the waveforms across multiple extraction radii reveals an outgoing curvature transient whose amplitude exceeds the numerical noise floor by two orders of magnitude. However, because the ad-hoc extrinsic-curvature kick introduces both Hamiltonian and momentum constraint violations, the CCZ4 constraint-cleaning modes inevitably contaminate the Weyl scalar; a definitive claim of physical gravitational radiation requires constraint-satisfying initial data, which we leave to future work. In the super-critical regime, the extracted \(\Psi_4\) is dominated by junk radiation and early transient contamination. These results show that understanding unsupported wormhole dynamics requires tracking the geometry and its diagnostics together: the fate of the topological defect is robustly encoded in the collapse indicators, while the interpretation of \(\Psi_4\) depends sensitively on the regime and on the contamination introduced by the initial constraint defect.

\section{Acknowledgements}
This research was supported by Gravity Frontiers (\texttt{ic@gravityfrontiers.org}) under a research grant for the development of software for mathematical modeling of wormholes.

The author also thanks Ilya Nachevsky for generously providing the high-performance computing resources necessary for this work. These simulations would not have been possible without his support. 

\appendix*
\section{Numerical details, validation, and convergence}
We exploit Cartesian reflection symmetries and evolve only the positive octant, \(x\ge 0\), \(y\ge 0\), \(z\ge 0\), imposing parity conditions on the inner faces and Sommerfeld radiative conditions at the outer boundary. For the validation diagnostics shown here, the constraints are computed directly during the evolution on AMR level 0, after filling two ghost cells of the updated state. The code then evaluates the standard CCZ4 constraint operator into a temporary four-component field containing the Hamiltonian constraint and the three momentum-constraint components,
\[
(\mathrm{Ham},\mathrm{Mom}_x,\mathrm{Mom}_y,\mathrm{Mom}_z),
\]
and forms volume-weighted root-mean-square norms over the valid level-0 cells. Since the cell volume is constant on a given AMR level, the reported diagnostics reduce to
\begin{align}
L_2(\mathrm{Ham}) &=
\left(\frac{\sum_i \mathrm{Ham}_i^2\,V_i}{\sum_i V_i}\right)^{1/2},\\
L_2(\mathrm{Mom}) &=
\left(\frac{\sum_i \left(\mathrm{Mom}_{x,i}^2+\mathrm{Mom}_{y,i}^2+\mathrm{Mom}_{z,i}^2\right)V_i}{\sum_i V_i}\right)^{1/2}.
\end{align}
These two scalars are written to the small-data file as \texttt{L2\_Ham} and \texttt{L2\_Mom}. Because the unsupported or partially unsupported evolutions begin from data that are not exactly on the target vacuum constraint surface, these norms are generically nonzero at \(t=0\). The relevant stability test is therefore not whether the initial slice is constraint-free, but whether the CCZ4 evolution keeps the subsequent violation bounded. Figure~\ref{fig:constraints_plot} shows that both the Hamiltonian and momentum RMS norms exhibit an initial transient and then remain bounded throughout the evolution. We note that the saturated Hamiltonian norm (\(\sim 2.5\times 10^{-2}\)) is large by the standards of constraint-satisfying numerical relativity, reflecting the fact that the grid never fully settles onto the Einstein constraint surface---an expected consequence of the deliberately unphysical initial data and the persistent under-resolution of the compactified second asymptotic end near $\bar r=0$. The boundedness of the norms nevertheless confirms numerical stability and ensures that the CCZ4 evolution does not develop runaway violations that would invalidate the dynamical conclusions.
\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{constraints_plot.pdf}
\caption{Level-0, volume-weighted RMS norms of the Hamiltonian and momentum constraints written to \texttt{constraint\_norms.dat} as \texttt{L2\_Ham} and \texttt{L2\_Mom}. The initial slice carries a finite constraint defect, but the CCZ4 evolution damps the transient and keeps both diagnostics bounded during the subsequent evolution.}
\label{fig:constraints_plot}
\end{figure}

\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{constraints_unperturbed.pdf}
\caption{Volume-weighted RMS norms of the Hamiltonian and momentum constraints for the unperturbed evolution. The Hamiltonian norm grows slowly but remains at \(\mathcal{O}(10^{-2})\); the momentum norm saturates at \(\sim 8\times 10^{-3}\). These levels are large relative to constraint-satisfying evolutions, reflecting the deliberately unphysical initial data, but the norms remain bounded throughout the \(t\approx 18\) evolution, confirming numerical stability.}
\label{fig:constraints_unperturbed}
\end{figure}
% --- Bibliography ---
\clearpage
\bibliographystyle{apsrev4-2}
\bibliography{references}
\end{document}