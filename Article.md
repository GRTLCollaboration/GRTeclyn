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
\title{\MakeTextUppercase{Closing the Bridge: Gravitational signatures of wormhole collapse}}

\author{\MakeTextUppercase{Nikita M. Shirokov}}
\email{shirokov.nm@phystech.edu}
\affiliation{\textit{Independent Researcher, Moscow, Russia}}

% --- Version Date ---
\date{Version July 10, 2025}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- CORRECTED POSITION FOR THE ABSTRACT ---
% For REVTeX, the abstract MUST come BEFORE the \maketitle command.
% The class will then automatically format it correctly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\begin{abstract}
Wormholes are solutions to the Einstein field equations of General Relativity that represent topological bridges connecting distinct regions of spacetime. Extensive research has investigated stable wormhole metrics, which typically necessitate exotic matter violating the null energy condition. Consequently, unstable solutions devoid of exotic matter have received considerably less attention, despite the profound insights they offer into dynamic topological transitions. In this work, we investigate the fundamental nonlinear dynamics and near-field curvature emission of a collapsing traversable wormhole. By employing 3D numerical relativity simulations, we demonstrate that the transition from a topological bridge to a black hole emits a distinct gravitational signature, providing valuable theoretical insights into the dynamics of exotic compact objects.
\end{abstract}


% The \maketitle command generates the title block AFTER the abstract
\maketitle

% --- Main Content of the Proposal ---
\section{Introduction}
Wormholes—topological bridges connecting distinct regions of spacetime or separate universes—stand among the most fascinating predictions of General Relativity. The geometric possibility of such structures was first identified by Ludwig Flamm in 1916 \cite{flamm16}, shortly after the discovery of the Schwarzschild solution. The concept was formally introduced into physics by Einstein and Rosen in 1935 as an attempt to model elementary particles as non-singular field configurations \cite{einstein35}. However, it was not until 1957 that John Archibald Wheeler coined the term "wormhole" while investigating the topological foam of spacetime at the quantum scale \cite{misner57}.

For decades, these solutions were considered topological curiosities, primarily because the standard Einstein-Rosen bridge is non-traversable; the throat pinches off faster than even a photon can cross it. The modern era of wormhole physics began in 1988 when Morris and Thorne demonstrated that static, traversable wormholes are valid solutions to the Einstein equations, provided the throat is threaded by "exotic matter" that violates the null energy condition \cite{morris88}. Consequently, a significant portion of the literature has since focused on minimizing this exotic matter requirement or exploring modified theories of gravity to sustain stable wormholes without it \cite{visser95, lobo05}.

The motivation for the present work, however, diverges entirely from the search for stability. As envisaged by Wheeler's concept of spacetime foam, transient quantum wormholes might constantly form and vanish at the Planck scale. It is theorized that the early universe's rapid exponential expansion (inflation) could have provided a mechanism for enlarging these primordial topological defects to macroscopic scales \cite{kardashev07, garcia16}. Unlike the carefully fine-tuned metrics required for interstellar travel, these primordial structures would likely be inherently unstable. Indeed, recent theoretical work by \citet{dimaschko25} suggests that any bound state of a traversable wormhole is non-stationary. Rather than investigating the conditions required to stabilize them, this work embraces their inevitable dynamical collapse to search for observational signatures.

The dynamical fate of such traversable wormholes was first established in a pioneering numerical study by Shinkai and Hayward \cite{shinkai02}. Evolving the classic Morris-Thorne (Ellis-Bronnikov) metric in 1D spherical symmetry, they demonstrated that these structures are highly unstable against perturbations. They observed that the wormhole throat suffers a bifurcation of horizons, leading either to an inflationary expansion (given negative energy) or a rapid collapse into a black hole (given positive energy). Crucially, the introduction of even a minute pulse of normal matter—representing, for instance, a traveler attempting to traverse the bridge—inevitably forces the wormhole to collapse into a black hole, sealing off causal connection between the two regions. Their work elegantly proved that exotic topological structures can be rigorously studied using the same local, dynamical methods applied to standard black holes, effectively unifying the two phenomena. 

While Shinkai and Hayward established the ultimate fate of the wormhole topology, their strict assumption of 1D spherical symmetry inherently precluded the study of its gravitational radiation. According to Birkhoff's theorem, a purely spherically symmetric collapse cannot emit gravitational waves. In reality, the violent topological transition from a traversable bridge into disjoint black holes should theoretically produce a powerful, asymmetric gravitational wave (GW) outburst. To capture these critical dynamics, one must move beyond 1D restrictions and perform full 3D, non-linear dynamical simulations.

Simulating exotic spacetimes without the continuous support of exotic matter presents a massive computational challenge. However, a methodological blueprint for such endeavors was recently established by \citet{clough24}, who successfully simulated the 3D dynamical collapse of an Alcubierre warp drive bubble. Following a similar philosophy, we overcome the limitations of early 1D models by deploying state-of-the-art 3+1 Numerical Relativity (NR) solvers. Specifically, this project utilizes the \texttt{GRTeclyn} framework—built upon the highly successful \texttt{GRChombo} infrastructure \cite{Andrade2021}—to robustly handle the complex, multi-puncture topology of the wormhole on an adaptive 3D Cartesian grid.

The ultimate goal of this research is to investigate the fundamental nonlinear dynamics and near-field curvature emission of a topological collapse. The advent of GW astronomy by the LIGO-Virgo-KAGRA collaborations \citep{abbott16,abbott18,abbott21,abbott23} has moved the search for Exotic Compact Objects (ECOs) \cite{cardoso19} from pure theory to observational reality. However, understanding the dynamical signatures of ECOs requires rigorous theoretical modeling. Because our computational domain restricts extraction to the strong-field regime, we do not aim to provide asymptotic observational templates for GW observatories. Instead, by extracting the near-zone waveforms from a full 3D collapse, we aim to identify the unique near-field curvature dynamics that distinguish this topological transition from standard black hole formation. This remains a highly valuable theoretical result for the numerical relativity community and a critical step toward understanding the dynamical signatures of macroscopic topological defects.

\section{Formalism and Methodology}

\subsection{3+1 Decomposition and the CCZ4 Formulation}
We follow the standard numerical-relativity approach to solving the Einstein Equations,
\begin{equation}
G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2} R g_{\mu\nu} = 8\pi T_{\mu\nu},
\end{equation}
as an initial value problem, using the 3+1 decomposition and setting $G = c = 1$. The four-dimensional spacetime manifold is foliated into a family of spacelike hypersurfaces $\Sigma_t$, parametrized by a time coordinate $t$. The line element is written in terms of the lapse function $\alpha$, the shift vector $\beta^i$, and the spatial metric $\gamma_{ij}$ as:
\begin{equation}
    ds^2 = \left(-\alpha^2 + \beta_i \beta^i \right) dt^2 + 2 \beta_i dt dx^i + \gamma_{ij} dx^i dx^j.
\end{equation}
For the supported-collapse experiments studied here, the dynamical evolution is governed by the Einstein equations coupled to a phantom scalar field, implemented within the \texttt{GRTeclyn} framework using the conformal and covariant Z4 (CCZ4) formulation. Our objective is to start from a constraint-consistent Ellis-Bronnikov wormhole at $t=0$ and then drive the system away from equilibrium by reducing the exotic support that keeps the throat open. In this setup, the primary trigger is the controlled removal of the supporting stress-energy, not a large vacuum kick. Any extrinsic-curvature perturbation that is added is kept parametrically small and is used only to select the unstable branch and break exact spherical symmetry. This ``clean'' supported setup is therefore conceptually distinct from separate forced-collapse vacuum experiments, in which a strong kick dominates the early-time dynamics.

\subsection{Initial Data: The Morris-Thorne Wormhole in Isotropic Coordinates}

We initialize the simulation using the static, spherically symmetric Morris-Thorne metric:
\begin{equation}
ds^2 = -e^{2\Phi(r)} dt^2 + \frac{dr^2}{1 - b(r)/r} + r^2 (d\theta^2 + \sin^2\theta d\phi^2).
\end{equation}
To balance physical relevance with numerical stability, we adopt the specific functional forms corresponding to the Ellis-Bronnikov wormhole, a canonical massless scalar field solution. We set the redshift function to zero ($\Phi(r) = 0$), which ensures the absence of initial horizons, sets the proper time of Eulerian observers equal to coordinate time at $t=0$, and eliminates radial tidal forces at the throat. We choose the inverse-$r$ shape function $b(r) = b_0^2/r$, where $b_0$ represents the throat radius. This choice guarantees asymptotic flatness as $r \to \infty$.

In standard Schwarzschild-like coordinates, the metric component $g_{rr} = (1 - b_0^2/r^2)^{-1}$ diverges at the throat ($r=b_0$). To resolve this coordinate singularity and simultaneously map both "universes" connected by the bridge onto a single computational grid, we transform to \textbf{isotropic coordinates}. We introduce the isotropic radius $\bar{r}$, related to the proper distance $\ell = \pm \sqrt{r^2 - b_0^2}$ via the embedding $\ell = b_0 \sinh(\ln(2\bar{r}/b_0))$. 

This transforms the spatial line element into a manifestly conformally flat form:
\begin{equation}
    ds^2_{\rm spatial} = \psi^4 \left( d\bar{r}^2 + \bar{r}^2 d\Omega^2 \right) = \psi^4 \delta_{ij} dx^i dx^j,
\end{equation}
where the conformal factor $\psi$ encodes the non-trivial topology:
\begin{equation}
    \psi = 1 + \frac{b_0^2}{4\bar{r}^2}.
\end{equation}
In this chart, the primary universe extends to $\bar{r} \to \infty$ ($\psi \to 1$), the throat is located at $\bar{r} = b_0/2$, and the secondary universe is compactified such that its asymptotic infinity maps to the grid origin $\bar{r} \to 0$ ($\psi \to \infty$). 

\subsection{Transformation to CCZ4 Variables and Regularization}

For stable numerical evolution in the CCZ4 formulation, we utilize the conformal transverse-traceless decomposition. The physical spatial metric $\gamma_{ij}$ is represented by a conformal metric $\tilde{\gamma}_{ij}$ with unit determinant, and a regularized conformal factor $\chi$:
\begin{align}
    \tilde{\gamma}_{ij} &= \psi^{-4} \gamma_{ij} = \delta_{ij}, \\
    \chi &= \psi^{-4} = \left( 1 + \frac{b_0^2}{4\bar{r}^2} \right)^{-4}.
\end{align}
To strictly avoid division-by-zero errors at the grid origin ($\bar{r}=0$), the conformal factor is algebraically rewritten and implemented as:
\begin{equation}
    \chi = \left( \frac{4\bar{r}^2}{4\bar{r}^2 + b_0^2} \right)^{4}.
\end{equation}
This smoothly vanishes at the origin ($\chi \to 0$), accurately mirroring the "puncture" methodology used for black holes. The initial data passed to the solver therefore consists of a flat conformal metric ($\tilde{\gamma}_{ij} = \delta_{ij}$) and a coupled scalar field $\phi$ profile that entirely captures the wormhole geometry and stress-energy.

\subsection{Dynamical Triggering via Support Removal}

The initialized Ellis-Bronnikov geometry requires exotic matter violating the Null Energy Condition to remain static. Previous studies often relied on a ``sudden approximation,'' evolving this geometry in a pure vacuum ($T_{\mu\nu}=0$), which results in a massive, instantaneous violation of the Hamiltonian constraint and violent numerical transients.

To achieve a physically consistent, controlled collapse, we instead explicitly model the required exotic matter. We couple the CCZ4 evolution to a "phantom" scalar field $\phi$. The stress-energy tensor for this field is defined with a reversed overall sign relative to a canonical scalar field, providing the negative energy density necessary to support the wormhole throat and satisfy the constraints at $t=0$.

We systematically control the onset of the dynamics by reducing this support in a causal, throat-centered manner. Rather than switching the entire slice at once, we define a local support multiplier
\begin{equation}
    S(\mathbf{x},t) = S_0\!\left(t - \frac{r_{\rm th}(\mathbf{x})}{v_{\rm c}}\right),
\end{equation}
where $r_{\rm th}(\mathbf{x})$ measures distance from the \emph{throat surface} and $v_{\rm c}$ is the speed at which the support-removal front propagates outward. In the single-throat isotropic Ellis-Bronnikov chart used here, the coordinate origin $\bar r=0$ is not the throat itself but the compactified asymptotic end of the second universe, while the physical throat is the spherical shell at $\bar r=b_0/2$ around the grid center. Accordingly, the support-removal front is launched from that shell rather than from the coordinate origin. The throat-center schedule $S_0$ is held at unity during an initial settling stage and then reduced smoothly to zero over a duration $\Delta t$ using a cosine profile:
\begin{equation}
    S_0(t) =
    \begin{cases}
        1, & t \le t_{\rm start}, \\
        \frac{1}{2} \left[ 1 + \cos\left( \pi \frac{t - t_{\rm start}}{\Delta t} \right) \right], & t_{\rm start} < t < t_{\rm start} + \Delta t, \\
        0, & t \ge t_{\rm start} + \Delta t.
    \end{cases}
\end{equation}
This construction suppresses the unphysical, non-causal whole-domain response produced by a global quench and ensures that the unsupported region first forms near the throat and only later propagates outward through the rest of the slice. In the clean supported-collapse setup, we set the initial-data kick to zero at $t=0$ so that the spacetime begins in a supported equilibrium state without an artificial implosive impulse. After the support-removal front has passed through the throat and the geometry has entered the unsupported regime, we then apply a separate weak delayed pulse in the trace of the extrinsic curvature,
\begin{equation}
    K(\bar{r}, \theta, \phi) = \left[A_0 + A_2 Y_{20}(\theta, \phi)\right]
    \exp\left(-\frac{\bar{r}^2}{\sigma^2}\right),
\end{equation}
with $|A_0|, |A_2| \ll 1$, multiplied by a compact-in-time envelope and activated only after a prescribed delay $t_{\rm kick} > t_{\rm start} + \Delta t$. Operationally, this means that support removal remains the primary physical trigger, while the delayed $K$ pulse acts only as a branch selector and symmetry breaker that nudges the now-unsupported throat off the near-equilibrium branch. This delayed forcing cleanly separates the physical de-supporting process from the later perturbation used to probe whether the unsupported geometry relaxes, disperses, or enters genuine topological collapse.

\subsection{Apparent Horizon Detection}

To definitively distinguish true topological pinch-off and black hole formation from mere gauge effects (such as "lapse collapse"), we implemented a geometric trapped-surface diagnostic. At each time step, we perform a spherical-shell scan outward from the origin, calculating the expansion of outgoing null rays, $\theta_+$. In spherical symmetry, this is given by:
\begin{equation}
    \theta_+ = \frac{2\sqrt{\chi}}{r} + \frac{2}{3}K - \chi \tilde{A}_{rr}.
\end{equation}
The formation of a trapped surface is unambiguously identified when $\theta_+ \le 0$. The maximum radius where this condition holds provides a proxy for the apparent horizon radius, $r_{\rm AH}$. The sudden jump of $r_{\rm AH}$ from zero to a finite positive value during the simulation serves as our primary physical indicator of successful collapse.

\subsection{Gauge Conditions}

To ensure stability through the highly non-linear topological transition from a traversable wormhole into disjoint black holes, we employ robust gauge choices. 

For the shift vector, we initially set $\beta^i = 0$. Since our initial state possesses no tangential rotation or shift, and the physical mechanism under investigation is a radial pinch-off, a vanishing initial shift prevents gauge-induced artifacts from contaminating the early gravitational wave extraction.

For the lapse function, we employ the $1+\log$ slicing condition, which evolves as $\partial_t \alpha - \beta^i \partial_i \alpha = -2\alpha K$. In the clean supported runs we initialize the lapse as
\begin{equation}
    \alpha(t=0) = 1,
\end{equation}
so that the subsequent lapse evolution is generated dynamically by the matter removal and the geometry itself rather than being biased toward an already ``pre-collapsed'' gauge profile. During the evolution, if the throat truly enters a collapsing branch and physical singularities begin to form, the $1+\log$ condition drives the lapse toward zero (``lapse collapse''), effectively freezing the evolution near the singular region and preventing numerical breakdown.





\section{Numerical Setup}

Numerical evolutions were performed using the \texttt{GRTeclyn} codebase, which is currently under active development. Built upon the \texttt{AMReX} framework~\cite{amrex}, \texttt{GRTeclyn} implements block-structured Adaptive Mesh Refinement (AMR) and is explicitly optimized for high-performance GPU acceleration. Production runs targeted H100-class GPUs (corresponding to \texttt{CUDA\_ARCH=90} in the build configuration) and used one MPI rank per GPU, scaling up to 8 GPUs for the largest calculations. In practice, the exact runtime configuration depends on the local MPI stack: when CUDA-aware MPI was unavailable or unstable, we disabled GPU-aware MPI communication (\texttt{amrex.use\_gpu\_aware\_mpi = 0}) while retaining the same distributed multi-GPU layout.

The full physical computational domain spans a coordinate length of $L_{\text{full}} = 80.0$, covered by a coarse grid of $N_{\text{full}} = 160$ cells, yielding a coarse resolution of $dx_{\text{coarse}} = 0.5$. However, to drastically minimize memory consumption, we rigorously exploit the Cartesian reflection symmetries inherent in the spherically symmetric initial geometry. We evolve only the positive octant ($x \ge 0, y \ge 0, z \ge 0$) of the full domain, effectively modeling just 1/8th of the physical volume. Parity symmetry conditions are strictly applied at the inner reflection boundaries ($x=0, y=0, z=0$), while Sommerfeld radiative boundary conditions are enforced at the outer edges.

To capture the extreme curvature gradients developing during the collapse while maintaining tractability, we employ 5 levels of 2:1 adaptive mesh refinement (AMR). The mesh is dynamically regridded based on the gradients of the conformal factor ($\chi$), allowing the code to focus resolution precisely where the wormhole throat is pinching off. On the finest level of refinement (level 5), the grid spacing reaches $dx_{\text{fine}} \approx 0.0156$. Time stepping is handled using a 4th-order Runge-Kutta scheme, with a conservative Courant factor (dt\_multiplier) reduced to 0.05 to maintain stability during the long, weakly forced supported evolutions as well as the more abrupt collapse branches.

The primary objective of this work is to investigate the collapse dynamics subject to Gaussian extrinsic curvature perturbations. We explored two main geometric configurations:
\begin{enumerate}
    \item \textbf{Small Throat ($b_0 = 0.5$):} The baseline simulations successfully utilized a wormhole throat radius of 0.5 with a perturbation width of $\sigma = 1.0$. This configuration, even with 5 levels of AMR tracking the implosion, remained within the memory limits of our 8-GPU array, allowing us to evolve the spacetime through complete horizon formation and extract the entire post-merger ringdown waveform up to $t=30M$.
    
    \item \textbf{Large Throat ($b_0 = 1.0$):} To test geometric scaling, the throat radius was doubled to 1.0. As previously discussed, triggering collapse in this larger geometry required doubling the perturbation width to $\sigma = 2.0$. This volumetrically massive perturbation injected significantly more kinetic energy into the grid, resulting in a much more violent, spatially extended dynamical response. The AMR algorithm reacted by aggressively expanding the high-resolution grids (levels 3, 4, and 5) to track the massive strong-field gradients propagating outward. Ultimately, the memory footprint of these massively expanded refinement levels exceeded the capacity of the 8-GPU cluster, resulting in out-of-memory (OOM) failures before the final quasi-normal ringdown could be entirely resolved.
\end{enumerate}


\subsection{Initial data problem}
The existence of the constraints in general relativity implies that it is not possible to choose arbitrarily the 12 dynamical quantities $\{\gamma_{ij} , K_{ij}\}$ as initial data. The initial data must be chosen in such a way that the four constraint equations are satisfied on the initial spatial hypersurface $\Sigma_0$. 

This means that before starting a numerical evolution, it is necessary to solve the initial data problem to obtain adequate values of $\{\gamma_{ij} , K_{ij}\}$ that represent the specific physical situation under study, such as a traversable wormhole or a binary system of compact objects.

In the 3+1 formalism, the constraints for a vacuum spacetime ($T_{\mu\nu}=0$) consist of the Hamiltonian constraint:
\begin{equation}
R + K^2 - K_{ij}K^{ij} = 0,
\end{equation}
and the momentum constraints:
\begin{equation}
D_j (K^{ij} - \gamma^{ij} K) = 0,
\end{equation}
where $R$ is the Ricci scalar of the spatial metric $\gamma_{ij}$, $D_j$ is the covariant derivative associated with $\gamma_{ij}$, and $K_{ij}$ is the extrinsic curvature. 

To find solutions for the wormhole and black hole geometries discussed in this work, we utilize the \textit{conformal decomposition} approach. We assume the spatial metric is conformally flat ($\gamma_{ij} = \psi^4 \delta_{ij}$) and that the system is at a moment of time symmetry ($K_{ij} = 0$). Under these assumptions, the momentum constraints are identically satisfied, and the Hamiltonian constraint reduces to the Laplace equation for the conformal factor:
\begin{equation}
\nabla^2 \psi = 0.
\end{equation}
The specific topology of the spacetime (e.g., number of mouths or punctures) is then determined by the boundary conditions and the choice of the singular parts of the solution $\psi$.

\subsubsection*{Unsupported vacuum data should be interpreted as a constraint defect, not as a stationary wormhole}
The Ellis--Bronnikov conformal factor
\begin{equation}
    \psi = 1 + \frac{b_0^2}{4\bar r^2}
\end{equation}
describes a traversable throat only when it is accompanied by the stress-energy source that violates the Null Energy Condition. This is not merely a matter of taste: the classical singularity theorems of Penrose and Hawking are formulated precisely around such energy conditions, so a traversable wormhole unsupported by exotic matter is generically not expected to remain a regular stationary vacuum object. If one takes the same geometry and evolves it with $T_{\mu\nu}=0$, the vacuum Hamiltonian constraint
\begin{equation}
    R + K^2 - K_{ij}K^{ij} = 0
\end{equation}
is no longer satisfied by the supported wormhole data on the initial slice. In other words, the missing exotic support appears numerically as an initial Hamiltonian-constraint defect.

This observation provides the correct interpretation of the unsupported \texttt{WormholeCollapse} experiments. They should not be viewed as equilibrated astrophysical wormholes evolving in clean vacuum, but rather as macroscopic topological defects placed out of equilibrium and then allowed to relax under CCZ4. In that setting, the sensitivity to damping parameters, expansion versus collapse, and the presence or absence of an initial inward $K$ pulse are not simply coding pathologies. They are signatures of a genuine bifurcation in how the Einstein system resolves the unsupported geometry back toward the vacuum constraint surface.

Operationally, two branches are observed. With weak or no inward forcing, the code damps and propagates away the initial defect, and the throat relaxes by widening or dispersing instead of self-trapping. With a sufficiently strong localized negative-$K$ perturbation, the same unsupported geometry can be driven onto a supercritical implosive branch that forms an apparent horizon. We therefore interpret the imposed $K$ profile not as a claim that the unsupported wormhole is naturally stable until ``numerically kicked,'' but as a controlled momentum fluctuation used to test whether the defect disperses or collapses once its exotic support has been removed.

\section{Results}

\subsection{Vacuum evolution of an unsupported topology}
In classical General Relativity, a traversable wormhole requires exotic matter violating the Null Energy Condition to remain static. In particular, in the Ellis--Bronnikov geometry the spatial Ricci scalar on the initial hypersurface is negative in the strong-field throat region, \(R^{(3)} < 0\). Consequently, when this geometry is initialized and evolved in pure vacuum (\(T_{\mu\nu}=0\)) without exotic support, the initial slice contains a substantial Hamiltonian constraint defect. We interpret this defect physically as an out-of-equilibrium macroscopic topological bridge rather than as a stationary vacuum solution.

Within the CCZ4 formulation, the constraint-damping mechanism rapidly converts the initial defect into propagating transient modes. After this early-time relaxation, the subsequent dynamics exhibit a clear nonlinear threshold separating two regimes: \emph{topological dissipation} (sub-critical evolutions) and \emph{dynamical pinch-off} (super-critical collapse).

\subsection{Sub-critical branch: dissipation without trapping}
\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{collapse_diagnostics_plot.pdf} % Replace with your sub-critical diagnostic plot
\caption{Diagnostics of a sub-critical evolution. The minimum lapse \(\min(\alpha)\) and minimum conformal factor \(\min(\chi)\) do not collapse toward zero, and no trapped surface is detected (\(r_{\rm AH}=0\), equivalently \(\theta_+>0\) everywhere for the spherical proxy). This indicates dissipation/expansion rather than black hole formation.}
\label{fig:subcritical_diagnostics_plot}
\end{figure}
In the absence of a sufficiently strong implosive perturbation, the unsupported geometry follows a dispersive branch. After an initial gauge adjustment, both \(\alpha\) and \(\chi\) relax toward their asymptotic values and the trapped-surface proxy remains zero (\(r_{\rm AH}=0\)), confirming that no event horizon forms in this branch.

\subsection{Super-critical branch: trapping and black hole formation}
Because a full elliptic Apparent Horizon finder is currently under active development within the \texttt{GRTeclyn} framework, we implemented a local trapped-surface proxy to track the dynamical pinch-off of the wormhole. At each time step, we evaluate the outgoing null expansion \(\theta_+ = \frac{2\sqrt{\chi}}{r} + \frac{2}{3}K - \chi \tilde{A}_{rr}\) assuming spherical symmetry. Grid cells satisfying \(\theta_+ \le 0\) are flagged as trapped. We report the maximum coordinate radius of these flagged cells as \(r_{\rm AH}\). While this proxy does not measure the exact physical area of an apparent horizon, its sudden transition from zero to a positive value serves as a robust, unambiguous indicator of the onset of gravitational trapping and the formation of a black hole. The subsequent step-wise growth of \(r_{\rm AH}\) reflects the discrete sampling of the AMR grid as the interior moving-puncture gauge coordinates stretch.

\begin{figure}[h]
    \centering
    \includegraphics[width=\linewidth]{collapse_diagnostics_plot.pdf}
    \caption{Diagnostics of a super-critical collapse. The minimum lapse (\(\alpha\)) and conformal factor (\(\chi\)) plunge toward zero, demonstrating the time-freezing and infinite spatial stretching characteristic of the moving-puncture gauge. Concurrently, the minimum null expansion proxy \(\min(\theta_+)\) crosses zero and becomes negative, causing the trapped-surface radius proxy \(r_{\rm AH}\) to jump from zero to a finite value. The subsequent step-wise growth in \(r_{\rm AH}\) is an expected artifact of the discrete Cartesian AMR grid tracking the outward coordinate stretching of the trapped region. Together, these features definitively mark dynamical pinch-off and black hole formation.}
    \label{fig:supercritical_diagnostics}
\end{figure}
To model a strong external perturbation (e.g. an early-universe fluctuation) we seed a localized inward perturbation in the trace of the extrinsic curvature. Because a perfectly spherically symmetric vacuum collapse cannot radiate (Birkhoff's theorem), we allow an optional small quadrupolar deformation to break exact symmetry:
\begin{equation}
K(\bar r,\theta,\phi) =
\left[A_0 + A_2\,Y_{20}(\theta,\phi)\right]\exp\!\left(-\frac{\bar r^2}{\sigma^2}\right).
\end{equation}
When the compressive amplitude exceeds a critical threshold, the implosion overcomes the topological repulsion of the throat. In this super-critical branch, \(\min(\theta_+)\) crosses zero, \(r_{\rm AH}\) becomes positive, and the throat pinches off into a black hole.

\subsection{Gravitational-wave signatures (near-zone)}
The dynamical pinch-off produces a sharp outgoing curvature burst. Because the initial vacuum slice contains a constraint defect, the earliest signal is contaminated by gauge and constraint-cleaning transients. We therefore analyze the extracted Weyl scalar \(\Psi_4\) in terms of retarded time \(u \equiv t - R_{\rm ext}\) to separate early junk from the physical collapse burst.
\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{psi4_extracted_simulation.pdf}
\caption{Radius-scaled curvature waveform \(r\,\mathrm{Re}(\Psi_4^{2,0})\) for a super-critical collapse case, shown in simulation time \(t\) (top) and retarded time \(u=t-R_{\rm ext}\) (bottom). Alignment of the dominant packet in retarded time across radii supports its interpretation as an outward-propagating physical signal.}
\label{fig:psi4_extracted_simulation}
\end{figure}

\section{Conclusion}
We investigated the nonlinear dynamics of an unsupported wormhole-like topology evolved in vacuum CCZ4. The evolution exhibits threshold behavior: sub-critical configurations dissipate/expand without trapping, whereas sufficiently strong inward perturbations drive the throat into dynamical pinch-off and black hole formation, unambiguously identified by \(\theta_+ \le 0\) and \(r_{\rm AH}>0\) in the trapped-surface proxy. In the super-critical branch, the extracted \(\Psi_4\) shows a distinct near-zone curvature burst followed by ringdown consistent with a perturbed black hole.

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
These two scalars are written to the small-data file as \texttt{L2\_Ham} and \texttt{L2\_Mom}. Because the unsupported or partially unsupported evolutions begin from data that are not exactly on the target vacuum constraint surface, these norms are generically nonzero at \(t=0\). The relevant stability test is therefore not whether the initial slice is constraint-free, but whether the CCZ4 evolution damps the initial defect and keeps the subsequent violation bounded during the collapse. Figure~\ref{fig:constraints_plot} shows precisely this behavior: both the Hamiltonian and momentum RMS norms exhibit an initial transient and then remain under control throughout the evolution.
\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{constraints_plot.pdf}
\caption{Level-0, volume-weighted RMS norms of the Hamiltonian and momentum constraints written to \texttt{constraint\_norms.dat} as \texttt{L2\_Ham} and \texttt{L2\_Mom}. The initial slice carries a finite constraint defect, but the CCZ4 evolution damps the transient and keeps both diagnostics bounded during the subsequent evolution.}
\label{fig:constraints_plot}
\end{figure}
% --- Bibliography ---
\bibliographystyle{apsrev4-2}
\bibliography{references}
\end{document}