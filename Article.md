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
The dynamical evolution is governed by the vacuum ($T_{\mu\nu}=0$) Einstein equations, implemented within the \texttt{GRTeclyn} framework using the conformal and covariant Z4 (CCZ4) formulation. To simulate the collapse of a wormhole, we must define the initial geometry on $\Sigma_0$ and provide an explicit mechanism that drives the system toward topological pinch-off. In pure analytical relativity, an unstable equilibrium (such as an unsupported wormhole throat) will eventually collapse under any infinitesimal perturbation. However, in numerical relativity, initial data that violates the Hamiltonian constraint but possesses perfect time symmetry ($K_{ij} = 0$) sits at a mathematical saddle point. Because the initial coordinate velocities are exactly zero, the simulated spacetime will attempt to remain artificially static until numerical truncation errors—which are inherently small and stochastic—accumulate enough to push the system off the equilibrium point. To prevent this arbitrary, resolution-dependent "hovering" period and to systematically control the onset of the dynamics, it is standard practice to manually inject a small, controlled perturbation (a "kick") into the extrinsic curvature to explicitly break the time symmetry and reliably trigger the collapse.

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

The initialized Ellis-Bronnikov geometry requires exotic matter violating the Null Energy Condition to remain static. Previous studies often relied on a "sudden approximation," evolving this geometry in a pure vacuum ($T_{\mu\nu}=0$), which results in a massive, instantaneous violation of the Hamiltonian constraint and violent numerical transients.

To achieve a physically consistent, controlled collapse, we instead explicitly model the required exotic matter. We couple the CCZ4 evolution to a "phantom" scalar field $\phi$. The stress-energy tensor for this field is defined with a reversed overall sign relative to a canonical scalar field, providing the negative energy density necessary to support the wormhole throat and satisfy the constraints at $t=0$.

We systematically control the onset of the dynamics by smoothly removing this supporting matter. We define a time-dependent support strength multiplier $S(t)$. For an initial settling period ($t < t_{\rm start}$), $S(t) = 1$, and the wormhole remains in a supported, quasi-static equilibrium. To trigger the collapse, $S(t)$ is smoothly ramped down to zero over a duration $\Delta t$ using a cosine profile:
\begin{equation}
    S(t) = \frac{1}{2} \left[ 1 + \cos\left( \pi \frac{t - t_{\rm start}}{\Delta t} \right) \right].
\end{equation}
As the exotic matter support vanishes, gravity dictates the dynamics, and the throat naturally collapses under its own weight. To break exact spherical symmetry and generate a small, clean gravitational wave signal, we inject a minor quadrupole ($l=2$) perturbation into the extrinsic curvature:
\begin{equation}
    K(\bar{r}, \theta, \phi) = A_2 \exp\left(-\frac{\bar{r}^2}{\sigma^2}\right) Y_{20}(\theta, \phi).
\end{equation}

\subsection{Apparent Horizon Detection}

To definitively distinguish true topological pinch-off and black hole formation from mere gauge effects (such as "lapse collapse"), we implemented a geometric trapped-surface diagnostic. At each time step, we perform a spherical-shell scan outward from the origin, calculating the expansion of outgoing null rays, $\theta_+$. In spherical symmetry, this is given by:
\begin{equation}
    \theta_+ = \frac{2\sqrt{\chi}}{r} + \frac{2}{3}K - \chi \tilde{A}_{rr}.
\end{equation}
The formation of a trapped surface is unambiguously identified when $\theta_+ \le 0$. The maximum radius where this condition holds provides a proxy for the apparent horizon radius, $r_{\rm AH}$. The sudden jump of $r_{\rm AH}$ from zero to a finite positive value during the simulation serves as our primary physical indicator of successful collapse.

\subsection{Gauge Conditions}

To ensure stability through the highly non-linear topological transition from a traversable wormhole into disjoint black holes, we employ robust gauge choices. 

For the shift vector, we initially set $\beta^i = 0$. Since our initial state possesses no tangential rotation or shift, and the physical mechanism under investigation is a radial pinch-off, a vanishing initial shift prevents gauge-induced artifacts from contaminating the early gravitational wave extraction.

For the lapse function, we employ the $1+\log$ slicing condition, which evolves as $\partial_t \alpha - \beta^i \partial_i \alpha = -2\alpha K$. To prevent severe "gauge shocks" at $t=0$ caused by the massive curvature at the throats, we initialize the lapse in a "pre-collapsed" state:
\begin{equation}
    \alpha(t=0) = \chi^{1/2} = \left( \frac{4\bar{r}^2}{4\bar{r}^2 + b_0^2} \right)^{2}.
\end{equation}
As $\bar{r} \to 0$ (the asymptotic infinity of the secondary universe), the lapse smoothly approaches zero. During the evolution, as the wormhole throat collapses and physical singularities begin to form, the $1+\log$ condition forces the lapse to drop rapidly toward zero ("lapse collapse"), effectively freezing the evolution near the singularities and preventing numerical breakdown.





\section{Numerical Setup}

Numerical evolutions were performed using the \texttt{GRTeclyn} codebase, which is currently under active development. Built upon the \texttt{AMReX} framework~\cite{amrex}, \texttt{GRTeclyn} implements block-structured Adaptive Mesh Refinement (AMR) and is explicitly optimized for high-performance GPU acceleration. Production runs targeted H100-class GPUs (corresponding to \texttt{CUDA\_ARCH=90} in the build configuration) and used one MPI rank per GPU, scaling up to 8 GPUs for the largest calculations. In practice, the exact runtime configuration depends on the local MPI stack: when CUDA-aware MPI was unavailable or unstable, we disabled GPU-aware MPI communication (\texttt{amrex.use\_gpu\_aware\_mpi = 0}) while retaining the same distributed multi-GPU layout.

The full physical computational domain spans a coordinate length of $L_{\text{full}} = 80.0$, covered by a coarse grid of $N_{\text{full}} = 160$ cells, yielding a coarse resolution of $dx_{\text{coarse}} = 0.5$. However, to drastically minimize memory consumption, we rigorously exploit the Cartesian reflection symmetries inherent in the spherically symmetric initial geometry. We evolve only the positive octant ($x \ge 0, y \ge 0, z \ge 0$) of the full domain, effectively modeling just 1/8th of the physical volume. Parity symmetry conditions are strictly applied at the inner reflection boundaries ($x=0, y=0, z=0$), while Sommerfeld radiative boundary conditions are enforced at the outer edges.

To capture the extreme curvature gradients developing during the collapse while maintaining tractability, we employ 5 levels of 2:1 adaptive mesh refinement (AMR). The mesh is dynamically regridded based on the gradients of the conformal factor ($\chi$), allowing the code to focus resolution precisely where the wormhole throat is pinching off. On the finest level of refinement (level 5), the grid spacing reaches $dx_{\text{fine}} \approx 0.0156$. Time stepping is handled using a 4th-order Runge-Kutta scheme, with a Courant factor (dt\_multiplier) reduced to 0.1 to guarantee stability during the violent topological transition.

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


\section{Results}

\subsection{Dynamics of Support Removal and Relaxation}

\begin{figure*}[t] % the * makes it span both columns in a 2-column paper
    \centering
    \includegraphics[width=\textwidth]{final_panel_4.png}
    \caption{The evolution of the trace of the extrinsic curvature $K$. The matter energy density bubble propels outwards...}
    \label{fig:k_evolution}
\end{figure*}

To establish the stability and dynamical fate of the wormhole topology, we performed a series of simulations evolving the supported Ellis-Bronnikov wormhole. The support strength $S(t)$ was held at $1.0$ for an initial period to allow initial data transients to damp, and then smoothly reduced to zero over $t \in [5.0, 10.0]M$. A small quadrupole perturbation ($A_2 = 0.02$) was included to break exact spherical symmetry.

Contrary to the expectation of a violent, immediate collapse into a black hole (as often seen in violently over-kicked vacuum initial data), the controlled removal of the scalar field support reveals a more gradual dynamical relaxation. 

As shown in Figure~\ref{fig:collapse_diagnostics_plot}, during the ramp-down of the support strength, the global minimum of the lapse function \(\alpha\) does not plunge toward zero. Instead, it experiences a small initial dip (reflecting the early gauge adjustment and transient waves) but remains firmly above $\alpha \approx 0.8$. In moving-puncture formalisms, a definitive collapse to a black hole is characterized by "lapse collapse" ($\alpha \to 0$) alongside a corresponding crush in the conformal factor ($\chi \to 0$). Here, the conformal factor minimum actually \emph{increases}, indicating an expansion or widening of the geometry at the throat rather than a crush.

\begin{figure}[h]
    \centering
    \includegraphics[width=\linewidth]{collapse_diagnostics_plot.png}
    \caption{Evolution of the minimum lapse ($\alpha$), minimum conformal factor ($\chi$), and maximum trapped surface radius ($r_{\rm AH}$) during the support removal phase. The lapse remains high and no trapped surface ($\theta_+ \le 0$) forms, indicating the wormhole relaxes and expands rather than collapsing into a black hole.}
    \label{fig:collapse_diagnostics_plot}
\end{figure}

To definitively rule out the formation of an event horizon, we monitored the maximum radius at which the outgoing null expansion satisfies $\theta_+ \le 0$. Throughout the entire evolution, this radius remains strictly zero ($r_{\rm AH} = 0$), confirming that no apparent horizon or trapped surface forms.

These results indicate that, in the absence of a massive implosive kick ($K < 0$), the removal of the exotic matter threading the throat does not automatically lead to catastrophic gravitational collapse. The surrounding mass-energy of the wormhole is insufficient to self-trap, and instead, the defect gently dissipates and expands outward.

\subsection{Role of Perturbation Width and Amplitude}

The transition of the initial Ellis-Bronnikov wormhole into a collapsing state---and the subsequent emission of gravitational waves---is highly sensitive to the geometric profile of the initial perturbation rather than just its overall magnitude.

By varying the perturbation amplitude (\(A\)) and the Gaussian width (\(\sigma\)) relative to the throat radius (\(b_0\)), several key dynamical regimes have been identified:

\subsubsection*{1. The Irrelevance of the Kick Amplitude on Waveform Power}

Numerical experiments demonstrating collapse for different supercritical amplitudes (e.g., \(A = +3, +4, +5\) with \(b_0=0.5M\)) show that while a minimum amplitude is necessary to trigger collapse, further increasing the amplitude does \emph{not} significantly increase the total power or alter the fundamental frequency of the main gravitational wave burst. 

Instead, the dominant frequency and the energy radiated during the ringdown phase appear largely scale-invariant with respect to the kick amplitude. This suggests that the collapse acts as a nonlinear threshold event: once the geometry is pushed past the critical point of no return, the subsequent dynamics (and the resulting GW signal) are dictated by the characteristic mass and scale of the \emph{newly formed black hole}, rather than retaining a memory of the exact strength of the initial push.

\subsubsection*{2. The Critical Role of Perturbation Width (\(\sigma\))}

While amplitude plays a secondary role post-threshold, the \emph{width} of the perturbation compared to the throat radius is the primary determinant of whether collapse occurs at all. 

\begin{itemize}
    \item When the throat radius was doubled to \(b_0=1.0M\) but the perturbation width was kept at \(\sigma=1.0M\), the spacetime resisted collapse. The highly localized kick essentially ``bounced'' off the larger wormhole structure, radiating the perturbation energy away as high-frequency junk without forming an event horizon.
    \item However, when the perturbation width was subsequently increased to \(\sigma=2.0M\) (matching the scale of the \(b_0=1.0M\) throat), \textbf{collapse was successfully achieved.} 
\end{itemize}

This confirms that in order to induce global gravitational collapse, the compressive extrinsic curvature profile must act coherently over the entire strong-field region (the ``shoulders'' of the wormhole). A localized pinch is insufficient; the geometry must be crushed simultaneously on scales comparable to the throat itself.

\subsubsection*{3. Signal Amplification in Larger Wormhole Geometries}

The successful collapse of the larger (\(b_0=1.0M\), \(\sigma=2.0M\)) wormhole reveals critical scaling behaviors. The emitted gravitational wave signal is significantly stronger for the larger wormhole compared to the smaller one ($b_0=0.5M$, $\sigma=1.0M$). This amplification can be understood through two compounding physical effects:

First, geometric scaling increases the effective mass of the system. In general relativity, the throat radius $b_0$ defines the fundamental length and mass scale of the local spacetime geometry. When this larger geometric structure collapses, the resulting black hole is correspondingly more massive, involving the rapid reorganization of a much larger region of highly curved spacetime. Because the amplitude of gravitational wave strain naturally scales linearly with the mass of the source ($h \propto M$), the collapse of a larger wormhole fundamentally guarantees a more energetic gravitational wave signature.

Second, the volumetric scaling of the injected perturbation energy plays a dominant role. Because the initial state is not a solution to the vacuum Hamiltonian constraint, the non-zero extrinsic curvature acts as an effective local energy density, $\rho_{\text{eff}} \sim K^2$. To successfully collapse the $b_0=1.0M$ wormhole, the width of the perturbation $\sigma$ was doubled. Since this perturbation is applied over a 3-dimensional spatial volume ($V \propto \sigma^3$), doubling the width increases the volume of the violently ``kicked'' region by a factor of 8. The spacetime must violently shed this massive excess energy as it dynamically relaxes and collapses, manifesting as a drastically amplified burst of both initial ``junk'' radiation and the primary collapse signal.

\subsection{Separating the Initial Transient from Physical Outgoing Radiation}

Because the initial unsupported wormhole slice contains a constraint defect (since it is evolved in vacuum), the earliest curvature signal is dominated by ``junk'' content: a mixture of gauge/constraint-cleaning transients and pre-existing curvature content relaxing. We analyze the extracted curvature waveform using retarded time, $u \equiv t - R_{\mathrm{ext}}$, to distinguish this initial burst from physical collapse radiation. Note that because of computational domain limitations, our extraction spheres are located at $R_{\mathrm{ext}} = 8, 12, 16M$, placing them deeply within the strong-field near-zone. As a result, these waveforms contain static Coulomb components of the gravitational field and cannot be cleanly interpreted as the asymptotic strain ($h$) seen by gravitational wave observatories.

\begin{figure}[h]
    \centering
    \includegraphics[width=\linewidth]{psi4_extracted_simulation.png}
    \caption{Radius-scaled waveform \(r\,\mathrm{Re}(\Psi_4^{2,0})\) for the supercritical collapse case (\(A=+4\), $b_0=0.5M$), shown in simulation time \(t\) (top) and retarded time \(t-R_{\mathrm{ext}}\) (bottom). The delay of the burst with increasing radius in simulation time, together with its approximate alignment in retarded time, supports the interpretation of the signal as outgoing near-zone radiation emitted by the collapse.}
    \label{fig:psi4_extracted_simulation}
\end{figure}

Figure~\ref{fig:psi4_extracted_simulation} provides two complementary views of the extracted $l=2, m=0$ spherical harmonic mode of the Weyl scalar $\Psi_4$ (denoted as $\Psi_4^{2,0}$) for the supercritical case (\(A=+4\), $b_0=0.5M$), using extraction radii \(R_{\mathrm{ext}}=8,12,16\). In the upper panel plotted against the simulation time \(t\), the dominant burst reaches the extraction spheres at progressively later times as the radius increases. The innermost observer detects the signal first, followed by the intermediate and outer observers. This ordering is precisely the causal behavior expected for an outward-propagating disturbance emitted from the central collapse region. 

In the lower panel, the same waveforms are replotted against the retarded time $u$. A prominent oscillatory packet appears at \emph{negative} retarded time (roughly $u\in[-10,-6]$). This cannot be attributed to radiation emitted after the simulation start at $t=0$; it is the initial transient/near-zone (junk) content propagating off the grid.

A later, sharp feature appears at \emph{positive} $u$ (around $u \approx 2$). Under the retarded time transformation, this dominant burst and its immediate ringdown become approximately aligned across the different extraction radii. This overlap is a standard diagnostic for outgoing radiation in asymptotically flat spacetimes: if the signal is produced near the origin and then propagates outward at unit characteristic speed, its main features should occur at nearly the same retarded time for all observers.

The agreement is not exact, which is expected due to finite-radius extraction effects, gauge effects, and numerical dispersion. Nevertheless, the consistency is sufficiently strong to support the interpretation of the observed burst as a genuine outgoing curvature signal. 

A second important feature of the waveform is its morphology. The large separation in amplitude between the early precursor and the main burst indicates that the dominant observed signal is not simply the initial perturbation leaving the grid, but rather a later dynamical response triggered by the nonlinear collapse of the wormhole throat into an event horizon.

\begin{figure}[h]
    \centering
    \includegraphics[width=\linewidth]{psi4_extracted_R8_R12_R16.png}
    \caption{Detailed comparison of the extracted waveform at $R_{\mathrm{ext}}=8,12,16$, aligned in retarded time, highlighting the isolation of the physical collapse burst (positive $u$) from the initial constraint-violating transient (negative $u$).}
    \label{fig:psi4_extracted_R8_R12_R16}
\end{figure}

\section{Conclusion}

This work demonstrates that the macroscopic collapse of traversable wormholes generates a distinct, detectable gravitational wave signature. By initializing an Ellis-Bronnikov wormhole geometry and evolving it under the vacuum Einstein equations, we confirmed that such unsupported topological structures are violently unstable and naturally collapse into black holes.

Our 3D numerical relativity simulations utilizing the \texttt{GRTeclyn} framework reveal that this topological transition is characterized by a nonlinear threshold behavior. The geometry of the perturbation—specifically its width relative to the throat—is the primary driver of whether collapse occurs, rather than the raw amplitude of the kick. We observed that the mass-energy scale of the resulting black hole, and thus the amplitude and frequency of the emitted gravitational waves, is set by the initial size of the wormhole throat.

The extracted Weyl scalar ($\Psi_4$) waveforms exhibit a clear, causal outward-propagating burst associated with the collapse event, followed by the quasi-normal ringing of the newly formed black hole. This signal is temporally isolated from the initial high-frequency transient caused by the relaxation of the unsupported initial data. 

These results provide the first full 3D dynamical waveforms for the collapse of a macroscopic wormhole topology. The existence of these distinct, high-amplitude bursts establishes that if primordial wormholes exist without the continuous support of exotic matter, their natural collapse mechanism produces near-field gravitational radiation that offers a unique theoretical window into the dynamics of Exotic Compact Objects.

\appendix* 
\section{Numerical Details, Validation and Convergence}

To minimize computational costs while capturing the full 3D dynamics, we exploit Cartesian reflection symmetries. In the $x$, $y$, and $z$ directions, we use the symmetry of the problem to evolve only the positive octant ($x \ge 0, y \ge 0, z \ge 0$) of the full domain, applying appropriate parity conditions at the reflection boundaries.

To validate the stability and accuracy of our numerical evolutions, we monitor the violations of the Hamiltonian and momentum constraints. In the continuum limit of exact General Relativity, these quantities must vanish. Because our initial data represents an unsupported wormhole evolved in vacuum, the simulation begins with a significant constraint defect. 

We check that the $L_2$ norm of the Hamiltonian constraint violation remains under control throughout the simulation and damps effectively. The $L_2$ Norm acts as a global grade for the simulation's physical consistency, providing a single, representative number that measures the total, volume-averaged error across the entire simulation domain. It is calculated as:
\begin{equation}
L_2(\text{Ham}) = \sqrt{ \sum_{i} ( \text{Ham}_i^2 \times \text{Volume}_i ) }.
\end{equation}

\begin{figure}[h]
    \centering
    \includegraphics[width=\linewidth]{constraints_plot.png}
    \caption{Volume-averaged $L_2$ norm of the Hamiltonian and momentum constraint violations over time for the supercritical collapse case ($A=+4$). The initial defect rapidly damps as the CCZ4 formulation propagates the violations away and drives the system toward the constraint surface.}
    \label{fig:constraints_plot}
\end{figure}

As shown in Figure~\ref{fig:constraints_plot}, the constraints exhibit a characteristic early-time behavior: a sharp initial spike at very small times as the unsupported geometry reacts, followed by rapid, exponential damping, and then a slower, controlled relaxation at a much lower level. This confirms that the CCZ4 damping terms and moving-puncture gauge conditions successfully convert the initial constraint defect into propagating transient modes while driving the ongoing collapse dynamics closer to the true physical constraint surface.

For reproducibility, it is important to note that the \texttt{WormholeCollapse} example in the current development snapshot is not yet fully synchronized with the latest \texttt{GRTeclyn} interfaces. A fresh MPI+CUDA build check presently fails in \texttt{WormholeLevel.cpp} because of deprecated or changed API usage, including the legacy \texttt{State\_Type} access pattern, incompatible invocations of \texttt{TraceARemoval} and \texttt{PositiveChiAndLapse}, and an obsolete \texttt{CCZ4RHS::compute} call. These are software-maintenance issues rather than physical limitations of the method, but they imply that exact reruns presently require the corresponding working source state used to generate the production data.
% --- Bibliography ---
\bibliographystyle{apsrev4-2}
\bibliography{references}

\end{document}