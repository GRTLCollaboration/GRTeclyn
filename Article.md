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
Wormholes are topological bridges permitted by General Relativity, but the best-known traversable solutions require exotic matter that may not survive beyond the earliest stages of the universe. If primordial wormhole-like defects were seeded from quantum foam and subsequently stretched to macroscopic scales during inflation, their later evolution would be governed not by a sustained exotic source but by the unsupported relaxation of the geometry itself. In this work, we investigate the nonlinear dynamics of such unsupported wormhole-like topological defects under localized extrinsic-curvature perturbations using full 3D numerical relativity with \texttt{GRTeclyn}. We find a bifurcation between sub-critical relaxation without detected trapping and super-critical self-trapping, in which a strong implosive kick drives topological pinch-off while the extracted near-zone curvature is dominated by transient junk radiation. These results clarify the vacuum dynamical fate of unsupported exotic geometries and provide a concrete framework for studying relic topological defects without assuming persistent exotic matter.
\end{abstract}


% The \maketitle command generates the title block AFTER the abstract
\maketitle

% --- Main Content of the Proposal ---
\section{Introduction}
Wormholes—topological bridges connecting distinct regions of spacetime or separate universes—stand among the most fascinating predictions of General Relativity. The geometric possibility of such structures was first identified by Ludwig Flamm in 1916 \cite{flamm16}, shortly after the discovery of the Schwarzschild solution. The concept was formally introduced into physics by Einstein and Rosen in 1935 as an attempt to model elementary particles as non-singular field configurations \cite{einstein35}. Later work on the Ellis-Bronnikov class of geometries established explicit traversable wormhole solutions supported by non-standard stress-energy, and in 1957 John Archibald Wheeler coined the term "wormhole" while investigating the topological foam of spacetime at the quantum scale \cite{misner57}.

For decades, these solutions were considered topological curiosities, primarily because the standard Einstein-Rosen bridge is non-traversable; the throat pinches off faster than even a photon can cross it. The modern era of wormhole physics accelerated when Morris and Thorne demonstrated that static, traversable wormholes are valid solutions to the Einstein equations, provided the throat is threaded by ``exotic matter'' that violates the null energy condition \cite{morris88}. Consequently, much of the literature has focused on how such stress-energy might be engineered, minimized, replaced in modified gravity \cite{visser95, lobo05}, or sourced by astrophysical dark matter profiles \cite{garattini2026}. Additionally, recent stationary models demonstrate that rapid rotation can arbitrarily reduce the amount of exotic matter required to keep the throat open \textbf{\cite{uemichi2026}}, pointing toward spin as a potential stabilizing mechanism. That perspective is natural if one is attempting to engineer or stabilize a wormhole for transport, but it is not the only physically relevant question.

The motivation for the present work, however, diverges entirely from the search for stability. In Wheeler's picture of spacetime foam, transient microscopic wormholes may be produced at very early times, and inflation may in principle stretch some topological defects to much larger scales \cite{kardashev07, garcia16}. If such relic objects ever existed, the exotic component that originally sustained them need not remain present at late times; it may have been tied to an inflationary sector, a transient field configuration, or some other early-universe mechanism that subsequently decayed away. In that case the relevant physical problem is not how to maintain a supported traversable wormhole indefinitely, but how an unsupported wormhole-like geometry relaxes, disperses, or collapses once left to evolve under ordinary vacuum gravity. Indeed, recent theoretical work by \citet{dimaschko25} suggests that any bound state of a traversable wormhole is non-stationary. Rather than trying to keep the bridge open, this work investigates the subsequent dynamics of unsupported topological defects under external perturbation.

The dynamical fate of such traversable wormholes was first established in a pioneering numerical study by Shinkai and Hayward \cite{shinkai02}. Evolving the classic Morris-Thorne (Ellis-Bronnikov) metric in 1D spherical symmetry, they demonstrated that these structures are highly unstable against perturbations. They observed that the wormhole throat suffers a bifurcation of horizons, leading either to an inflationary expansion (given an injection of negative energy) or a rapid collapse into a black hole (given positive energy), as conceptually illustrated in Fig.~\ref{fig:wormhole_collapse_stages}. Crucially, the introduction of even a minute pulse of normal matter—representing, for instance, a traveler attempting to traverse the bridge—inevitably forces the wormhole to collapse into a black hole, sealing off causal connection between the two regions. Their work elegantly proved that exotic topological structures can be rigorously studied using the same local, dynamical methods applied to standard black holes, effectively unifying the two phenomena. 

While Shinkai and Hayward established the ultimate instability and evolutionary paths (contraction vs. expansion) of the wormhole topology, their strict assumption of 1D spherical symmetry inherently precluded the study of non-spherical dynamics and any extracted curvature signal. By Birkhoff's theorem, a purely spherically symmetric vacuum evolution cannot radiate. To determine how unsupported wormhole-like defects behave once symmetry is broken, one must therefore move beyond 1D restrictions and perform full 3D, non-linear dynamical simulations.

Simulating an unsupported exotic geometry is computationally challenging, but it is also the physically relevant limit if the sustaining exotic source is absent. A useful methodological precedent is the recent 3D study of unsupported Alcubierre-bubble collapse by \citet{clough24}, which showed that numerical relativity can be used to follow the vacuum relaxation of an exotic spacetime after its source is removed. Following a similar philosophy, we study unsupported wormhole dynamics with the \texttt{GRTeclyn} framework—built upon the \texttt{GRChombo} infrastructure \cite{Andrade2021}—which robustly evolves the geometry on an adaptive 3D Cartesian grid. Our aim is not to assume away the missing exotic matter, but to understand the nonlinear fate of the residual topological geometry once that support is gone.

\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{wormhole_collapse_stages.pdf}
\caption{Conceptual 3D embedding diagram illustrating the two dynamical evolutionary paths of an unstable traversable wormhole upon perturbation. \textbf{Top Row (a--c):} The super-critical branch ($K<0$), representing an injection of positive-energy matter. The throat constricts (b) until its topology dynamically pinches off, resulting in the bifurcation of horizons and the formation of a black hole (c) that seals the two universes. \textbf{Bottom Row (d--f):} The sub-critical branch ($K \ge 0$). Evolved as a bare geometry in pure vacuum without a continuous ghost field, the resolved exterior geometry expands and weakens (e). In the unsupported experiment studied here, we only claim relaxation toward a flatter weak-curvature state on the resolved grid, not a proof that the full two-ended topology has become exactly Minkowski (f).}
\label{fig:wormhole_collapse_stages}
\end{figure}

The ultimate goal of this research is to investigate the nonlinear dynamics of unsupported wormhole-like defects and to determine which perturbations lead to dispersal, which lead to trapping, and what sort of near-zone curvature signal accompanies each regime. The advent of GW astronomy by the LIGO-Virgo-KAGRA collaborations \citep{abbott16,abbott18,abbott21,abbott23} has moved the search for Exotic Compact Objects (ECOs) \cite{cardoso19} from pure theory to observational reality, but the first theoretical task is to understand the dynamics of candidate geometries themselves. Because our computational domain restricts extraction to the strong-field regime, we do not aim to provide asymptotic waveform templates for GW observatories. Instead, we use extracted \(\Psi_4\) primarily as a diagnostic of the near-zone vacuum response of the geometry, distinguishing cleaner outgoing signals from junk-dominated transients while mapping the underlying bifurcation between dispersal and pinch-off.

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
For the experiments studied here, the dynamical evolution is governed by the pure vacuum Einstein equations, implemented within the \texttt{GRTeclyn} framework using the conformal and covariant Z4 (CCZ4) formulation. Our objective is to start from the geometry of an Ellis-Bronnikov wormhole at $t=0$, subject to an initial extrinsic-curvature perturbation. The evolution is then driven by the nonlinear adjustment of this unsupported geometry via the CCZ4 constraints.

\subsection{Initial Data: The Morris-Thorne Wormhole in Isotropic Coordinates}

We initialize the simulation using the static, spherically symmetric Morris-Thorne metric:
\begin{equation}
ds^2 = -e^{2\Phi(r)} dt^2 + \frac{dr^2}{1 - b(r)/r} + r^2 (d\theta^2 + \sin^2\theta d\phi^2).
\end{equation}
To balance physical relevance with numerical stability, we adopt the specific functional forms corresponding to the Ellis-Bronnikov wormhole, a canonical massless scalar field solution. We set the redshift function to zero ($\Phi(r) = 0$), which ensures the absence of initial horizons, sets the proper time of Eulerian observers equal to coordinate time at $t=0$, and eliminates radial tidal forces at the throat. We choose the inverse-$r$ shape function $b(r) = b_0^2/r$, where $b_0$ represents the throat radius. This choice guarantees asymptotic flatness as $r \to \infty$.

In standard Schwarzschild-like coordinates, the metric component $g_{rr} = (1 - b_0^2/r^2)^{-1}$ diverges at the throat ($r=b_0$). To resolve this coordinate singularity and simultaneously map both "universes" connected by the bridge onto a single computational grid, we transform to \textbf{isotropic coordinates}. We introduce the isotropic radius $\bar{r}$, related to the proper distance $\ell = \pm \sqrt{r^2 - b_0^2}$ via the embedding $\ell = b_0 \sinh(\ln(2\bar{r}/b_0))$. 

To initialize the wormhole on a 3D Cartesian grid, we transform the standard metric into isotropic coordinates. As demonstrated by Nandi et al.~\cite{nandi2008energetics}, the spatial geometry of the Ellis-Bronnikov wormhole admits a conformally flat representation where the BSSN/CCZ4 conformal factor is strictly given by:
\begin{equation}
    \psi = \sqrt{1 + \frac{b_0^2}{4\bar{r}^2}}.
\end{equation}
This differs fundamentally from the Brill-Lindquist black hole puncture, ensuring the initial slice genuinely captures the topological bridge. In this chart, the primary universe extends to $\bar{r} \to \infty$ ($\psi \to 1$), the throat is located at $\bar{r} = b_0/2$, and the secondary universe is compactified such that its asymptotic infinity maps to the grid origin $\bar{r} \to 0$ ($\psi \to \infty$). 

\subsection{Transformation to CCZ4 Variables and Regularization}

For stable numerical evolution in the CCZ4 formulation, we utilize the conformal transverse-traceless decomposition. The physical spatial metric $\gamma_{ij}$ is represented by a conformal metric $\tilde{\gamma}_{ij}$ with unit determinant, and a regularized conformal factor $\chi$:
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

The classic Ellis-Bronnikov traversable wormhole is an exact solution to the Einstein Equations sourced by a massless ghost scalar field, which provides the negative energy density required to violate the Null Energy Condition and support the throat. In this work, we investigate the dynamical fate of this topology in the limit where this exotic ghost field is suddenly removed. To achieve this, we employ a "sudden approximation" approach: rather than defining and dynamically depleting an active matter field, we initialize the exact spatial geometry of the Ellis-Bronnikov wormhole at $t=0$ but deliberately set the stress-energy tensor to zero ($T_{\mu\nu}=0$).

By starting with this "empty" spacetime, we are directly simulating the bare residual topology the exact moment after the supporting matter has been removed. Consequently, the initial slice purposefully carries a Hamiltonian constraint defect that corresponds exactly to the missing exotic matter. This allows the CCZ4 formulation to naturally evolve the vacuum relaxation of the unsupported geometry.

Because the wormhole throat possesses strong negative spatial curvature ($R^{(3)}<0$), removing the matter that previously balanced it causes the bare geometry to naturally undergo a volumetric expansion. To rigorously map the full phase space of this topological defect—specifically the threshold between this natural dispersion and complete gravitational collapse—we must introduce a competing implosive effect. We therefore inject a localized, tunable compressive perturbation in the trace of the extrinsic curvature ($K$) at $t=0$. This momentum kick models a violent external perturbation (such as an incoming gravitational wave or matter pulse) impacting the unsupported throat, allowing us to successfully push the geometry across the critical threshold into dynamical pinch-off and black hole formation.

Specifically, the localized perturbation in the trace of the extrinsic curvature $K$ at $t=0$ is given by:
\begin{equation}
    K(\bar{r}, \theta, \phi) = \left[A_0 + A_2 Y_{20}(\theta, \phi)\right]
    \exp\left(-\frac{\bar{r}^2}{\sigma^2}\right).
\label{eq:kick_profile}
\end{equation}
By tuning the amplitude $A_0$, this setup serves as a tunable critical-collapse experiment, explicitly mapping the bifurcation between the natural sub-critical dispersion (where the negative curvature dominates) and the super-critical dynamical pinch-off (where the implosive momentum dominates, forcing the formation of an apparent horizon).

Furthermore, the quadrupolar amplitude $A_2$ explicitly breaks the exact spherical symmetry of the initial spacetime. Because Birkhoff's theorem dictates that a perfectly spherically symmetric vacuum evolution cannot radiate, this \(Y_{20}\) deformation—modeling an asymmetric tidal encounter—allows us to inspect the extracted Weyl curvature \(\Psi_4\) and distinguish cleaner outgoing structure from purely gauge- or constraint-driven transients.

The dynamics are therefore driven entirely by the interplay between the unsupported vacuum geometry and this initial momentum fluctuation, cleanly dictating whether the topological defect disperses or undergoes genuine gravitational collapse.

\subsection{Apparent Horizon Detection}

During the period in which this work was carried out, \texttt{GRTeclyn} did not yet include a production-ready full apparent-horizon finder for these runs. To distinguish genuine topological pinch-off and black-hole-like self-trapping from gauge effects such as lapse collapse, we therefore implemented a geometric trapped-surface proxy on coordinate spheres centered on the throat. At each time step, we scan outward from the origin and estimate the expansion of outgoing null rays, \(\theta_+\). For a conformally flat spatial metric \(\gamma_{ij}=\chi^{-1}\delta_{ij}\), the spherical expression used in the code is
\begin{equation}
    \theta_+ = \frac{2\sqrt{\chi}}{r} - \frac{\partial_r \chi}{\sqrt{\chi}} + \tilde{A}_{rr} - \frac{2}{3}K.
\end{equation}
The formation of a trapped surface is identified when \(\theta_+ \le 0\). The maximum radius where this condition holds provides a proxy for the apparent-horizon radius, \(r_{\rm AH}\). Although this quantity does not replace a full elliptic horizon finder, it gives a practical collapse diagnostic for the parameter study reported here.

\subsection{Gauge Conditions}

To ensure stable evolution across both dispersive and collapsing branches, we employ robust gauge choices. 

For the shift vector, we initially set $\beta^i = 0$. Since our initial state possesses no tangential rotation or shift, and the physical mechanism under investigation is a radial pinch-off, a vanishing initial shift prevents gauge-induced artifacts from contaminating the early gravitational wave extraction.

For the lapse function, we employ the $1+\log$ slicing condition, which evolves as $\partial_t \alpha - \beta^i \partial_i \alpha = -2\alpha K$. We initialize the lapse as
\begin{equation}
    \alpha(t=0) = 1,
\end{equation}
so that the subsequent lapse evolution is generated dynamically by the geometry itself rather than being biased toward an already ``pre-collapsed'' gauge profile. In collapsing runs, the \(1+\log\) condition drives the lapse toward zero (``lapse collapse''), effectively freezing the evolution near the singular region and preventing numerical breakdown; in sub-critical runs, by contrast, the lapse need not collapse and can instead relax back toward unity as the geometry disperses.





\section{Numerical Setup}

Numerical evolutions were performed using the \texttt{GRTeclyn} codebase, which is currently under active development. Built upon the \texttt{AMReX} framework~\cite{amrex}, \texttt{GRTeclyn} implements block-structured Adaptive Mesh Refinement (AMR) and is explicitly optimized for high-performance GPU acceleration. Production runs targeted H100-class GPUs (corresponding to \texttt{CUDA\_ARCH=90} in the build configuration) and used one MPI rank per GPU, scaling up to 8 GPUs for the largest calculations. In practice, the simulation workflow follows the same broad philosophy as other recent vacuum relaxations of exotic geometries, including the unsupported Alcubierre-bubble study of \citet{clough24}: evolve the bare geometry directly and diagnose whether it disperses, self-traps, or converts its initial defect into propagating transients.

The full physical computational domain spans a coordinate length of \(L_{\text{full}} = 80.0\), covered by a coarse grid of \(N_{\text{full}} = 160\) cells, yielding a coarse resolution of \(dx_{\text{coarse}} = 0.5\). To drastically reduce the computational cost, we exploit the Cartesian reflection symmetries of the single-throat initial data and evolve only the positive octant (\(x \ge 0, y \ge 0, z \ge 0\)), effectively modeling just one eighth of the physical volume. Parity symmetry conditions are imposed at the inner reflection boundaries, while Sommerfeld radiative boundary conditions are enforced at the outer edges.

To capture the extreme curvature gradients developing during the collapse while maintaining tractability, we employ 5 levels of 2:1 adaptive mesh refinement (AMR). The mesh is dynamically regridded based on the gradients of the conformal factor ($\chi$), allowing the code to focus resolution precisely where the wormhole throat is pinching off. On the finest level of refinement (level 5), the grid spacing reaches $dx_{\text{fine}} \approx 0.0156$. Time stepping is handled using a 4th-order Runge-Kutta scheme, with a conservative Courant factor (dt\_multiplier) reduced to 0.05 to maintain stability during the constraint relaxation and the abrupt collapse branches.

The primary objective of this work is to investigate the collapse dynamics subject to Gaussian extrinsic curvature perturbations (Eq.~\ref{eq:kick_profile}). In this manuscript we focus on two representative evolutions that illustrate the sub-critical and super-critical branches: a sub-critical dispersive run with throat radius \(b_0=0.5\) and kick width \(\sigma_K=1.0\), and a super-critical collapsing run with \(b_0=0.2\) and \(\sigma_K=0.2\). For reproducibility, we also report the CCZ4 constraint-damping strength \(\kappa_1\) and the Kreiss--Oliger dissipation coefficient \(\sigma_{\rm KO}\), which we tune to control early-time junk radiation in the extracted diagnostics.


\section{Results}

\subsection{Vacuum evolution of an unsupported topology}
In classical General Relativity, a traversable wormhole inherently requires exotic matter violating the Null Energy Condition to remain static. Recent no-go theorems confirm that this requirement is a deeply fundamental consequence of the geometric flaring-out condition, persisting even in modified frameworks such as Unimodular Gravity \cite{cataldo2026}. In particular, in the Ellis--Bronnikov geometry the spatial Ricci scalar on the initial hypersurface is negative in the strong-field throat region, \(R^{(3)} < 0\). Consequently, when this geometry is initialized and evolved in pure vacuum (\(T_{\mu\nu}=0\)) without exotic support, the initial slice contains a substantial Hamiltonian constraint defect, and any nontrivial localized \(K\) kick introduces additional momentum-constraint violation. We therefore interpret the unsupported evolution as an out-of-equilibrium vacuum-relaxation problem rather than as a stationary vacuum solution.

Within the CCZ4 formulation, the constraint-damping mechanism rapidly converts the initial defect into propagating transient modes. After this early-time relaxation, the subsequent dynamics exhibit a clear nonlinear threshold separating two regimes: \emph{topological dissipation} (sub-critical evolutions) and \emph{dynamical pinch-off} (super-critical collapse).

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
\caption{Extracted radius-scaled curvature signal \(r\,\mathrm{Re}(\Psi_4^{2,0})\) for the sub-critical evolution. In this regime the retarded-time alignment across extraction radii is consistent with a modest outgoing near-zone curvature transient emitted during the vacuum relaxation of the unsupported throat.}
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

\noindent\textbf{Run parameters.} For the successful collapse shown in Fig.~\ref{fig:supercritical_diagnostics}, we used \(b_0=0.2\) with a strong monopolar compressive kick \(A_0=5.0\), a small quadrupolar deformation \(A_2=0.05\), and width \(\sigma_K=0.2\). In contrast to the sub-critical case, we did not use large damping/dissipation; we set \(\kappa_1=0.1\) and \(\sigma_{\rm KO}=0.5\) to avoid over-damping the nonlinear collapse while still controlling numerical noise.

\begin{figure*}[t]
    \centering
    \includegraphics[width=0.92\textwidth]{collapse_diagnostics_supercritical.pdf}
    \caption{Diagnostics of a super-critical collapse. The minimum lapse (\(\alpha\)) and conformal factor (\(\chi\)) plunge toward zero, demonstrating the time-freezing and infinite spatial stretching characteristic of the moving-puncture gauge. Concurrently, the minimum null expansion proxy \(\min(\theta_+)\) crosses zero and becomes negative, causing the trapped-surface radius proxy \(r_{\rm AH}\) to jump from zero to a finite value. The subsequent step-wise growth in \(r_{\rm AH}\) is an expected artifact of the discrete Cartesian AMR grid tracking the outward coordinate stretching of the trapped region. Together, these features are consistent with dynamical self-trapping and pinch-off in the unsupported experiment.}
    \label{fig:supercritical_diagnostics}
\end{figure*}
To model a strong external perturbation (e.g. an early-universe fluctuation) we seed a localized inward perturbation in the trace of the extrinsic curvature. Because a perfectly spherically symmetric vacuum collapse cannot radiate (Birkhoff's theorem), we utilize the quadrupolar deformation defined in Eq.~\ref{eq:kick_profile} to break exact symmetry.
When the compressive amplitude exceeds a critical threshold, the implosion overcomes the topological repulsion of the throat. In this super-critical branch, \(\min(\theta_+)\) crosses zero, \(r_{\rm AH}\) becomes positive, and the evolution becomes consistent with black-hole-like self-trapping and throat pinch-off.

\subsection{Gravitational-wave signatures (near-zone)}
The extracted Weyl signal is strongly regime-dependent. In the sub-critical case, the retarded-time waveform in Fig.~\ref{fig:psi_sub_critical} shows the clearest indication of a modest outgoing curvature pulse associated with the relaxation of the unsupported geometry. In the super-critical case, by contrast, the extracted \(\Psi_4\) is dominated by early-time junk radiation and constraint-cleaning transients, making it much harder to interpret as a clean physical outgoing burst. We therefore use the super-critical waveform primarily as a diagnostic of contamination in the near zone rather than as evidence for a robust observable emission channel.
\begin{figure}[h]
\centering
\includegraphics[width=\linewidth]{psi_super_critical.pdf}
\caption{Radius-scaled curvature waveform \(r\,\mathrm{Re}(\Psi_4^{2,0})\) for the super-critical run, shown in simulation time \(t\) (top) and retarded time \(u=t-R_{\rm ext}\) (bottom). In this regime the extracted signal is dominated by junk radiation and early transient contamination, so the plot is best interpreted as a near-zone diagnostic of the violent vacuum relaxation rather than as a clean outgoing gravitational-wave template.}
\label{fig:psi_super_critical}
\end{figure}

\section{Discussion: 3D Dynamics versus 1D Instability Proofs}

Our 3D numerical results are consistent with the qualitative bifurcation identified by Shinkai and Hayward \cite{shinkai02}, but the emphasis here is different. Rather than revisiting the 1D instability proof itself, we use full 3D numerical relativity to study how an unsupported wormhole-like geometry relaxes once its sustaining source is absent. Their seminal work showed that the Morris-Thorne (Ellis-Bronnikov) wormhole lies near a bifurcation between contraction and expansion; our calculations revisit that vacuum relaxation problem in a setting where non-spherical dynamics and extracted curvature can also be monitored:
\begin{itemize}
    \item \textbf{Positive Energy (Collapse to a Black Hole):} If perturbed by normal, positive-energy matter (such as a traveling observer, a normal scalar field, or a gravitational wave), gravity dominates and the throat instantly collapses into a black hole, sealing off the two universes.
    \item \textbf{Negative Energy (Inflationary Expansion):} If more exotic negative-energy matter (a ghost field) is injected into the throat, repulsive anti-gravity dominates, causing the throat to explode and driving the universe into an exponential, inflationary expansion.
\end{itemize}

Our simulation results provide a modern 3D picture of this bifurcation, with essential updates resulting from evolving the bare geometry in pure vacuum rather than with an active scalar field:
\begin{itemize}
    \item \textbf{Super-critical trapping:} In the super-critical branch, a sufficiently large compressive kick drives the unsupported geometry through a trapping transition in the spherical proxy. The throat pinches off, \(\theta_+\) becomes negative, and the trapped-surface proxy turns on. This is the regime in which the vacuum defect is strong enough to self-trap rather than disperse.
    \item \textbf{Sub-critical relaxation:} In the sub-critical branch, the unsupported geometry does not self-trap. Instead, the negative-curvature region relaxes, the lapse and conformal factor recover, and the resolved exterior spacetime tends toward a flatter weak-curvature configuration. Without a continuously maintained exotic source, the unsupported throat does not remain stationary.
\end{itemize}

Crucially, our 3D approach reveals behavior that 1D codes inherently cannot capture. Because Shinkai and Hayward utilized a 1D spherically symmetric simulation, Birkhoff's theorem prevented any nontrivial curvature signal from being extracted. By moving the problem to full 3D and explicitly breaking exact symmetry with a quadrupolar perturbation, we can inspect \(\Psi_4\) directly. The key lesson of the present simulations is not that every pinch-off yields a clean gravitational-wave burst, but that the extracted curvature signal depends strongly on the dynamical branch: the sub-critical regime shows the clearest outgoing structure, while the super-critical regime is dominated by junk radiation in the near zone.

\section{Conclusion}
We have investigated the nonlinear vacuum dynamics of unsupported wormhole-like topological defects using full 3D numerical relativity. The central question in this work is not how to sustain a traversable wormhole with exotic matter, but what becomes of the residual geometry once that support is absent and the object is perturbed.

The evolution exhibits a clear threshold behavior. Sub-critical perturbations drive relaxation toward a flatter geometry without detected trapping, while sufficiently strong kicks produce a super-critical branch in which the throat pinches off and the trapped-surface proxy becomes positive. In this sense, unsupported wormholes in vacuum CCZ4 behave as dynamical topological defects with a bifurcation between relaxation and self-trapping rather than as stationary exotic objects. In particular, the unsupported branch studied here does not remain as a static traversable throat on the resolved grid. While quantum backreaction may theoretically stabilize microscopic topological defects at the Planck scale \cite{mehulic2026}, our results demonstrate that once stretched to macroscopic scales where classical vacuum dynamics dominate, they lack the sustained exotic support required to survive. This provides a natural cosmological mechanism explaining why inflationary wormholes do not abundantly persist in the late universe.

The extracted curvature signal is likewise branch-dependent. The clearest outgoing near-zone structure appears in the sub-critical regime, whereas the super-critical waveform is dominated by junk radiation and early transient contamination. These results show that understanding unsupported wormhole dynamics requires tracking the geometry and its diagnostics together: the fate of the topological defect is robustly encoded in the collapse indicators, while the interpretation of \(\Psi_4\) depends sensitively on the regime and on the contamination introduced by the initial constraint defect.

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
\clearpage
\bibliographystyle{apsrev4-2}
\bibliography{references}
\end{document}