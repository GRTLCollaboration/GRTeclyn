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
Wormholes are topological bridges permitted by General Relativity, but the best-known traversable solutions require exotic matter that may not survive beyond the earliest stages of the universe. If primordial wormhole-like defects were seeded from quantum foam and subsequently stretched to macroscopic scales during inflation, their later evolution would be governed not by a sustained exotic source but by the unsupported relaxation of the geometry itself. In this work, we investigate the nonlinear dynamics of such unsupported wormhole-like topological defects under localized initial extrinsic-curvature deformations using full 3D numerical relativity with \texttt{GRTeclyn}. We find that in all cases the unsupported geometry ultimately collapses into a black hole: (i)~in the unperturbed case ($A_0{=}A_2{=}0$), the bare topology spontaneously collapses due to the inherent linear instability of the unsupported geometry, with numerical truncation error providing a sufficient seed perturbation; (ii)~under a compressive initial deformation ($A_0>0$), the collapse is accelerated and accompanied by a near-zone curvature transient whose propagation speed ($v\approx c$) and oscillatory structure are consistent with physical gravitational radiation---in sharp contrast to the superluminal constraint-damping modes that dominate the unperturbed signal. These results demonstrate that unsupported wormhole-like topological defects are generically unstable and provide a concrete framework for studying relic topological defects without assuming persistent exotic matter.
\end{abstract}


% The \maketitle command generates the title block AFTER the abstract
\maketitle

% --- Main Content of the Proposal ---
\section{Introduction}
Wormholes—topological bridges connecting distinct regions of spacetime or separate universes—stand among the most fascinating predictions of General Relativity. The geometric possibility of such structures was first identified by Ludwig Flamm in 1916 \cite{flamm16}, shortly after the discovery of the Schwarzschild solution. The concept was formally introduced into physics by Einstein and Rosen in 1935 as an attempt to model elementary particles as non-singular field configurations \cite{einstein35}, and in 1957 John Archibald Wheeler coined the term ``wormhole'' while investigating the topological foam of spacetime at the quantum scale \cite{misner57}. The first explicit traversable wormhole solutions supported by a phantom scalar field were found independently by Ellis~\cite{ellis73} and Bronnikov~\cite{bronnikov73} in 1973; this Ellis--Bronnikov geometry remains the canonical example studied in the present work.

For decades, these solutions were considered topological curiosities, primarily because the standard Einstein-Rosen bridge is non-traversable; the throat pinches off faster than even a photon can cross it. The modern era of wormhole physics accelerated when Morris and Thorne demonstrated that static, traversable wormholes are valid solutions to the Einstein equations, provided the throat is threaded by ``exotic matter'' that violates the null energy condition \cite{morris88}. Consequently, much of the literature has focused on how such stress-energy might be engineered, minimized, replaced in modified gravity \cite{visser95, lobo05}, or sourced by astrophysical dark matter profiles \cite{garattini2026}. Additionally, recent stationary models demonstrate that rapid rotation can arbitrarily reduce the amount of exotic matter required to keep the throat open \textbf{\cite{uemichi2026}}, pointing toward spin as a potential stabilizing mechanism. That perspective is natural if one is attempting to engineer or stabilize a wormhole for transport, but it is not the only physically relevant question.

The motivation for the present work, however, diverges entirely from the search for stability. In Wheeler's picture of spacetime foam, transient microscopic wormholes may be produced at very early times, and inflation may in principle stretch some topological defects to much larger scales \cite{kardashev07, garcia16}. If such relic objects ever existed, the exotic component that originally sustained them need not remain present at late times; it may have been tied to an inflationary sector, a transient field configuration, or some other early-universe mechanism that subsequently decayed away. The residual geometry would then inherit a primordial volumetric deformation from its cosmological environment---the cumulative effect of inflationary stretching, reheating dynamics, and the macroscopic amplification of quantum vacuum fluctuations. In that case the relevant physical problem is not how to maintain a supported traversable wormhole indefinitely, but how an unsupported wormhole-like geometry, initialized with this primordial geometric deformation, relaxes, disperses, or collapses once left to evolve under ordinary vacuum gravity. Indeed, recent theoretical work by \citet{dimaschko25} suggests that any bound state of a traversable wormhole is non-stationary. Rather than trying to keep the bridge open, this work investigates the subsequent dynamics of unsupported topological defects under such initial geometric deformations.

The dynamical fate of such traversable wormholes was first established in a pioneering numerical study by Shinkai and Hayward \cite{shinkai02}. Evolving the classic Ellis--Bronnikov metric in 1D spherical symmetry, they demonstrated that these structures are highly unstable against perturbations. They observed that the wormhole throat suffers a bifurcation of horizons, leading either to an inflationary expansion (when the initial deformation is rarefactive, $K<0$) or a rapid collapse into a black hole (when the initial deformation is compressive, $K>0$), as conceptually illustrated in Fig.~\ref{fig:wormhole_collapse_stages}. Crucially, even a minute compressive deformation---representing, for instance, a traveler attempting to traverse the bridge---inevitably forces the wormhole to collapse into a black hole, sealing off causal connection between the two regions. Their work elegantly proved that exotic topological structures can be rigorously studied using the same local, dynamical methods applied to standard black holes, effectively unifying the two phenomena. 

While Shinkai and Hayward established the ultimate instability and evolutionary paths (contraction vs. expansion) of the wormhole topology, their strict assumption of 1D spherical symmetry inherently precluded the study of non-spherical dynamics and any extracted curvature signal. By Birkhoff's theorem, a purely spherically symmetric vacuum evolution cannot radiate. To determine how unsupported wormhole-like defects behave once symmetry is broken, one must therefore move beyond 1D restrictions and perform full 3D, non-linear dynamical simulations. The ultimate goal of this research is to investigate the nonlinear dynamics of unsupported wormhole-like defects---both with and without an initial compressive deformation---and to determine whether these geometries collapse, and what near-zone curvature signal accompanies the process.

Recent 3D studies have begun exploring the nonlinear fate of exotic geometries in numerical relativity \cite{clough24}. In this work, we adopt a ``sudden approximation'' in which the exotic source is removed entirely at $t{=}0$, leaving a pure vacuum relaxation problem with an initial Hamiltonian constraint defect. Our aim is not to assume away the missing exotic matter, but to understand the nonlinear fate of the residual topological geometry once that support is gone.

Moreover, any primordial wormhole that survived inflationary stretching would not remain perfectly spherically symmetric: tidal interactions with surrounding density perturbations, the stochastic nature of reheating, and the expansion of the universe all introduce both coherent and aspherical volume changes. The extrinsic-curvature deformation we impose at $t{=}0$ therefore models this physically expected post-inflationary state---a snapshot of the deformed geometry at the moment its exotic support decays---rather than being an externally applied dynamical "kick".

Simulating such an unsupported exotic geometry in full 3D is computationally demanding: the extreme curvature gradients at the wormhole throat require deep adaptive mesh refinement, and the long evolution times needed to distinguish dispersal from collapse push memory and throughput requirements well beyond what CPU-only codes can deliver in a practical timeframe. We therefore study these dynamics with the \texttt{GRTeclyn} framework—a GPU-native numerical-relativity code built upon the \texttt{AMReX} infrastructure \cite{amrex} and descended from \texttt{GRChombo} \cite{Andrade2021}—which implements block-structured adaptive mesh refinement explicitly optimized for high-performance GPU acceleration on H100-class hardware. This capability is essential: the parameter surveys reported here require dozens of full 3D evolutions, each utilizing up to 8 GPUs, a workload that would be prohibitive without dedicated GPU acceleration.

\begin{figure}[htp]
\centering
\includegraphics[width=\linewidth]{wormhole_collapse_stages.pdf}
\caption{Conceptual 3D embedding diagram illustrating the dynamical collapse of an unstable traversable wormhole. \textbf{Top Row (a--c):} A compressive deformation ($K>0$) drives rapid throat constriction (b) and topological pinch-off into a black hole (c). \textbf{Bottom Row (d--f):} Without any applied deformation, the bare topology follows the same qualitative pathway on a slower timescale: the inherent linear instability of the unsupported Ellis--Bronnikov geometry triggers spontaneous collapse (e), ultimately producing a black-hole remnant (f). In both cases, the unsupported wormhole throat cannot remain static once the exotic matter is absent.}
\label{fig:wormhole_collapse_stages}
\end{figure}

The advent of GW astronomy by the LIGO-Virgo-KAGRA collaborations \citep{abbott16,abbott18,abbott21,abbott23} has moved the search for Exotic Compact Objects (ECOs) \cite{cardoso19} from pure theory to observational reality, but the first theoretical task is to understand the dynamics of candidate geometries themselves. We extract the Weyl scalar \(\Psi_4\) at multiple coordinate radii in the near zone; propagation speed analysis across these radii allows us to distinguish near-zone curvature transients consistent with gravitational radiation (propagating at $v=c$) from CCZ4 constraint-damping modes (propagating superluminally). Because the initial extrinsic-curvature deformation does not satisfy the momentum constraint, the extracted signals inevitably contain a contribution from constraint-cleaning modes; a definitive separation of physical radiation from this contamination awaits constraint-satisfying initial data.

Our primary discriminant is the contrast between the unperturbed and perturbed runs, which share identical numerical settings (gauge, resolution, constraint damping, conformal-factor floor). In the unperturbed case---which theoretically produces zero radiation---the extracted \(\Psi_4\) is weak and propagates superluminally, establishing the constraint-noise baseline. In the perturbed case, the signal is two orders of magnitude stronger and propagates at \(v\approx c\). Because both runs start from the same type of constraint-violating initial data and form the same type of remnant, any feature that appears \emph{only} in the perturbed run cannot be attributed to the constraint defect alone. This built-in null test provides the strongest evidence that the perturbed signal contains a physical component.

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
For the experiments studied here, the dynamical evolution is governed by the pure vacuum Einstein equations, implemented within the \texttt{GRTeclyn} framework using the conformal and covariant Z4 (CCZ4) formulation. Our objective is to start from the geometry of an Ellis--Bronnikov wormhole at $t=0$, subject to an initial extrinsic-curvature deformation. The evolution is then driven by the nonlinear adjustment of this unsupported geometry via the CCZ4 constraints.

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

\subsection{Dynamical Triggering via Primordial Geometric Deformation}

The classic Ellis--Bronnikov traversable wormhole is an exact solution to the Einstein equations sourced by a massless phantom scalar field that violates the Null Energy Condition and thereby supports the throat. In this work, we investigate the dynamical fate of this topology in the limit where the phantom source is suddenly removed. Unlike \citet{clough24}, who modelled the collapse of an Alcubierre warp bubble by explicitly evolving the exotic matter as a dynamical stiff fluid with an ad-hoc equation of state, we adopt a ``sudden approximation'': we initialize the exact spatial geometry of the Ellis--Bronnikov wormhole at $t=0$ but deliberately set the stress-energy tensor to zero ($T_{\mu\nu}=0$), rather than defining and dynamically depleting an active matter field.

By starting with this "empty" spacetime, we are directly simulating the bare residual topology the exact moment after the supporting matter has been removed. Consequently, the initial slice purposefully carries a Hamiltonian constraint defect that corresponds exactly to the missing exotic matter. This allows the CCZ4 formulation to naturally evolve the vacuum relaxation of the unsupported geometry.

We emphasize that the sudden approximation is an idealized limiting case---the instantaneous-removal limit of what in reality would be a gradual matter decay---and is not intended as a claim of physical completeness. The full Ellis--Bronnikov spacetime has vanishing ADM mass ($M_{\rm ADM}=0$) when both asymptotic ends are included; on our computational domain, which resolves only one asymptotic end on a finite grid, the constraint defect carries effective gravitational energy that the CCZ4 system redistributes dynamically. The remnant mass that emerges ($M_{\rm BH}\sim 0.1$ in code units) therefore originates from this constraint-energy redistribution rather than from a conserved ADM integral. We regard this as a consistent---though not fully rigorous---energy accounting: the constraint defect is the numerical proxy for the exotic matter's gravitational footprint, and its conversion into a black-hole remnant is the dynamical consequence of removing that support. A fully rigorous treatment would require evolving the exotic matter explicitly (either as a phantom scalar field or a phenomenological fluid), thereby satisfying the Hamiltonian and momentum constraints at $t{=}0$ and tracking the energy budget through the matter sector. Such an extension is an important target for future work; here we instead exploit the sudden approximation to map the qualitative dynamical phase space---collapse vs.\ no collapse, silent vs.\ radiating---with minimal implementation overhead, and use the unperturbed run as a built-in null test to discriminate constraint artifacts from physical signals.

In the ADM convention used throughout this paper, the trace of the extrinsic curvature $K$ encodes the local rate of change of spatial volume: $K>0$ corresponds to volumetric contraction (compression) and $K<0$ to volumetric expansion. A nonzero $K$ at $t{=}0$ means the spatial geometry is already expanding or contracting at the instant we begin the evolution, modeling the residual velocity field inherited from the cosmological environment.

When evolved with no deformation ($A_0{=}A_2{=}0$), the unsupported geometry spontaneously collapses due to the inherent linear instability of the Ellis--Bronnikov solution (Sec.~IV\,A). To study how an initial compressive deformation modifies the collapse dynamics and whether it produces a near-zone curvature signal consistent with gravitational radiation, we introduce a localized, tunable initial deformation ($K>0$) at $t{=}0$:
\begin{equation}
    K(\bar{r}, \theta, \phi) = \left[A_0 + A_2 Y_{20}(\theta, \phi)\right]
    \exp\left(-\frac{\bar{r}^2}{\sigma^2}\right).
\label{eq:kick_profile}
\end{equation}
The monopolar amplitude $A_0$ controls the overall coherent compressive strength (a breathing mode), analogous to an initial implosive velocity. A positive $A_0$ accelerates the collapse by supplying additional inward momentum that reinforces the intrinsic instability.

The quadrupolar amplitude $A_2$ explicitly breaks the exact spherical symmetry of the initial data, modeling a tidal deformation. As spherically symmetric spacetimes cannot radiate, this $Y_{20}$ deformation is required to produce a physical $\Psi_4$ signal. Physically, perfect spherical symmetry is unrealistic for any post-inflationary relic, so the quadrupolar component models the asphericity that tidal interactions and frozen-out vacuum fluctuations would naturally impart.

\subsection{Apparent Horizon Detection}

During the period in which this work was carried out, \texttt{GRTeclyn} did not yet include a production-ready full apparent-horizon finder for these runs. To distinguish genuine topological pinch-off and black-hole-like self-trapping from gauge effects such as lapse collapse, we therefore implemented a geometric trapped-surface proxy on coordinate spheres centered on the throat. At each time step, we scan outward from the origin and estimate the expansion of outgoing null rays, \(\theta_+\). For a conformally flat spatial metric \(\gamma_{ij}=\chi^{-1}\delta_{ij}\), the spherical expression used in the code is
\begin{equation}
    \theta_+ = \frac{2\sqrt{\chi}}{r} - \frac{\partial_r \chi}{\sqrt{\chi}} + \tilde{A}_{rr} - \frac{2}{3}K.
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

\subsection{Collapse diagnostics: physical meaning and computation}
\label{sec:diagnostics_method}

The nine-panel diagnostic figure tracks the geometric state of the wormhole throat at each time step. Each quantity is evaluated on the finest available AMR level to ensure the extreme gradients near the throat are fully resolved.

\noindent\textbf{(a) Minimum lapse \(\boldsymbol{\min(\alpha)}\).}
In the 3+1 decomposition, the lapse function \(\alpha\) measures the rate at which proper time \(\tau\) elapses relative to coordinate time \(t\) for an Eulerian (normal) observer: \(d\tau = \alpha\,dt\). A plunging lapse \(\alpha \to 0\) indicates that proper time is effectively freezing at that location---the hallmark of black-hole formation in the moving-puncture gauge. A lapse that remains nonzero confirms the absence of a singularity (or, more precisely, the trumpet-slice equilibrium in which the lapse asymptotes to a small constant at the puncture). We track \(\min_{\mathbf{x}}[\alpha(t,\mathbf{x})]\) over the computational domain.

\noindent\textbf{(b) Minimum conformal factor \(\boldsymbol{\min(\chi)}\).}
In the BSSN/CCZ4 conformal decomposition, the physical spatial metric is \(\gamma_{ij} = \chi^{-1}\tilde\gamma_{ij}\). Consequently, \(\chi \to 0\) means the physical distances diverge---the spatial geometry is being infinitely stretched, which corresponds to the interior of a black hole in the moving-puncture representation. The conformal factor is related to the conformal weight by \(\chi = \psi^{-4}\), where \(\psi\) is the standard conformal factor. Tracking \(\min(\chi)\) reveals when and how the throat interior stretches toward the trumpet geometry.

\noindent\textbf{(c) Maximum curvature \(\boldsymbol{\max|K|}\).}
The trace of the extrinsic curvature \(K = \gamma^{ij} K_{ij}\) encodes the local rate of change of the spatial volume element:
\begin{equation}
    \partial_t (\sqrt{\gamma}) = -\alpha K \sqrt{\gamma} + \nabla_i \beta^i \sqrt{\gamma}.
    \label{eq:K_volume}
\end{equation}
In the ADM sign convention used here, \(K>0\) means the spatial volume is contracting (compression) and \(K<0\) means it is expanding. A \(\max|K|\) that spikes and then settles to a constant indicates the geometry has reached a stationary state---for a collapsed black hole in the moving-puncture gauge, this constant is the trumpet-slice equilibrium value. An exponentially decaying \(\max|K|\) would instead indicate relaxation toward flat space.

\noindent\textbf{(d) Trapped-surface radius proxy \(\boldsymbol{r_{\rm AH}}\).}
This is the maximum coordinate radius at which the outgoing null expansion \(\theta_+\) is non-positive (Eq.~8). A nonzero \(r_{\rm AH}\) indicates that outgoing light rays are converging---a trapped surface has formed and the topology has pinched off. In the moving-puncture gauge, \(r_{\rm AH}\) typically shrinks over time as the coordinate system stretches the interior, even though the physical apparent horizon area may be growing; the diagnostic therefore tracks the coordinate representation of the trapped region.

\noindent\textbf{(e) Minimum null expansion proxy \(\boldsymbol{\min(\theta_+)}\).}
The expansion of outgoing null geodesics, \(\theta_+\), measures whether a congruence of outgoing light rays is diverging (\(\theta_+>0\)) or converging (\(\theta_+<0\)). On a marginally trapped surface, \(\theta_+=0\). A dramatic plunge of \(\min(\theta_+)\) to large negative values signals the onset of strong gravitational trapping---the geometry is compressing outgoing light cones so severely that not even photons can escape.

\noindent\textbf{(f) Coordinate radius at \(\boldsymbol{\min(\theta_+)}\).}
This reports where on the grid the strongest trapping occurs. In a moving-puncture evolution, this radius settles at the puncture location (\(\bar r\approx 0\)), confirming that the extremal trapping is associated with the black-hole interior rather than a finite-radius feature.

\noindent\textbf{(g) Throat areal radius \(\boldsymbol{R_{\rm areal,min}}\).}
The areal radius is a coordinate-independent measure of the physical size of a sphere. For our conformally flat spatial metric, it is computed as
\begin{equation}
    R_{\rm areal}(\bar r) = \frac{\bar r}{\sqrt{\chi(\bar r)}},
    \label{eq:R_areal}
\end{equation}
where \(\bar r\) is the isotropic coordinate radius. At each time step, we locate the minimum of \(R_{\rm areal}\) along the \(x\)-axis, which corresponds to the physical size of the wormhole throat. A decreasing \(R_{\rm areal,min}\) indicates physical contraction of the throat; a late-time plateau indicates the trumpet radius of the remnant black hole. This is the most direct, gauge-independent indicator of collapse.

\noindent\textbf{(h) Expansion velocity \(\boldsymbol{dR_{\rm areal}/dt}\).}
The time derivative of the throat areal radius gives the coordinate velocity of the throat contraction or expansion, expressed in units of \(c\). Superluminal values (\(|dR_{\rm areal}/dt| > c\)) during early transients are coordinate-velocity artifacts of the gauge adjustment (the proper velocity of the throat, measured by local observers, remains subluminal). The velocity approaching zero at late times confirms a static remnant.

\noindent\textbf{(i) \(\boldsymbol{|K|}\) decay and lifetime \(\boldsymbol{\tau}\).}
To quantify how rapidly the geometry settles after the initial transient, we fit an exponential decay plus constant:
\begin{equation}
    \max|K|(t) = A\,e^{-t/\tau} + C
    \label{eq:K_decay_fit}
\end{equation}
to the late-time portion of the \(\max|K|\) curve. A large \(\tau\) indicates slow relaxation toward flat space; \(\tau\approx 0\) with a nonzero \(C\) indicates the system has settled into a stationary collapsed state (the trumpet slice) rather than dispersing. The offset \(C\) then gives the equilibrium value of \(|K|\) on the trumpet.

\subsection{Gravitational-wave extraction: methodology and observational connection}
\label{sec:gw_method}

The six-panel gravitational-wave analysis figure connects the numerical simulation to observable gravitational-wave physics. The extraction pipeline begins with the Weyl scalar and ends with a detectability estimate against Advanced LIGO.

\noindent\textbf{(a,\,b) Radius-scaled waveform \(\boldsymbol{r\,\mathrm{Re}(\Psi_4^{2,0})}\).}
The Weyl scalar \(\Psi_4\) encodes the outgoing transverse-traceless curvature perturbation, and in the wave zone satisfies \(\Psi_4 = \ddot{h}_+ - i\ddot{h}_\times\), where \(h_{+,\times}\) are the two polarizations of the gravitational-wave strain. We decompose \(\Psi_4\) into spin-weighted spherical harmonics; the \((\ell,m)=(2,0)\) mode dominates for our axisymmetric deformation. Because \(\Psi_4\) falls off as \(1/r\) in the wave zone, multiplying by the extraction radius \(r\) yields a quantity that is approximately radius-independent for genuine gravitational waves.

Panel~(a) plots this radius-scaled quantity in simulation time \(t\), showing the raw signal arrival at each extraction radius. Panel~(b) shifts to \emph{retarded time} \(u = t - R_{\rm ext}\), which removes the light-travel delay between the source and the extraction sphere. For a genuine outgoing wavefront, the waveforms at all extraction radii should align in retarded time; misalignment indicates either near-zone effects or non-radiative (gauge/constraint) modes.

\noindent\textbf{(c) Energy Spectral Density (ESD) of \(\boldsymbol{\Psi_4}\).}
Because the extracted signal is a short, single-burst transient rather than continuous stationary noise, we compute the one-sided energy spectral density using a Tukey-windowed discrete Fourier transform (with taper parameter \(\alpha=0.25\)) rather than a segment-averaged periodogram.  Both gravitational-wave polarizations are retained:
\begin{equation}
    S_{\Psi_4}(f) = \frac{\Delta t^2}{T}\left(\left|\widetilde{\mathrm{Re}(\Psi_4)}\right|^2 + \left|\widetilde{\mathrm{Im}(\Psi_4)}\right|^2\right),
    \label{eq:esd_psi4}
\end{equation}
where \(\widetilde{\cdot}\) denotes the DFT.  Peaks in the ESD identify the characteristic oscillation frequencies of the spacetime.  For a black-hole remnant, these correspond to quasi-normal modes (QNMs); a monotonically falling spectrum with no pronounced peak indicates a broadband transient rather than a ringdown.

\noindent\textbf{(d) Propagation speed analysis.}
This panel provides the critical test for distinguishing physical gravitational radiation from numerical artifacts. The procedure is:
\begin{enumerate}
    \item For each extraction radius \(R_i\), find the time \(t_i^{\rm peak}\) of the dominant peak in the envelope \(|r\,\Psi_4|\).
    \item Compute the propagation speed between successive radii:
    \begin{equation}
        v = \frac{R_{i+1} - R_i}{t_{i+1}^{\rm peak} - t_i^{\rm peak}}.
        \label{eq:propagation_speed}
    \end{equation}
\end{enumerate}
In geometrized units (\(G=c=1\)), physical gravitational waves propagate at exactly \(v=1\). CCZ4 constraint-damping modes, by contrast, propagate at speeds set by the numerical damping parameters---typically superluminal (\(v>1\)). A measurement of \(v\approx 1.0\) within the numerical accuracy of peak identification therefore constitutes strong evidence for a physical gravitational wave, while \(v\gg 1\) identifies the signal as a constraint-cleaning artifact. The red dots on the plot mark the identified peaks.

\noindent\textbf{(e) Strain PSD (code units).}
Gravitational-wave detectors measure strain \(h\), not \(\Psi_4\). Since \(\Psi_4 = \ddot h\) in the frequency domain, the strain PSD is obtained by dividing by \((2\pi f)^4\):
\begin{equation}
    S_h(f) = \frac{S_{\Psi_4}(f)}{(2\pi f)^4}.
    \label{eq:strain_psd}
\end{equation}
Because numerical relativity data always contain low-frequency drift, dividing by \(f^4\) as \(f\to 0\) would cause unphysical divergence.  We therefore apply a 4th-order Butterworth-style high-pass roll-off below \(f_{\rm low} = 0.05\,f_{\rm max}\) to suppress this artifact.  The result is in code units (\(G=c=M=1\)), with dimensions of \(M^3\).

\noindent\textbf{(f) Characteristic strain vs.\ Advanced LIGO sensitivity.}
To connect the simulation to actual detectors, the code-unit strain PSD must be rescaled to physical units by specifying a total mass \(M\) and luminosity distance \(D\):
\begin{align}
    f_{\rm phys} &= \frac{f_{\rm code}}{M \cdot M_{\odot,\rm sec}}, \label{eq:f_phys}\\
    S_h^{\rm phys}(f) &= S_h^{\rm code}(f) \cdot M_{\odot,\rm sec} \cdot \left(\frac{M\cdot M_{\odot,\rm m}}{D}\right)^2, \label{eq:Sh_phys}
\end{align}
where \(M_{\odot,\rm sec} = GM_\odot/c^3 \approx 4.93\times 10^{-6}\)~s and \(M_{\odot,\rm m} = GM_\odot/c^2 \approx 1477\)~m. The characteristic strain \(h_{\rm char}(f) = \sqrt{f\,S_h(f)}\) is then overlaid on the Advanced LIGO design sensitivity noise PSD \(S_n(f)\)~\cite{ajith11}. Where the signal curve lies above the noise curve, the signal is in principle detectable. The optimal matched-filter signal-to-noise ratio is
\begin{equation}
    \mathrm{SNR}^2 = 4\int_0^\infty \frac{|{\tilde h}(f)|^2}{S_n(f)}\,df.
    \label{eq:snr}
\end{equation}
An SNR~$\gtrsim 8$ is conventionally required for confident detection.


\section{Numerical Setup}

Numerical evolutions were performed using the \texttt{GRTeclyn} codebase, which is currently under active development. Built upon the \texttt{AMReX} framework~\cite{amrex}, \texttt{GRTeclyn} implements block-structured Adaptive Mesh Refinement (AMR) and is explicitly optimized for high-performance GPU acceleration. Production runs targeted H100-class GPUs (corresponding to \texttt{CUDA\_ARCH=90} in the build configuration) and used one MPI rank per GPU, scaling up to 8 GPUs for the largest calculations. In practice, the simulation workflow evolves the bare geometry directly under the sudden approximation (Sec.~II\,D), and diagnoses whether it disperses, self-traps, or converts its initial constraint defect into propagating transients.

The full physical computational domain spans a coordinate length of \(L_{\text{full}} = 40\), covered by a coarse grid of \(N_{\text{full}} = 160\) cells, yielding a coarse resolution of \(dx_{\text{coarse}} = 0.25\). To drastically reduce the computational cost, we exploit the Cartesian reflection symmetries of the single-throat initial data and evolve only the positive octant (\(x \ge 0, y \ge 0, z \ge 0\)), effectively modeling just one eighth of the physical volume. Parity symmetry conditions are imposed at the inner reflection boundaries, while Sommerfeld radiative boundary conditions are enforced at the outer edges.

To capture the extreme curvature gradients developing during the collapse while maintaining tractability, we employ up to 6 levels of 2:1 adaptive mesh refinement (AMR). The mesh is dynamically regridded based on the gradients of the conformal factor ($\chi$), allowing the code to focus resolution precisely where the wormhole throat is pinching off. The finest grid spacing ranges from $dx_{\rm fine}\approx 7.8\times10^{-3}$ (5 levels) to $dx_{\rm fine}\approx 3.9\times10^{-3}$ (6 levels), depending on the run. Time stepping is handled using a 4th-order Runge-Kutta scheme, with a conservative Courant factor (dt\_multiplier) reduced to 0.05 to maintain stability during the constraint relaxation and the abrupt collapse branches.

We focus on two representative evolutions---unperturbed and perturbed (compressive deformation)---whose full parameter sets are given at the beginning of Sec.~IV.


\section{Results}

All evolutions of the unsupported Ellis--Bronnikov geometry reveal a universal outcome: regardless of whether an initial compressive deformation is applied, the wormhole collapses into a black hole. The deformation modulates the collapse timescale and the character of the accompanying curvature signal, but does not alter the final state. Below we summarize the two representative configurations and the key numerical choices that differ between them.

\noindent\textbf{Diagnostic sampling.}
All collapse diagnostics---$\min(\alpha)$, $\min(\chi)$, $\max|K|$, the trapped-surface proxy $r_{\rm AH}$, and $\min(\theta_+)$---are evaluated on the \emph{finest available AMR level} at each coarse time step. This ensures that the extreme gradients developing in the throat region are fully resolved in the reported quantities, rather than being smoothed by the coarse base grid. Note that $\max|K|$ reports the absolute value; the signed spatial structure of $K$ (distinguishing compressive $K>0$ from expansive $K<0$ regions) is shown separately in the evolution snapshots.

\noindent\textbf{Perturbed configuration (compressive deformation).}
For the collapse with an applied compressive deformation we use the same throat radius as the unperturbed case, $b_0=0.5$, with a moderate monopolar compressive deformation $A_0=2.0$, a small quadrupolar component $A_2=0.05$ (to break spherical symmetry and allow $\ell=2$ radiation), and width $\sigma_K=0.5$. The initial lapse is $\alpha(t{=}0)=\sqrt{\chi}$ (Eq.~\ref{eq:lapse_precollapsed}), and the same conformal-factor floor $\chi_{\rm min}=10^{-8}$ is used. We employ $\kappa_1=3.0$, $\sigma_{\rm KO}=1.0$, $L_{\rm full}=40$, $N_{\rm full}=160$, and 5 levels of AMR ($\text{max\_level}=5$, $dx_{\rm fine}\approx 7.8\times10^{-3}$). Using identical gauge, floor, and grid settings for both the perturbed and unperturbed configurations ensures a clean comparison of the physics.

\noindent\textbf{Unperturbed configuration.}
To isolate the natural vacuum fate of the unsupported geometry we also evolve a run with $A_0=A_2=0$ (identically vanishing extrinsic curvature at $t{=}0$, meaning the initial spatial slice is momentarily static) and throat radius $b_0=0.5$. The evolution is therefore driven entirely by the Hamiltonian constraint defect inherent in the unsupported wormhole. The initial lapse is $\alpha(t{=}0)=\sqrt{\chi}$ (Eq.~\ref{eq:lapse_precollapsed}), and the conformal-factor floor is set to $\chi_{\rm min}=10^{-8}$, low enough to avoid artificially masking the collapse dynamics near the compactified second asymptotic end at $\bar r=0$. We use $\kappa_1=3.0$, $\sigma_{\rm KO}=1.0$, $L_{\rm full}=40$, $N_{\rm full}=160$, and 5 levels of AMR ($\text{max\_level}=5$, $dx_{\rm fine}\approx 7.8\times10^{-3}$).

\subsection{Unperturbed evolution: spontaneous collapse of the bare topology}

The most fundamental question about an unsupported wormhole is what happens to the bare geometry when no initial deformation is applied at all. Using the unperturbed configuration defined above, we evolve the Ellis--Bronnikov initial data, which is momentarily static and driven entirely by the inherent Hamiltonian constraint defect.


The diagnostics in Fig.~\ref{fig:unperturbed_diagnostics} reveal a geometry that spontaneously \emph{collapses}. The most direct indicator is the throat areal radius (panel~g), which falls from the initial \(R_{\rm areal}=0.50\) to a plateau of \(\sim 0.10\) within \(t\approx 5\). This contraction occurs in two phases: (i)~a violent gauge/constraint transient (\(t\lesssim 2\)), during which the lapse plunges, the null expansion proxy dives to \(\sim -250\), and the throat contracts at coordinate speeds approaching \(c\); and (ii)~a slower, monotonic contraction that stabilizes into a static remnant by \(t\approx 8\). The remnant state exhibits all the hallmarks of a trumpet-sliced black hole: the lapse settles at a small but nonzero value (\(\sim 4\times 10^{-2}\)), \(\max|K|\) converges to a constant (\(\sim 2.5\)), and the areal radius asymptotes to a fixed plateau.

The physical mechanism driving this collapse is the linear instability of the Ellis--Bronnikov wormhole identified by Gonzalez~et~al.~\cite{gonzalez09} and confirmed in spherical symmetry by Shinkai and Hayward~\cite{shinkai02}. Removing the exotic support field leaves the throat unsupported against perturbations, and the truncation error of the numerical scheme provides a seed perturbation that is sufficient to trigger the unstable mode. The choice of \(\chi_{\rm min}=10^{-8}\) is critical: with a higher floor (\(\chi_{\rm min}=10^{-4}\)), the artificial plateau near \(\bar r=0\) effectively stabilizes the inner region, masking the physical instability and producing an apparent expansion that is in fact an artifact of the protective floor. Only when the floor is lowered sufficiently does the true vacuum fate---collapse---emerge.


As theoretically expected from spherical symmetry, the unperturbed run should yield \(\Psi_4=0\) in the continuum limit. The small signal extracted in Fig.~\ref{fig:psi4_unperturbed} (amplitude \(\sim 4\times 10^{-4}\)) represents the numerical noise floor of our setup, driven by the Cartesian grid and AMR refinement boundaries breaking exact spherical symmetry during the violent constraint-damping phase. Crucially, the propagation speed analysis (panel~d) confirms that this signal is entirely non-physical: the dominant peaks travel at \(v\approx 1.4\)--\(1.5c\), well above the speed of light, identifying them as CCZ4 constraint-damping modes rather than gravitational waves.

\noindent\textbf{Absence of ringdown despite black-hole formation.}
Although the diagnostics in Fig.~\ref{fig:unperturbed_diagnostics} unambiguously confirm that a black hole has formed---the lapse freezes, \(\theta_+\) becomes negative, and the areal radius settles onto a trumpet plateau---the extracted \(\Psi_4\) in Fig.~\ref{fig:psi4_unperturbed} contains \emph{no quasi-normal-mode ringdown}. The ESD (panel~c) falls off monotonically with no peaked spectral feature, and the time-domain waveform (panels~a,\,b) shows erratic, noise-like oscillations rather than a coherent, exponentially decaying sinusoid. This is the expected consequence of Birkhoff's theorem applied to the remnant itself: because the collapse proceeds with exact spherical symmetry (to machine precision), the black hole forms in a perfectly round state with no \(\ell=2\) deformation. A QNM is a resonance that must be \emph{excited} by a non-spherical perturbation; when none is present, the black hole simply sits silently. This result provides a critical null test for the perturbed analysis in Sec.~IV\,C: the oscillatory late-time structure observed in the perturbed \(\Psi_4\) waveform (Fig.~\ref{fig:psi_perturbed}) is entirely absent here, confirming that it is driven by the \(A_2\) quadrupolar deformation rather than being an intrinsic feature of the collapse dynamics. The two-order-of-magnitude amplitude ratio (\(\sim 3\times 10^{-2}\) vs.\ \(\sim 4\times 10^{-4}\)) further establishes that the perturbed signal sits well above the numerical noise floor.

The constraint norms (Appendix, Fig.~\ref{fig:constraints_unperturbed}) remain bounded throughout, confirming that the CCZ4 evolution is well-controlled.


\subsection{Perturbed evolution: compressive deformation and accelerated collapse}

When an initial compressive deformation is applied, the collapse is dramatically accelerated and accompanied by a near-zone curvature signal that is absent in the unperturbed case.

Using the perturbed configuration defined above, we apply a compressive deformation to explicitly break spherical symmetry and enable \(\ell=2\) curvature radiation.


The diagnostics in Fig.~\ref{fig:perturbed_diagnostics} confirm that the compressive deformation drives black-hole formation on the same timescale as the unperturbed case but with substantially stronger initial trapping. The most direct comparison is the throat areal radius (panel~g): it falls from \(R_{\rm areal}=0.50\) to the same \(\sim 0.10\) trumpet-radius plateau observed in the unperturbed run, but the initial contraction is noticeably steeper. The null expansion proxy (panel~e) plunges to \(\sim -700\), nearly three times deeper than the unperturbed value of \(\sim -250\), reflecting the additional compressive momentum. At late times, however, both runs converge to the same equilibrium: \(\min(\alpha)\approx 6\times 10^{-2}\), \(\max|K|\approx 2.5\), and \(R_{\rm areal,min}\approx 0.10\)---confirming that the deformation modulates the collapse dynamics but not the final remnant.

\subsection{Gravitational-wave signatures (near-zone)}

The extracted Weyl signal differs sharply between the two configurations. As established in Sec.~IV\,A, the unperturbed signal represents purely numerical noise and propagates superluminally.

\begin{figure*}[tp]
\centering
\includegraphics[width=0.92\textwidth]{psi4_analysis_perturbed.pdf}
\caption{Gravitational-wave analysis for the perturbed evolution (\(A_0{=}2.0\), \(A_2{=}0.05\), \(b_0{=}0.5\)).
\textbf{(a)}~Waveform \(r\,\mathrm{Re}(\Psi_4^{2,0})\) in simulation time: the signal amplitude (\(\sim 3\times 10^{-2}\)) is two orders of magnitude above the unperturbed noise floor (\(\sim 4\times 10^{-4}\)), confirming a physical origin. The waveforms at all three extraction radii (\(R=8,12,16\)) show a coherent outgoing burst followed by oscillatory structure.
\textbf{(b)}~Waveform in retarded time \(t - R_{\rm ext}\) with QNM fit overlay (red dashed): the waveforms at different extraction radii align in retarded time, consistent with a coherent outgoing wavefront. The late-time tail is fit by a damped sinusoid \(Ae^{-t/\tau}\sin(2\pi f_{\rm QNM}t+\phi)\), yielding \(f_{\rm QNM}\approx 0.53\,M^{-1}\) and \(\tau\approx 1.35\,M\) (quality factor \(Q\approx 2.2\)). This oscillatory structure is entirely absent in the unperturbed case (Fig.~\ref{fig:psi4_unperturbed}b, where the fitter produces a flat line), confirming it is driven by the \(A_2\) quadrupolar deformation.
\textbf{(c)}~Energy spectral density of \(\Psi_4\): the ESD shows broadband power peaking near \(f\sim 1\,M^{-1}\), in contrast to the monotonically falling ESD of the unperturbed run.
\textbf{(d)}~Propagation speed analysis: the red dots mark the time-of-arrival of the dominant envelope peak \(|r\Psi_4|_{\max}\) at each extraction radius. Both measurements are consistent with the speed of light: \(v\approx 1.01c\) (\(R{=}8\to 12\)) and \(v\approx 1.00c\) (\(R{=}12\to 16\)). This contrasts sharply with the unperturbed run (Fig.~\ref{fig:psi4_unperturbed}d), where the same analysis yields superluminal speeds (\(v\approx 1.4\text{--}1.5c\)), identifying those signals as constraint-damping modes. The luminal propagation here is consistent with a physical curvature transient rather than a numerical artifact.
\textbf{(e)}~Strain PSD in code units (\(R{=}16\)): peaks near \(f\sim 1\,M^{-1}\) with a clear spectral feature, unlike the featureless unperturbed spectrum.
\textbf{(f)}~Characteristic strain vs.\ Advanced LIGO design sensitivity: for a hypothetical \(M{=}30\,M_\odot\) source at \(D{=}10\)~Mpc the signal peaks at \(\sim 5700\)~Hz, above the most sensitive LIGO band, yielding SNR\(\approx 0\). Detection would require a much more massive progenitor to shift the signal into the LIGO band.}
\label{fig:psi_perturbed}
\end{figure*}

In the perturbed case (Fig.~\ref{fig:psi_perturbed}), the extracted \(\Psi_4\) amplitude is substantially larger, and the propagation speed across all three extraction radii is consistent with the speed of light. This is the central observational result of the perturbed evolution: the compressive deformation not only accelerates the collapse but also produces a near-zone curvature transient that travels at the physically expected speed, in stark contrast to the superluminal modes that dominate the unperturbed signal. Because both runs share the same constraint defect, numerical settings, and remnant type, features appearing only in the perturbed case cannot be attributed to the constraint violation alone.

\noindent\textbf{Ringdown structure and the role of the quadrupolar deformation.}
The retarded-time waveform (Fig.~\ref{fig:psi_perturbed}b) exhibits a clear oscillatory late-time tail, which we quantitatively characterize via a damped-sinusoid fit in Sec.~IV\,D. This oscillatory structure is \emph{completely absent} in the unperturbed case. The contrast is definitive: because both runs form the same type of black-hole remnant (same \(\max|K|\) plateau, same areal radius, same lapse equilibrium), the presence of oscillatory structure only when \(A_2\neq 0\) proves it is seeded by the quadrupolar symmetry breaking. A perfectly spherical black hole has no \(\ell=2\) deformation and therefore no quasi-normal mode to ring down.

\noindent\textbf{Caveats.}
An important caveat applies. The initial extrinsic-curvature deformation (Eq.~\ref{eq:kick_profile}) modifies \(K\) without simultaneously solving for a consistent traceless part \(\tilde A_{ij}\), thereby introducing a momentum-constraint violation in addition to the Hamiltonian defect. The CCZ4 evolution damps these violations by propagating them outward at approximately the speed of light, and these constraint-cleaning modes inevitably couple to the Weyl scalar. Consequently, the extracted \(\Psi_4\) contains a superposition of any physical curvature radiation and constraint-driven transients. A definitive separation of physical radiation from constraint contamination would require constraint-satisfying initial data---an important target for future work.

Furthermore, all extractions are performed in the near zone (\(R_{\rm ext}=8\text{--}16\) with \(L_{\rm full}=40\)); the signals reported here should therefore be understood as near-zone curvature diagnostics rather than asymptotic gravitational-wave templates. The constraint norms for this run (Appendix, Fig.~\ref{fig:constraints_perturbed}) remain bounded throughout, confirming numerical stability. The spatial structure of the \(K\) field during the perturbed evolution is shown in Fig.~\ref{fig:perturbed_kz_evolution} (Appendix).


\subsection{Black-hole remnant characterisation}

The universal collapse outcome permits several quantitative characterisations of the newly formed black hole.

\noindent\textbf{Irreducible mass from the throat areal radius.}
The gauge-independent areal radius (Eq.~\ref{eq:R_areal}) plateaus at late times at \(R_{\rm areal,min}\approx 0.10\) in the unperturbed run. In the moving-puncture 1+log gauge the minimum areal radius corresponds to the trumpet throat. Using the standard relation \(R_{\rm trumpet} \approx 1.31\,M_{\rm BH}\)~\cite{hannam07}, this yields an estimated remnant mass \(M_{\rm BH}\approx 0.077\) in code units. The lapse settles at \(\alpha_{\rm min}\approx 0.034\), consistent with the trumpet-slice puncture value for a BH of this mass.

\noindent\textbf{Total radiated energy.}
Integrating \(|r\Psi_4|^2/(16\pi)\) over retarded time at the innermost extraction radius (\(R=8\)) gives \(E_{\rm rad}\approx 1.6\times 10^{-5}\,M\)---an extremely small fraction of the system energy, consistent with the weakness of the quadrupolar deformation (\(A_2=0.05\)).

\noindent\textbf{QNM ringdown fit.}
To test whether the late-time oscillation in the perturbed \(\Psi_4\) waveform is a black-hole quasi-normal mode, we fit a damped sinusoid \(A e^{-t/\tau}\sin(2\pi f_{\rm QNM} t + \phi)\) to the ringdown tail at \(R=8\). The fit yields \(f_{\rm QNM}\approx 0.53\,M^{-1}\) and \(\tau\approx 1.35\,M\), giving a quality factor \(Q = \pi f\tau \approx 2.2\). Interpreting this frequency as the fundamental Schwarzschild \(\ell=2\) QNM (\(\omega M = 0.3737\)) gives \(M_{\rm BH}^{\rm QNM} = 0.3737/(2\pi f_{\rm QNM}) \approx 0.11\,M\). As a consistency check, the imaginary part of the same mode (\(\omega_I M = 0.0890\)) predicts \(\tau_{\rm pred} = M_{\rm BH}/0.0890 \approx 1.26\,M\), which agrees with the observed damping time to within \(\sim 7\%\). This internal consistency between the frequency- and damping-derived mass estimates supports the identification of the late-time signal with black-hole ringdown, though contamination from CCZ4 constraint-cleaning modes is expected to broaden the residual scatter.

\noindent\textbf{Mass discrepancy between geometric and QNM estimates.}
The QNM-inferred mass (\(M_{\rm BH}^{\rm QNM}\approx 0.11\)) exceeds the trumpet-radius geometric mass (\(M_{\rm BH}^{\rm geom}\approx 0.077\)) by roughly \(45\%\). This discrepancy is expected given the contaminated background: the trumpet relation \(R_{\rm trumpet}\approx 1.31\,M\)~\cite{hannam07} assumes a pristine Schwarzschild spacetime in isotropic coordinates, whereas our background carries a persistent Hamiltonian constraint defect of order \(\sim 3\times 10^{-2}\) (effectively a residual phantom-energy bath from the removed exotic matter). Curvature waves propagating through this contaminated near-zone experience a modified effective potential, shifting the QNM frequencies and damping times away from the analytic predictions. Given that neither estimate is performed on a constraint-satisfying background, the \(\sim 45\%\) disagreement reflects the level of nonlinear contamination rather than a methodological error; both estimates independently confirm a compact remnant of order \(M_{\rm BH}\sim 0.1\) in code units.

\noindent An important cross-check is provided by the unperturbed run: despite forming an identical type of black-hole remnant (Sec.~IV\,A), its \(\Psi_4\) waveform contains \emph{no} oscillatory ringdown structure whatsoever---only noise-level, superluminally propagating constraint modes. The automated QNM fitter, applied to the unperturbed waveform, produces a flat line at zero amplitude (\(\tau\approx 0.02\,M\), \(Q\approx 0.02\))---it finds nothing to lock onto. Because the unperturbed collapse is spherically symmetric, the resulting black hole has no \(\ell=2\) deformation and therefore no QNM to ring down. The fact that the late-time oscillatory structure appears \emph{only} when a quadrupolar deformation is present (\(A_2\neq 0\)) confirms that the fitted oscillation is seeded by the initial \(A_2 Y_{20}\) symmetry breaking. However, the characteristic frequency of the initial Gaussian perturbation (set by \(\sigma_K\)) is comparable to the fitted QNM frequency for the parameters used here; a definitive separation of the genuine QNM ringdown from the spectral content of the initial perturbation would require repeating the simulation with a significantly different \(\sigma_K\) while keeping the remnant mass fixed. If the late-time frequency remains unchanged under such a variation, it is a genuine QNM; if it tracks \(\sigma_K\), it reflects the initial-data spectral content.

\noindent\textbf{Instability growth rate.}
The early-time departure of the throat areal radius from its initial value in the unperturbed run is well fit by an exponential: \(|R_{\rm areal}(t) - R_{\rm areal}(0)| \sim e^{\lambda t}\) with \(\lambda \approx 0.50\,M^{-1}\), corresponding to an e-folding time of \(\sim 2.0\,M\). This provides a numerical measurement of the linear instability rate of the unsupported Ellis--Bronnikov geometry, complementing the analytical predictions of Gonzalez~et~al.~\cite{gonzalez09}.

\noindent\textbf{Detectability.}
For a hypothetical \(M=30\,M_\odot\) remnant at \(D=10\)~Mpc, the physical QNM frequency is \(f_{\rm phys}\approx 3500\)~Hz, placing it at the upper edge of the Advanced LIGO band but with a characteristic strain \(h_{\rm char}\sim 10^{-28}\), far below the detector noise floor (\(\mathrm{SNR}\approx 0\)). Detection would require either a much more massive wormhole progenitor---an intermediate-mass wormhole (\(\sim 10^4\,M_\odot\)) would shift the signal into the most sensitive LIGO band (\(\sim 200\)~Hz), while a supermassive wormhole (\(\sim 10^8\,M_\odot\)) would be a target for the space-based LISA observatory in the millihertz band---or a much closer source.



\section{Discussion: 3D Dynamics versus 1D Instability Proofs}

Our 3D numerical results are consistent with the collapse pathway identified by Shinkai and Hayward \cite{shinkai02}, but extend the picture in two important ways. First, by lowering the conformal-factor floor to \(\chi_{\rm min}=10^{-8}\), we demonstrate that the unperturbed geometry spontaneously collapses---a result that was masked in earlier runs with a higher floor (\(\chi_{\rm min}=10^{-4}\)), which artificially stabilized the inner region and produced an apparent expansion. Second, by performing the simulations in full 3D, we can extract the Weyl scalar \(\Psi_4\) and use propagation speed analysis to distinguish near-zone curvature transients consistent with gravitational radiation from constraint-damping artifacts.

Crucially, our 3D approach reveals behavior that 1D codes inherently cannot capture. Because Shinkai and Hayward utilized a 1D spherically symmetric simulation, any nontrivial curvature signal was identically zero. By moving the problem to full 3D and explicitly breaking exact symmetry with a quadrupolar deformation component, we can inspect \(\Psi_4\) directly and employ propagation speed analysis to assess its physical origin. 

Beyond the collapse dynamics, the 3D simulations yield several additional quantitative characterizations of the remnant that were inaccessible in 1D, as detailed in Sec.~IV\,D. We stress that all near-zone curvature quantities reported here should be regarded as diagnostics of the local vacuum response rather than as asymptotic gravitational-wave templates.

\section{Conclusion}
We have investigated the nonlinear vacuum dynamics of unsupported wormhole-like topological defects using full 3D numerical relativity. The central question in this work is not how to sustain a traversable wormhole with exotic matter, but what becomes of the residual geometry once that support is absent.

The central result is that the unsupported Ellis--Bronnikov geometry collapses into a black hole in all cases studied. With no applied initial deformation, the bare topology spontaneously collapses due to the well-known linear instability of the Ellis--Bronnikov solution~\cite{gonzalez09}. When a compressive initial deformation is applied, the collapse is dramatically accelerated. In both cases, the unsupported geometry does not remain as a static traversable throat on the resolved grid. While quantum backreaction may theoretically stabilize microscopic topological defects at the Planck scale \cite{mehulic2026}, our results demonstrate that once stretched to macroscopic scales where classical vacuum dynamics dominate, they lack the sustained exotic support required to survive. This provides a natural cosmological mechanism explaining why inflationary wormholes do not abundantly persist in the late universe.

The extracted near-zone curvature signal provides a powerful diagnostic tool. The contrast between the unperturbed and perturbed runs---identical numerical settings and constraint defect, but dramatically different signals---is the primary evidence that the perturbed curvature transient contains a physical component. As previously noted, the inevitable contamination from constraint-cleaning modes means we characterize the perturbed signal as a near-zone curvature transient consistent with gravitational radiation, rather than claiming a definitive detection of physical gravitational waves.

Several natural extensions would place the results on firmer ground. First, constructing constraint-satisfying initial data---by solving the Hamiltonian and momentum constraints with the desired \(K\) deformation, using an elliptic solver---would eliminate the initial constraint defect and remove the dominant source of contamination in the extracted \(\Psi_4\). Second, evolving the exotic matter explicitly, either as a phantom scalar field (the exact Ellis--Bronnikov source) or as a phenomenological anisotropic fluid with a prescribed decay law (in the spirit of \citealt{clough24}), would provide a fully rigorous energy accounting: the ADM mass would be conserved, the matter's negative energy would be tracked through the fluid sector, and the remnant mass would be a calculable prediction rather than an inferred estimate. Third, extracting \(\Psi_4\) at larger radii (or using Cauchy-characteristic extraction) would move the measurement into the wave zone and enable the construction of asymptotic waveform templates suitable for matched-filter searches. These extensions are the targets of our ongoing development of the \texttt{GRTeclyn} framework.

\section{Acknowledgements}
This research was supported by Gravity Frontiers (\texttt{ic@gravityfrontiers.org}) under a research grant for the development of software for mathematical modeling of wormholes.

The author also thanks Ilya Nachevsky for generously providing the high-performance computing resources necessary for this work. These simulations would not have been possible without his support. 

\clearpage

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
These two scalars are written to the small-data file as \texttt{L2\_Ham} and \texttt{L2\_Mom}. Because the unsupported or partially unsupported evolutions begin from data that are not exactly on the target vacuum constraint surface, these norms are generically nonzero at \(t=0\). The relevant stability test is therefore not whether the initial slice is constraint-free, but whether the CCZ4 evolution keeps the subsequent violation bounded. Figures~\ref{fig:constraints_unperturbed} and~\ref{fig:constraints_perturbed} show that for both the unperturbed and perturbed evolutions, the Hamiltonian and momentum RMS norms exhibit an initial transient and then remain bounded throughout the evolution. We note that the saturated Hamiltonian norm (\(\sim 2.5\times 10^{-2}\)) is large by the standards of constraint-satisfying numerical relativity, reflecting the fact that the grid never fully settles onto the Einstein constraint surface---an expected consequence of the deliberately unphysical initial data and the persistent under-resolution of the compactified second asymptotic end near $\bar r=0$. The boundedness of the norms nevertheless confirms numerical stability and ensures that the CCZ4 evolution does not develop runaway violations that would invalidate the dynamical conclusions.

\begin{figure*}[tp]
\centering
\includegraphics[width=0.92\textwidth]{diagnostic_unperturbed.pdf}
\caption{Diagnostics of the unperturbed evolution (\(A_0{=}A_2{=}0\), \(b_0=0.5\), \(\chi_{\rm min}=10^{-8}\)).
\textbf{(a)}~Minimum lapse \(\min(\alpha)\): the lapse plunges during the initial gauge/constraint transient (\(t\lesssim 2\)), then partially recovers to \(\sim 4\times 10^{-2}\) and stabilizes at a small but nonzero value consistent with trumpet slicing of the newly formed black hole.
\textbf{(b)}~Minimum conformal factor \(\min(\chi)\): drops sharply during the transient, then slowly rises to \(\sim 10^{-3}\) as the moving-puncture gauge stretches the interior.
\textbf{(c)}~Maximum curvature \(\max|K|\): spikes to \(\sim 10\) during constraint cleaning, then settles to a constant \(\sim 2.5\)---the trumpet-slice equilibrium value for the formed black hole. The absence of further decay confirms that the collapsed state is stationary.
\textbf{(d)}~Trapped-surface radius proxy \(r_{\rm AH}\): starts at \(\sim 0.28\) (reflecting the initial wormhole throat topology) and decreases monotonically to \(\sim 0.14\) as the throat pinches off and the coordinate representation of the trapped region shrinks under the moving-puncture gauge.
\textbf{(e)}~Minimum null expansion proxy \(\min(\theta_+)\): the dramatic plunge to \(\sim -250\) at \(t\approx 1.3\) signals the onset of strong self-trapping; subsequent recovery toward zero reflects the gauge settling into the trumpet slice.
\textbf{(f)}~Coordinate radius at \(\min(\theta_+)\): settles at \(\sim 0.02\), indicating the puncture location on the grid.
\textbf{(g)}~Throat areal radius \(R_{\rm areal,min}\): falls from \(0.50\) (the bare wormhole throat) to \(\sim 0.10\) within \(t\approx 5\), confirming physical contraction and topological pinch-off. The late-time plateau corresponds to the trumpet radius of the remnant.
\textbf{(h)}~Expansion velocity \(dR_{\rm areal}/dt\): the brief superluminal spike (\(\sim -0.9c\)) at \(t\lesssim 2\) is a coordinate-velocity artifact of the gauge transient (proper velocity remains subluminal). The velocity subsequently approaches zero, consistent with a static remnant.
\textbf{(i)}~\(|K|\) decay and lifetime \(\tau\): the fit yields \(\tau\approx 0\), reflecting the fact that \(\max|K|\) converges to a nonzero constant rather than decaying---the geometry has collapsed into a black hole whose trumpet slice maintains a fixed \(K\) rather than relaxing toward flat space.}
\label{fig:unperturbed_diagnostics}
\end{figure*}

\begin{figure*}[tp]
\centering
\includegraphics[width=0.92\textwidth]{psi4_analysis_unperturbed.pdf}
\caption{Gravitational-wave analysis for the unperturbed evolution (\(A_0{=}A_2{=}0\)).
\textbf{(a)}~Waveform \(r\,\mathrm{Re}(\Psi_4^{2,0})\) in simulation time: the signal amplitude (\(\sim 4\times 10^{-4}\)) is two orders of magnitude below the deformed runs, consistent with purely numerical noise.
\textbf{(b)}~Waveform in retarded time \(t - R_{\rm ext}\) with QNM fit overlay (red dashed): aligning by retarded time reveals that the waveforms at different extraction radii do not collapse onto a single curve, indicating the absence of a coherent outgoing wavefront. The automated QNM fitter attempts a damped-sinusoid fit to the late-time tail, but the result is a flat line at zero amplitude (\(\tau\approx 0.02\,M\), \(Q\approx 0.02\))---the fitter finds no oscillatory ringdown structure to lock onto, confirming the complete absence of quasi-normal-mode ringing despite the confirmed black-hole formation.
\textbf{(c)}~Energy spectral density of \(\Psi_4\): the ESD at all three extraction radii falls off steeply and monotonically, with no pronounced quasi-normal-mode peak---in stark contrast to the perturbed case, which exhibits spectral structure near the QNM frequency.
\textbf{(d)}~Propagation speed analysis: the red dots mark the time-of-arrival of the dominant envelope peak \(|r\Psi_4|_{\max}\) at each extraction radius; the propagation speed is computed as \(v = \Delta R / \Delta t_{\rm peak}\) between consecutive pairs. Both measurements yield significantly superluminal values (\(v\approx 1.39c\) for \(R{=}8\to 12\) and \(v\approx 1.49c\) for \(R{=}12\to 16\)). Physical gravitational waves travel at \(v=c\); the superluminal speeds confirm that the extracted signal is dominated by CCZ4 constraint-damping modes rather than physical radiation.
\textbf{(e)}~Strain PSD in code units (\(R{=}16\)): monotonically decreasing, with no identifiable spectral feature.
\textbf{(f)}~Characteristic strain vs.\ Advanced LIGO design sensitivity: for a hypothetical \(M{=}30\,M_\odot\) source at \(D{=}10\)~Mpc the signal falls entirely below the detector noise floor, yielding SNR\(\approx 0\). This establishes the unperturbed run as a null baseline for the gravitational-wave analysis of the deformed configurations.}
\label{fig:psi4_unperturbed}
\end{figure*}

\begin{figure*}[tp]
    \centering
    \includegraphics[width=0.92\textwidth]{collapse_diagnostics_perturbed.pdf}
    \caption{Diagnostics of the perturbed collapse (\(A_0{=}2.0\), \(A_2{=}0.05\), \(b_0{=}0.5\), \(\chi_{\rm min}=10^{-8}\)).
\textbf{(a)}~Minimum lapse \(\min(\alpha)\): drops sharply to \(\sim 10^{-9}\) during the initial compressive transient, then recovers and stabilizes at \(\sim 6\times 10^{-2}\)---a small but nonzero value consistent with the trumpet-slice puncture equilibrium of the newly formed black hole.
\textbf{(b)}~Minimum conformal factor \(\min(\chi)\): plunges to the floor (\(10^{-8}\)) during the transient, then slowly rises to \(\sim 10^{-3}\) as the moving-puncture gauge stretches the interior.
\textbf{(c)}~Maximum curvature \(\max|K|\): spikes to \(\sim 10\) at \(t\approx 0\) (reflecting the applied \(A_0=2.0\) deformation), then settles to a constant \(\sim 2.5\)---the same trumpet-slice equilibrium value observed in the unperturbed case, confirming that both runs produce the same type of remnant.
\textbf{(d)}~Trapped-surface radius proxy \(r_{\rm AH}\): starts at \(\sim 0.32\) and decreases to \(\sim 0.15\) as the coordinate representation of the trapped region shrinks under the moving-puncture gauge.
\textbf{(e)}~Minimum null expansion proxy \(\min(\theta_+)\): plunges to \(\sim -700\) at \(t\approx 1\), indicating even stronger self-trapping than the unperturbed case (\(\sim -250\)), driven by the additional compressive momentum.
\textbf{(f)}~Coordinate radius at \(\min(\theta_+)\): settles at \(\sim 0.02\), consistent with the puncture location.
\textbf{(g)}~Throat areal radius \(R_{\rm areal,min}\): falls from \(0.50\) to \(\sim 0.10\) by \(t\approx 7\), matching the unperturbed remnant's trumpet radius plateau. The faster initial contraction compared to the unperturbed case reflects the additional compressive momentum.
\textbf{(h)}~Expansion velocity \(dR_{\rm areal}/dt\): the initial spike exceeds \(-c\) in coordinate velocity during the violent transient. The velocity approaches zero at late times, confirming a static remnant.
\textbf{(i)}~\(|K|\) decay and lifetime: the fit yields \(\tau\approx 0.24\), reflecting the rapid settling of \(\max|K|\) to a constant---the geometry reaches the trumpet equilibrium without relaxing toward flat space.}
    \label{fig:perturbed_diagnostics}
\end{figure*}
\begin{figure}[htp]
\centering
\includegraphics[width=\linewidth]{constraints_unperturbed.pdf}
\caption{Volume-weighted RMS norms of the Hamiltonian and momentum constraints for the unperturbed evolution (\(A_0{=}A_2{=}0\)). The Hamiltonian norm grows slowly but remains at \(\mathcal{O}(10^{-2})\); the momentum norm saturates at \(\sim 8\times 10^{-3}\). These levels are large relative to constraint-satisfying evolutions, reflecting the deliberately unphysical initial data, but the norms remain bounded throughout the \(t\approx 18\) evolution, confirming numerical stability.}
\label{fig:constraints_unperturbed}
\end{figure}

\begin{figure}[htp]
\centering
\includegraphics[width=\linewidth]{constraints_plot_perturbed.pdf}
\caption{Volume-weighted RMS norms of the Hamiltonian and momentum constraints for the perturbed evolution (\(A_0{=}2.0\), \(A_2{=}0.05\)). The Hamiltonian norm rises during the initial transient to \(\sim 3\times 10^{-2}\) and then slowly decays; the momentum norm peaks at \(\sim 1.5\times 10^{-2}\) during the constraint-cleaning phase before settling to \(\sim 9\times 10^{-3}\). Both norms remain bounded throughout, confirming numerical stability.}
\label{fig:constraints_perturbed}
\end{figure}

\begin{figure*}[tp]
\centering
\includegraphics[width=\textwidth]{evolution_panel_perturbed.pdf}
\caption{Evolution snapshots of the \(K\) field in the \(z{=}0\) plane for the perturbed run (\(A_0{=}2.0\), \(A_2{=}0.05\), \(b_0{=}0.5\)). The initial Gaussian compressive deformation (bright central region at \(t{=}0\)) drives rapid infall. By \(t\approx 2\), the throat has pinched off and the \(K\) field develops the trumpet-slice profile: a localized positive-\(K\) core at the puncture surrounded by a negative-\(K\) halo, with outward-propagating ripples representing the constraint-cleaning transient and the near-zone curvature signal. At late times the exterior relaxes toward \(K\to 0\) as required for asymptotic flatness.}
\label{fig:perturbed_kz_evolution}
\end{figure*}

% --- Bibliography ---
\clearpage
\bibliographystyle{apsrev4-2}
\bibliography{references}
\end{document}