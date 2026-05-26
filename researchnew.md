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
We outline a practical baseline for turning \texttt{GRTeclyn} into a closed-loop search engine for exotic spacetime metrics. The starting point is the previous 3D numerical-relativity study of the Ellis--Bronnikov wormhole, where \texttt{GRTeclyn} evolved the full Einstein--phantom-scalar system, resolved the collapse/expansion instability, extracted gravitational waves, and identified a phantom bounce. The next goal is to automate the search over initial data. A Python optimizer proposes compact recipes for CCZ4 fields, \texttt{GRTeclyn} evolves each candidate in C++, and diagnostics decide whether the geometry is stable, useful, or pathological. The target observables are effective faster-than-light transport, faster-than-light communication, portal-like shortcuts, wormhole throats, and new gravitational-wave-producing exotic dynamics. The immediate deliverable is not a claimed working warp drive, but a reproducible numerical-relativity pipeline that tests candidate metrics by nonlinear evolution rather than static ansatz inspection.
\end{abstract}

\maketitle

\section{Introduction}

Most proposed interstellar metrics are built by hand. One writes an ansatz, computes the required stress-energy tensor, and then checks whether it has interesting properties. This is useful, but it is also biased. We only test the geometries that humans know how to write down.

The previous \texttt{GRTeclyn} wormhole work gives a stronger starting point. It showed that the code can evolve an exotic spacetime in full 3D numerical relativity, track its nonlinear instability, identify horizon formation, and extract gravitational radiation. This means we can now ask a larger question: can we search the space of initial data automatically?

The proposed project treats metric discovery as an optimization problem. The optimizer proposes initial data. \texttt{GRTeclyn} evolves the spacetime. The diagnostics measure what happened. Failed runs are useful, because they reveal why a candidate does not work: collapse, dispersion, constraint violation, gauge failure, excessive exotic matter, or horizon formation.

Throughout this draft, ``FTL'' means effective faster-than-light behavior caused by spacetime geometry. It does not mean local motion faster than the local light cone.

\section{Goal}

The goal is to build an automated numerical-relativity laboratory for exotic spacetime metrics. The search targets four classes of behavior:

\begin{enumerate}
    \item \textbf{FTL transport:} a compact passenger region moves effectively faster than a light signal in the reference flat spacetime, while local tidal forces remain controlled.
    \item \textbf{FTL communication:} a null ray or scalar pulse reaches a detector earlier than it would in flat spacetime with the same asymptotic separation.
    \item \textbf{Portals:} two mouth-like regions create a stable shortcut between distant regions.
    \item \textbf{Wormholes:} a throat-like bridge remains open for a measurable time without immediately becoming an unusable horizon.
\end{enumerate}

The first success criterion is not to find a perfect stable metric. The first success criterion is to build the pipeline and classify known and unknown candidates dynamically.

\section{Search Variables}

\texttt{GRTeclyn} evolves the Einstein equations using CCZ4 variables. Therefore the search should operate directly on those fields. The first candidate recipes should generate:

\begin{equation}
    \{\chi,\alpha,\beta^i,K,\tilde A_{ij},\phi,\Pi\}_{t=0}.
\end{equation}

Here $\chi$ controls spatial stretching and throat-like structure, $\alpha$ is the lapse, $\beta^i$ is the shift, $K$ controls local expansion or contraction, $\tilde A_{ij}$ contains trace-free extrinsic curvature, and $\phi,\Pi$ describe the scalar matter sector.

The optimizer should not generate every grid point directly. Instead, it should generate a small parameter vector $\theta$. A decoder turns $\theta$ into smooth initial fields on the AMReX grid.

\section{Candidate Parameterization}

\subsection{Level 1: Radial Recipes}

The first implementation should use simple radial profiles,
\begin{equation}
    q(r) = q_\infty + \sum_n a_n B_n(r),
\end{equation}
where $q$ may be $\chi$, $\alpha$, $\beta^i$, $K$, $\phi$, or $\Pi$. This level is simple, cheap, and sufficient for testing throat-like or bubble-like structures.

\subsection{Level 2: Angular Recipes}

After the radial pipeline works, angular modes can be added:
\begin{equation}
    q(r,\theta,\varphi)
    =
    q_\infty
    + \sum_{n\ell m} a_{n\ell m}B_n(r)Y_{\ell m}(\theta,\varphi).
\end{equation}
This allows asymmetric wormholes, portal-like mouths, rotating structures, and quadrupolar gravitational-wave sources.

\subsection{Level 3: Neural Decoder}

Only after the simple recipes work should the project use a neural field decoder,
\begin{equation}
    (x,y,z)
    \rightarrow
    \{\chi,\alpha,\beta^i,K,\tilde A_{ij},\phi,\Pi\}.
\end{equation}
This is the most flexible option, but it is also the hardest to debug. It should come after the episode loop, diagnostics, and rewards are stable.

\section{Connection to \texttt{GRTeclyn}}

The implementation should keep a clean division of labor:

\begin{itemize}
    \item \textbf{C++ / \texttt{GRTeclyn}:} build initial data, evolve the spacetime, compute diagnostics, and write outputs.
    \item \textbf{Python:} choose parameters, write input files, launch runs, read outputs, compute scores, and update the optimizer.
\end{itemize}

At first, Python can run \texttt{GRTeclyn} as a subprocess. A complicated API is not needed for the baseline.

\subsection{Minimal C++ Additions}

The first useful changes to \texttt{GRTeclyn} are small:

\begin{enumerate}
    \item \textbf{Recipe initial-data class.} A new C++ initial-data module reads coefficients from the input file and fills $\chi$, $\alpha$, $\beta^i$, $K$, $\tilde A_{ij}$, $\phi$, and $\Pi$.
    \item \textbf{Machine-readable diagnostics.} Each run writes constraint norms, minimum lapse, minimum $\chi$, maximum $|K|$, throat radius, horizon proxy, scalar-field extrema, and termination reason.
    \item \textbf{Early stopping.} A run stops if constraints explode, curvature becomes unresolved, lapse behavior is pathological, or the wall-time limit is reached.
    \item \textbf{Probe module.} Null rays, scalar pulses, or timelike tracer particles are evolved to test communication and transport.
    \item \textbf{Optional grid-field loader.} Later, Python can write full grid fields to disk for neural-decoder candidates.
\end{enumerate}

\subsection{Anatomy of a GRTeclyn Example}

Every runnable \texttt{GRTeclyn} problem is an \texttt{Examples/<Name>/} directory with roughly ten tightly coupled files. The Python wrapper can template \texttt{params.txt}, launch the binary, and parse diagnostics, but it cannot invent new physics without these C++ pieces.

\begin{itemize}
    \item \textbf{\texttt{GNUmakefile}.} Build entry point. Selects AMReX packages and the \texttt{Source/} modules (\texttt{CCZ4}, \texttt{Matter}, \texttt{GRTeclynCore}, \texttt{Tagging}, \texttt{AMRInterpolator}, etc.) that are compiled into the executable.
    \item \textbf{\texttt{Make.package}.} Registers local \texttt{.cpp} sources and \texttt{.hpp} headers into \texttt{GRTECLYN\_CEXE\_*} variables consumed by the build system.
    \item \textbf{\texttt{Main\_*.cpp}.} Standard \texttt{GRAMR} time loop: load parameters, construct a \texttt{DefaultLevelFactory<Level>}, initialize AMR, and advance until \texttt{stop\_time}.
    \item \textbf{\texttt{*Level.hpp/.cpp}.} Extends \texttt{GRAMRLevel}. Wires initial data, RHS evolution, AMR tagging, and post-timestep diagnostics through virtual hooks such as \texttt{initData}, \texttt{specificEvalRHS}, and \texttt{specificPostTimeStep}.
    \item \textbf{\texttt{SimulationParameters.hpp}.} Problem-specific parameter class. Must live in the example directory because \texttt{GRTeclynCore} includes it at compile time.
    \item \textbf{\texttt{StateVariables.hpp}.} Defines \texttt{NUM\_VARS}, variable names, and boundary-condition parities. Matter examples append \texttt{c\_phi} and \texttt{c\_Pi} after the CCZ4 indices.
    \item \textbf{\texttt{*InitialData.hpp}.} GPU initial-data functor. Writes all CCZ4 and matter components at $t=0$.
    \item \textbf{Potential / matter header.} Example-local potential (e.g.\ \texttt{PhantomDecayPotential.hpp}) or a shared \texttt{DefaultPotential} from \texttt{Source/Matter/}.
    \item \textbf{\texttt{params.txt}.} Runtime inputs: grid, BCs, evolution, output paths, and physics parameters read by \texttt{GRParmParse}.
\end{itemize}

\textbf{Compile-time header injection.} \texttt{SimulationParameters.hpp} and \texttt{StateVariables.hpp} are \emph{not} in \texttt{Source/}. The example directory is first on \texttt{INCLUDE\_LOCATIONS}, so each example supplies its own version. This is why copying an example is the standard way to start a new run.

\textbf{Template type consistency.} If the example uses matter, the exact matter type (e.g.\ \texttt{ScalarField<DefaultPotential>} or \texttt{ExoticScalarField<PhantomDecayPotential>}) must match in \texttt{variableSetUp}, \texttt{specificEvalRHS}, constraint diagnostics, and Weyl4 setup. A mismatch often compiles in one place and fails elsewhere.

\textbf{Parameter loading chain.} \texttt{AMReXParameters} reads grid, BC, and output settings. \texttt{SimulationParametersBase} reads CCZ4 gauge, floors, and optional extraction geometry. The derived \texttt{SimulationParameters} class reads problem-specific keys and validates them.

\textbf{Wrapper vs C++ boundary.} The \texttt{grteclyn-wrapper} package automates episode directories, \texttt{params.txt} templating, \texttt{check\_params=1}, execution, log capture, diagnostic parsing, scoring, and atlas batching. C++ work is still required for new initial-data recipes, new matter models, new runtime diagnostics, and new extraction hooks.

\textbf{Recommended starting point.} \texttt{Examples/RadialRecipe/} is a minimal CCZ4+matter template that reads Gaussian radial recipe coefficients from \texttt{params.txt}. It uses \texttt{ScalarField<DefaultPotential>} and the same collapse diagnostics as \texttt{SupportedWormholeCollapse}, but replaces the hand-written Ellis--Bronnikov ansatz with optimizer-friendly basis coefficients (\texttt{recipe\_chi\_coeff\_0}, \texttt{recipe\_phi\_coeff\_0}, etc.).

\subsection{C++ Interpolation and Extraction Options}

Several paths exist for turning evolved grid data into machine-readable observables:

\begin{enumerate}
    \item \textbf{Volume reductions (current baseline).} Implemented in \texttt{specificPostTimeStep} via \texttt{ParallelFor}, \texttt{ReduceOps}, and \texttt{SmallDataIO}. This is what \texttt{SupportedWormholeCollapse} and \texttt{RadialRecipe} use today for \texttt{collapse\_diagnostics.dat} and \texttt{constraint\_norms.dat}. It works now and is ideal for atlas scoring, but it only produces aggregate quantities (min/max/norms, horizon proxy, scalar extrema).
    \item \textbf{Surface extraction via \texttt{SphericalExtraction}.} Samples variables on spherical shells and can project onto spin-weighted spherical harmonics (used by \texttt{WeylExtraction} for $\Psi_4$ modes). The API lives in \texttt{Source/AMRInterpolator/}, but the AMReX-ported \texttt{AMRInterpolator::refresh()} path is incomplete, so runtime GW extraction is not yet wired up in current examples.
    \item \textbf{\texttt{CylindricalExtraction}.} Same surface-extraction framework on cylindrical shells. Useful for axial throat profiles and jet-like structures.
    \item \textbf{\texttt{ParticleInterpolator}.} AMReX-native interpolation at arbitrary points. A lighter-weight option for probe lines or tracer particles without the full surface-extraction stack.
    \item \textbf{Plotfile post-processing.} AMReX tools (\texttt{fextract}, \texttt{fcompare}) or Python (\texttt{yt}, AMReX bindings) can read plotfiles after the run. No C++ changes are needed, but this conflicts with long-run data retention unless plotfiles are kept for promoted candidates only.
\end{enumerate}

\noindent\textbf{Recommendation.} Use volume reductions for the immediate atlas pipeline; restore spherical extraction as a milestone for $\Psi_4$ regression against the previous wormhole paper; reserve plotfile post-processing for one-off analysis of promoted candidates.

\section{Baseline Workflow}

One search episode is:

\begin{enumerate}
    \item Python proposes parameters $\theta$.
    \item Python performs fast sanity checks.
    \item Python writes a \texttt{GRTeclyn} input file.
    \item \texttt{GRTeclyn} builds the initial data.
    \item \texttt{GRTeclyn} runs a short low-resolution evolution.
    \item \texttt{GRTeclyn} writes diagnostics.
    \item Python reads the diagnostics and computes a score.
    \item The optimizer proposes the next candidate.
\end{enumerate}

The first version should be optimizer-agnostic:
\begin{verbatim}
theta = optimizer.ask()
if not sanity_check(theta):
    return bad_score

write_grteclyn_input(theta)
run_grteclyn()
diagnostics = read_outputs()
score = compute_score(diagnostics)
optimizer.tell(theta, score)
\end{verbatim}

\section{Scoring}

The score should be simple at first. Good behavior includes survival, bounded constraints, nontrivial geometry, an open throat, early signal arrival, low tidal forces, and reduced exotic matter compared with known examples.

Bad behavior includes constraint blow-up, trivial flat spacetime, immediate horizon formation, gauge crash, unresolved AMR gradients, huge negative energy, or disappearance of the effect at higher resolution.

The important rule is that a candidate is not rewarded for coordinate tricks. A large shift vector alone is not enough. The target property must survive operational checks.

\section{Avoiding Fake FTL}

The main scientific risk is false positives. A candidate can look superluminal because of coordinates, gauge motion, boundary effects, or constraint violation. Therefore promoted candidates must pass operational tests:

\begin{enumerate}
    \item send a null ray or scalar pulse from an emitter to a detector;
    \item compare arrival time with a flat-spacetime reference;
    \item verify that the signal path does not cross a trapped surface;
    \item verify that constraints remain bounded before arrival;
    \item rerun at higher resolution;
    \item rerun with slightly different gauge parameters;
    \item perturb the initial data and check that the effect survives.
\end{enumerate}

This converts the project from a search for coordinate artifacts into a search for reproducible dynamical behavior.

\section{Benchmarks}

Before searching for new metrics, the pipeline must reproduce known cases.

\begin{table}[h]
\caption{Initial benchmark suite for the search pipeline.}
\label{tab:benchmarks}
\begin{ruledtabular}
\begin{tabular}{lll}
Case & Expected result & Purpose \\
\hline
Minkowski & Stable, no shortcut & Null test \\
Ellis--Bronnikov & Unstable throat & Wormhole detector \\
Perturbed EB & Collapse + GW signal & Previous-paper regression \\
Alcubierre-like & High exoticity & Fake-FTL test \\
Teo-like & Quadrupole dynamics & Rotation benchmark \\
Random data & Mostly failures & Failure database \\
\end{tabular}
\end{ruledtabular}
\end{table}

\section{Step-by-Step Implementation Plan}

\subsection{Step 1: Reproduce the Previous Runs}

Create one script that runs the existing Ellis--Bronnikov setup and the perturbed collapse setup. Success means the throat radius, constraint behavior, horizon proxy, and $\Psi_4$ extraction match the previous study. This becomes the regression test.

\subsection{Step 2: Create Episode Logging}

Every run should create one episode folder containing input parameters, grid settings, git commit hash, runtime, termination reason, diagnostics, and final score. Failed runs are kept because they teach the optimizer what not to do.

\subsection{Step 3: Add Radial Recipe Initial Data}

Implement the first new C++ initial-data class. It reads radial coefficients and generates smooth profiles for $\chi$, $\alpha$, $\beta^i$, $K$, $\phi$, and $\Pi$. At this stage there is no optimization yet; hand-picked coefficients are enough.

\noindent\textbf{Implemented status.}
\texttt{Examples/RadialRecipe/} provides a buildable CCZ4+matter template with Gaussian radial basis functions. Coefficients are read from indexed \texttt{params.txt} keys such as \texttt{recipe\_chi\_coeff\_0}. The wrapper supports this example via \texttt{--example RadialRecipe}.

\subsection{Step 4: Add the Python Wrapper}

Write a Python script that creates input files, launches \texttt{GRTeclyn}, waits for completion, reads diagnostics, and computes a simple score. This creates the first complete episode loop.

\noindent\textbf{Implemented status.}
A standalone wrapper repository now exists in \texttt{grteclyn-wrapper}. It creates isolated episode directories, templates \texttt{params\_2gpu.txt}, launches \texttt{SupportedWormholeCollapse}, runs the \texttt{check\_params=1} preflight, writes \texttt{metadata.json}, streams \texttt{run.log}, parses \texttt{collapse\_diagnostics.dat} and \texttt{constraint\_norms.dat}, and writes \texttt{score.json}. It also includes shell scripts for medium, big, full, and atlas runs on a selected GPU.

\subsection{Step 5: Run Random Search}

Sample many radial candidates at low resolution. The goal is not discovery yet. The goal is to build a failure atlas: which parameters crash, which become flat, which form horizons, which create throats, and which diagnostics are useful.

\noindent\textbf{Implemented status.}
The wrapper now includes an \texttt{atlas} mode. It samples existing wormhole parameters, runs one isolated episode per sample, classifies the outcome, and writes \texttt{atlas.csv}, \texttt{atlas.jsonl}, and \texttt{summary.json}. Initial labels include \texttt{completed}, \texttt{missing\_diagnostics}, \texttt{constraint\_blowup}, \texttt{lapse\_collapse}, \texttt{horizon\_formed}, \texttt{trivial\_geometry}, and \texttt{solver\_failed}. A single-GPU smoke atlas completed successfully, and a larger \texttt{N\_full=128}, \texttt{max\_level=4} atlas test confirmed that the data pipeline works even when sampled candidates fail early with horizon formation and NaNs.

\subsection{Step 5b: Production Data Retention}

Final production and atlas runs must not keep all raw AMReX output indefinitely. Each episode should extract the useful machine-readable information first, then delete heavy files that are not needed for scoring or later analysis. The default retained set should be:

\begin{itemize}
    \item \texttt{params.txt}, \texttt{metadata.json}, \texttt{run.log}, and \texttt{score.json};
    \item \texttt{data/collapse\_diagnostics.dat} and \texttt{data/constraint\_norms.dat};
    \item selected \texttt{small\_data/*.dat} products such as \texttt{psi4\_mode\_l2m0.dat} and \texttt{areal\_radius.dat};
    \item the batch-level \texttt{atlas.csv}, \texttt{atlas.jsonl}, and \texttt{summary.json}.
\end{itemize}

After extraction, the wrapper should optionally delete heavy episode-local outputs such as \texttt{SupportedWormholePlt*}, \texttt{SupportedWormholeChk*}, \texttt{hdf5/}, \texttt{pout/}, and generated frame folders. This keeps long atlas runs from overwhelming disk storage while preserving the diagnostics needed for scoring, optimizer training, and failure analysis.

\subsection{Step 6: Add a Simple Optimizer}

Use CMA-ES or Bayesian optimization before reinforcement learning. These methods are simpler, cheaper, and easier to debug. RL should come only after the rewards and episode loop are stable.

\subsection{Step 7: Add Probe Tests}

Add null rays, scalar pulses, or timelike tracer particles. This is the transition from metric-shape scoring to real transport and communication tests.

\subsection{Step 8: Add Angular Modes}

Extend the recipe from radial profiles to angular modes. This opens the search to asymmetric wormholes, portal-like mouths, rotating structures, and gravitational-wave-producing candidates.

\subsection{Step 9: Promote Candidates}

High-scoring candidates must pass a resolution ladder,
\begin{equation}
    N=64 \rightarrow 128 \rightarrow 256,
\end{equation}
and also survive longer runtime, small perturbations, changed gauge parameters, and stricter constraint thresholds.

\subsection{Step 10: Add Neural Decoders and RL}

Only after the previous steps work should neural decoders and reinforcement learning be added. By then the project will already have a working episode loop, a diagnostics database, calibrated rewards, known benchmarks, and failure statistics.

\section{Near-Term Milestones}

The near-term milestones are:

\begin{enumerate}
    \item one-command reproduction of the previous wormhole runs;
    \item machine-readable diagnostics from \texttt{GRTeclyn};
    \item first radial recipe initial-data class in C++; \textbf{implemented as \texttt{Examples/RadialRecipe/}}
    \item Python wrapper for one episode; \textbf{implemented}
    \item random search over current wormhole parameters; \textbf{implemented as the first atlas mode}
    \item first failure-mode database; \textbf{implemented as \texttt{atlas.csv}/\texttt{atlas.jsonl}}
    \item production cleanup mode that extracts useful data and deletes heavy raw outputs;
    \item first optimizer-driven search;
    \item first null-ray or scalar-pulse validation.
\end{enumerate}

\section{Expected First Result}

A good first result does not need to be a stable FTL metric. A good first result is a working search loop connected to \texttt{GRTeclyn}, reproduction of the previous wormhole paper, automatic classification of known metrics, a database showing how random exotic initial data fail, and one or two nontrivial candidates promoted to higher resolution.

This is already publishable as a baseline. The larger goal is to expand the search space until the system either finds new viable candidates or gives strong evidence that broad classes of exotic metrics fail dynamically.

\section{Conclusion}

The project should be presented as an automated numerical-relativity laboratory for searching and testing exotic spacetime metrics. The power is not just AI. The power is AI plus \texttt{GRTeclyn}: the optimizer proposes, but full 3D Einstein evolution decides.

\end{document}
\documentclass[12pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{geometry}
\usepackage{hyperref}
\usepackage{enumitem}
\usepackage{booktabs}
\geometry{a4paper, margin=2cm}

\title{\textbf{Neural Spacetime Search with GRTeclyn: Practical Baseline Plan}}
\author{Nikita M. Shirokov}
\date{\today}

\begin{document}
\maketitle

\section*{Core Idea}

The previous wormhole paper showed that \texttt{GRTeclyn} can evolve exotic spacetime geometries in full 3D numerical relativity. It evolved the Ellis--Bronnikov wormhole, saw the collapse/expansion instability, extracted gravitational waves, and identified the phantom bounce.

The next step is to turn this into a search engine.

Instead of manually inventing one metric, checking it, and moving to the next one, we let an optimizer generate many candidate initial data sets. \texttt{GRTeclyn} evolves each candidate. The pipeline then scores what happened.

The goal is to search for metrics with useful interstellar properties:

\begin{itemize}[leftmargin=*]
    \item effective faster-than-light transport;
    \item faster-than-light communication through geometry;
    \item portal-like shortcuts;
    \item wormholes or throat-like bridges;
    \item new exotic compact-object dynamics and gravitational-wave signatures.
\end{itemize}

Important: the claim is not ``we will build a working warp drive.'' The claim is more careful and stronger scientifically:

\begin{quote}
We will build a closed-loop numerical-relativity search system that tests candidate spacetime metrics by their nonlinear time evolution, not only by static formulas.
\end{quote}

Even if most candidates fail, the result is useful. We learn how they fail: collapse, dispersion, horizon formation, constraint violation, gauge crash, or excessive exotic matter.

\section*{What The Agent Searches}

\texttt{GRTeclyn} evolves the spacetime in CCZ4 variables. So the search should work directly with those variables.

At the first stage, the optimizer generates simple recipes for:

\begin{itemize}[leftmargin=*]
    \item $\chi$: the conformal factor; useful for throats, stretching, and wormhole-like regions;
    \item $\alpha$: lapse;
    \item $\beta^i$: shift; useful for warp-like motion;
    \item $K$: expansion or contraction of space;
    \item $\tilde A_{ij}$: trace-free extrinsic curvature;
    \item $\phi$ and $\Pi$: scalar-field matter sector.
\end{itemize}

The optimizer should not generate every grid point directly. That is too many parameters. It should generate a small parameter vector, and a decoder turns it into smooth grid fields.

\section*{Simple Search Levels}

\subsection*{Level 1: Radial Recipes}

Start with spherical or nearly spherical fields:

\[
q(r) = q_\infty + \sum_n a_n B_n(r).
\]

This is the easiest first step. It can reproduce known wormhole-like throats and simple bubble profiles. It is also cheap enough for many test runs.

\subsection*{Level 2: Angular Recipes}

Add angular structure:

\[
q(r,\theta,\varphi) = q_\infty + \sum_{n\ell m} a_{n\ell m}B_n(r)Y_{\ell m}(\theta,\varphi).
\]

This allows asymmetric mouths, quadrupoles, rotating-like structures, and gravitational-wave-producing candidates.

\subsection*{Level 3: Neural Decoder}

Later, use a small neural network:

\[
(x,y,z) \rightarrow \{\chi,\alpha,\beta^i,K,\tilde A_{ij},\phi,\Pi\}.
\]

This should not be the first step. It is too flexible and too hard to debug. Use it only after the simpler pipeline works.

\section*{How We Connect This To GRTeclyn}

The clean split is:

\begin{itemize}[leftmargin=*]
    \item \textbf{C++ / GRTeclyn:} evolves the spacetime, computes diagnostics, writes outputs.
    \item \textbf{Python:} chooses candidate parameters, writes input files, launches runs, reads outputs, computes reward, updates optimizer.
\end{itemize}

At first, Python can run \texttt{GRTeclyn} as a subprocess. We do not need a complicated API immediately.

\section*{Minimal GRTeclyn Changes}

We need small, practical C++ additions.

\begin{enumerate}[leftmargin=*]
    \item \textbf{Recipe initial-data class.}
    Add a new initial-data module in \texttt{GRTeclyn}. It reads coefficients from the input file and fills $\chi$, $\alpha$, $\beta^i$, $K$, $\tilde A_{ij}$, $\phi$, and $\Pi$ on the AMReX grid.

    \item \textbf{Standard diagnostics output.}
    Write simple machine-readable files for:
    constraints, minimum lapse, minimum $\chi$, maximum $|K|$, throat radius, horizon proxy, scalar-field extrema, and run termination reason.

    \item \textbf{Early-stop conditions.}
    Stop a run if constraints explode, lapse collapses everywhere, curvature becomes unresolved, or wall time is exceeded. This saves GPU time.

    \item \textbf{Probe module.}
    Add simple test particles, null rays, or scalar pulses. These check whether communication or transport is real, not just a coordinate artifact.

    \item \textbf{Optional field loader.}
    Later, allow Python to write full grid fields to disk and let \texttt{GRTeclyn} load them. This is useful for neural decoders.
\end{enumerate}

\section*{Anatomy of a GRTeclyn Example}

Every runnable \texttt{GRTeclyn} problem is an \texttt{Examples/<Name>/} directory with roughly ten tightly coupled files. The Python wrapper can template \texttt{params.txt}, launch the binary, and parse diagnostics, but it cannot invent new physics without these C++ pieces.

\begin{itemize}[leftmargin=*]
    \item \textbf{\texttt{GNUmakefile}.} Build entry point. Selects AMReX packages and the \texttt{Source/} modules (\texttt{CCZ4}, \texttt{Matter}, \texttt{GRTeclynCore}, \texttt{Tagging}, \texttt{AMRInterpolator}, etc.).
    \item \textbf{\texttt{Make.package}.} Registers local \texttt{.cpp} sources and \texttt{.hpp} headers for the build.
    \item \textbf{\texttt{Main\_*.cpp}.} Standard \texttt{GRAMR} time loop.
    \item \textbf{\texttt{*Level.hpp/.cpp}.} Extends \texttt{GRAMRLevel} and wires initial data, RHS evolution, tagging, and diagnostics.
    \item \textbf{\texttt{SimulationParameters.hpp}.} Problem-specific parameters. Must live in the example directory because core headers include it at compile time.
    \item \textbf{\texttt{StateVariables.hpp}.} Defines \texttt{NUM\_VARS}, names, and BC parities.
    \item \textbf{\texttt{*InitialData.hpp}.} GPU functor that fills CCZ4 and matter fields at $t=0$.
    \item \textbf{Potential / matter header.} Example-local or shared from \texttt{Source/Matter/}.
    \item \textbf{\texttt{params.txt}.} Runtime inputs read by \texttt{GRParmParse}.
\end{itemize}

\textbf{Compile-time header injection.} \texttt{SimulationParameters.hpp} and \texttt{StateVariables.hpp} are resolved from the example directory, not from \texttt{Source/}. Copying an example is the standard way to start a new run.

\textbf{Template type consistency.} The matter type must match everywhere: RHS, constraints, Weyl4, and derived variables.

\textbf{Parameter loading chain.} \texttt{AMReXParameters} $\rightarrow$ \texttt{SimulationParametersBase} $\rightarrow$ derived \texttt{SimulationParameters}.

\textbf{Wrapper vs C++ boundary.} \texttt{grteclyn-wrapper} automates params templating, execution, parsing, scoring, and atlas batches. C++ is still required for new initial data, matter models, diagnostics, and extraction hooks.

\textbf{Recommended starting point.} \texttt{Examples/RadialRecipe/} reads Gaussian radial recipe coefficients from \texttt{params.txt} and is the template for optimizer-driven initial-data search.

\section*{C++ Interpolation and Extraction Options}

\begin{enumerate}[leftmargin=*]
    \item \textbf{Volume reductions (current baseline).} \texttt{ParallelFor} + \texttt{ReduceOps} + \texttt{SmallDataIO} in \texttt{specificPostTimeStep}. Works today; good for atlas scoring; limited to aggregate diagnostics.
    \item \textbf{\texttt{SphericalExtraction}.} Shell sampling and harmonic mode integration (GW / $\Psi_4$). API exists in \texttt{Source/AMRInterpolator/}, but runtime extraction is not fully wired in the current AMReX port.
    \item \textbf{\texttt{CylindricalExtraction}.} Same framework on cylindrical shells; useful for axial profiles.
    \item \textbf{\texttt{ParticleInterpolator}.} Point-wise GPU interpolation without the full surface stack.
    \item \textbf{Plotfile post-processing.} \texttt{fextract}, \texttt{yt}, etc. No C++ changes, but heavy files must be retained.
\end{enumerate}

\textbf{Recommendation.} Use volume reductions for atlas runs now; restore spherical extraction for GW regression; keep plotfile analysis for promoted candidates.

\section*{First Baseline Workflow}

One search episode should look like this:

\begin{enumerate}[leftmargin=*]
    \item Python picks parameters $\theta$.
    \item Python checks that the parameters are not obviously bad.
    \item Python writes a \texttt{GRTeclyn} input file.
    \item \texttt{GRTeclyn} builds the initial data from $\theta$.
    \item \texttt{GRTeclyn} runs a short low-resolution evolution.
    \item \texttt{GRTeclyn} writes diagnostics.
    \item Python reads diagnostics.
    \item Python computes a score.
    \item The optimizer proposes the next $\theta$.
\end{enumerate}

In pseudocode:

\begin{verbatim}
theta = optimizer.ask()
if not sanity_check(theta):
    return bad_score

write_grteclyn_input(theta)
run_grteclyn()
diagnostics = read_outputs()
score = compute_score(diagnostics)
optimizer.tell(theta, score)
\end{verbatim}

\section*{What We Score}

The score should be simple at first.

\subsection*{Good Things}

\begin{itemize}[leftmargin=*]
    \item The run survives until the requested final time.
    \item Constraints stay bounded.
    \item The geometry does something nontrivial.
    \item A throat stays open for some time.
    \item A signal arrives earlier than in flat spacetime.
    \item A passenger region has low tidal forces.
    \item Exotic matter is smaller than in known examples.
\end{itemize}

\subsection*{Bad Things}

\begin{itemize}[leftmargin=*]
    \item Hamiltonian or momentum constraints explode.
    \item The result is just flat spacetime.
    \item The candidate immediately collapses into a horizon.
    \item The effect disappears at higher resolution.
    \item The effect depends only on coordinates or gauge.
    \item The required negative energy is huge.
    \item AMR cannot resolve the gradients.
\end{itemize}

\section*{Avoiding Fake FTL}

We must be careful. A big shift vector or strange coordinates can look like faster-than-light motion, but may not be physical.

So a candidate only counts if it passes operational tests:

\begin{itemize}[leftmargin=*]
    \item send a null ray or scalar pulse from A to B;
    \item compare arrival time with flat spacetime;
    \item check the signal path does not cross a trapped surface;
    \item check constraints stay bounded before the signal arrives;
    \item rerun at higher resolution;
    \item rerun with slightly different gauge parameters;
    \item add small perturbations and check the effect survives.
\end{itemize}

This is the key scientific point. The system does not reward coordinate tricks. It rewards behavior that survives evolution and validation.

\section*{Benchmarks}

Before searching for anything new, we test the pipeline on known cases.

\begin{center}
\begin{tabular}{@{}lll@{}}
\toprule
\textbf{Case} & \textbf{Expected result} & \textbf{Purpose} \\
\midrule
Minkowski & Stable, no shortcut & Null test \\
Ellis--Bronnikov & Unstable throat & Wormhole detector test \\
Perturbed EB collapse & Collapse + GW signal & Reproduce previous paper \\
Alcubierre-like shift & Kinematic warp, high exoticity & Fake-FTL test \\
Teo-like rotation & Quadrupole dynamics & Rotating benchmark \\
Random smooth data & Mostly failures & Build failure database \\
\bottomrule
\end{tabular}
\end{center}

\section*{Step-by-Step Plan}

\subsection*{Step 1: Reproduce Previous Results}

Make one script that runs the existing Ellis--Bronnikov setup and the perturbed collapse setup.

Success means:

\begin{itemize}[leftmargin=*]
    \item throat radius matches the old run;
    \item constraint behavior matches the old run;
    \item horizon proxy appears in the collapse run;
    \item $\Psi_4$ extraction still works.
\end{itemize}

This becomes the regression test.

\subsection*{Step 2: Create Episode Logging}

Every run should create one episode folder:

\begin{itemize}[leftmargin=*]
    \item input parameters;
    \item grid settings;
    \item git commit hash;
    \item runtime;
    \item termination reason;
    \item diagnostics;
    \item final score.
\end{itemize}

This gives us a dataset. Failed runs are useful too.

\subsection*{Step 3: Add Radial Recipe Initial Data}

Add the first new C++ initial-data class.

It should read simple radial coefficients and generate smooth profiles for $\chi$, $\alpha$, $\beta^i$, $K$, $\phi$, and $\Pi$.

Do not optimize yet. Just run hand-picked coefficients and check that the system works.

\noindent\textbf{Implemented status.}
\texttt{Examples/RadialRecipe/} provides a buildable CCZ4+matter template with Gaussian radial basis functions. Coefficients are read from indexed \texttt{params.txt} keys such as \texttt{recipe\_chi\_coeff\_0}. The wrapper supports this example via \texttt{--example RadialRecipe}.

\subsection*{Step 4: Add Python Wrapper}

Create a Python script that:

\begin{itemize}[leftmargin=*]
    \item writes the input file;
    \item launches \texttt{GRTeclyn};
    \item watches for completion;
    \item reads diagnostics;
    \item computes a simple score.
\end{itemize}

At this stage, random sampling is enough.

\noindent\textbf{Implemented status.}
A standalone wrapper repository now exists in \texttt{grteclyn-wrapper}. It creates isolated episode directories, templates \texttt{params\_2gpu.txt}, launches \texttt{SupportedWormholeCollapse}, runs the \texttt{check\_params=1} preflight, writes \texttt{metadata.json}, streams \texttt{run.log}, parses \texttt{collapse\_diagnostics.dat} and \texttt{constraint\_norms.dat}, and writes \texttt{score.json}. It also includes shell scripts for medium, big, full, and atlas runs on a selected GPU.

\subsection*{Step 5: Run Random Search}

Sample many radial candidates at low resolution.

Goal: not discovery yet. Goal: build a map of failure modes.

Questions to answer:

\begin{itemize}[leftmargin=*]
    \item Which parameters always crash?
    \item Which ones make flat spacetime?
    \item Which ones make throats?
    \item Which ones form horizons?
    \item Which diagnostics are most useful?
\end{itemize}

\noindent\textbf{Implemented status.}
The wrapper now includes an \texttt{atlas} mode. It samples existing wormhole parameters, runs one isolated episode per sample, classifies the outcome, and writes \texttt{atlas.csv}, \texttt{atlas.jsonl}, and \texttt{summary.json}. Initial labels include \texttt{completed}, \texttt{missing\_diagnostics}, \texttt{constraint\_blowup}, \texttt{lapse\_collapse}, \texttt{horizon\_formed}, \texttt{trivial\_geometry}, and \texttt{solver\_failed}. A single-GPU smoke atlas completed successfully, and a larger \texttt{N\_full=128}, \texttt{max\_level=4} atlas test confirmed that the data pipeline works even when sampled candidates fail early with horizon formation and NaNs.

\subsection*{Step 5b: Production Data Retention}

Final production and atlas runs must not keep all raw AMReX output indefinitely. Each episode should extract the useful machine-readable information first, then delete heavy files that are not needed for scoring or later analysis. The default retained set should be:

\begin{itemize}[leftmargin=*]
    \item \texttt{params.txt}, \texttt{metadata.json}, \texttt{run.log}, and \texttt{score.json};
    \item \texttt{data/collapse\_diagnostics.dat} and \texttt{data/constraint\_norms.dat};
    \item selected \texttt{small\_data/*.dat} products such as \texttt{psi4\_mode\_l2m0.dat} and \texttt{areal\_radius.dat};
    \item the batch-level \texttt{atlas.csv}, \texttt{atlas.jsonl}, and \texttt{summary.json}.
\end{itemize}

After extraction, the wrapper should optionally delete heavy episode-local outputs such as \texttt{SupportedWormholePlt*}, \texttt{SupportedWormholeChk*}, \texttt{hdf5/}, \texttt{pout/}, and generated frame folders. This keeps long atlas runs from overwhelming disk storage while preserving the diagnostics needed for scoring, optimizer training, and failure analysis.

\subsection*{Step 6: Add A Simple Optimizer}

Use CMA-ES or Bayesian optimization before RL.

Reason: they are simpler, cheaper, and easier to debug.

Only after this works should we add PPO or another RL method.

\subsection*{Step 7: Add Probe Tests}

Add null rays, scalar pulses, or timelike tracer particles.

This is where the project becomes about real transport/communication tests instead of only metric shape.

Success means we can measure:

\begin{itemize}[leftmargin=*]
    \item signal travel time;
    \item passenger tidal forces;
    \item whether the path crosses a horizon;
    \item whether the result beats the flat-spacetime reference.
\end{itemize}

\subsection*{Step 8: Add Angular Modes}

Extend the recipe from radial profiles to angular modes.

This lets the search find:

\begin{itemize}[leftmargin=*]
    \item asymmetric wormholes;
    \item portal-like two-mouth structures;
    \item rotating or twisting geometries;
    \item quadrupolar gravitational-wave sources.
\end{itemize}

\subsection*{Step 9: Candidate Promotion}

Promising candidates must pass a ladder:

\[
N=64 \rightarrow 128 \rightarrow 256.
\]

They must also survive:

\begin{itemize}[leftmargin=*]
    \item longer runtime;
    \item small perturbations;
    \item changed gauge parameters;
    \item stricter constraint thresholds.
\end{itemize}

If a candidate fails here, it is probably a numerical artifact.

\subsection*{Step 10: Neural Decoder And RL}

Only after the previous steps work, add neural decoders and RL.

At that point we will already have:

\begin{itemize}[leftmargin=*]
    \item a working episode loop;
    \item a diagnostics database;
    \item calibrated rewards;
    \item known benchmarks;
    \item failure statistics.
\end{itemize}

Then the neural model can search a much larger space, and \texttt{GRTeclyn} remains the final judge.

\section*{Near-Term Milestones}

\begin{enumerate}[leftmargin=*]
    \item One-command reproduction of previous wormhole runs.
    \item Machine-readable diagnostics from \texttt{GRTeclyn}.
    \item First radial recipe initial-data class in C++. \textbf{Implemented as \texttt{Examples/RadialRecipe/}.}
    \item Python wrapper for one episode. \textbf{Implemented.}
    \item Random search over current wormhole parameters. \textbf{Implemented as the first atlas mode.}
    \item First failure-mode database. \textbf{Implemented as \texttt{atlas.csv}/\texttt{atlas.jsonl}.}
    \item Production cleanup mode that extracts useful data and deletes heavy raw outputs.
    \item First optimizer-driven search.
    \item First null-ray or scalar-pulse validation.
\end{enumerate}

\section*{What A Good First Result Looks Like}

A good first result does not need to be a stable FTL metric.

A good first result is:

\begin{itemize}[leftmargin=*]
    \item a working search loop connected to \texttt{GRTeclyn};
    \item reproduction of the previous wormhole paper;
    \item automatic classification of known metrics;
    \item a database showing how random exotic initial data fail;
    \item one or two nontrivial candidates promoted to higher resolution.
\end{itemize}

That is already a publishable baseline. The larger goal is to keep expanding the search space until the system either finds new viable candidates or gives strong evidence that broad classes of candidates fail dynamically.

\section*{Main Message}

The project should be presented as:

\begin{quote}
An automated numerical-relativity laboratory for searching and testing exotic spacetime metrics.
\end{quote}

The power is not just AI. The power is AI plus \texttt{GRTeclyn}: the optimizer proposes, but full 3D Einstein evolution decides.

\end{document}
