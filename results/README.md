# Publishable simulation results

Light, GitHub-friendly extracts of campaign outputs used by papers.
Full evolutions (frames, plotfiles, gridinit, metric stacks) stay in the
gitignored `/runs` tree on the machine that produced them.

| Paper pack | Article | Regenerator |
|---|---|---|
| [`matter-first-automated-discovery-of-transient-spacetime-shortcuts/`](matter-first-automated-discovery-of-transient-spacetime-shortcuts/) | [`research/neuralspacetime/article/research.tex`](../research/neuralspacetime/article/research.tex) | [`research/neuralspacetime/pack_publishable_results.sh`](../research/neuralspacetime/pack_publishable_results.sh) |
| [`bondi-dipole-runaway/`](bondi-dipole-runaway/) | in preparation — source material in the pack itself | [`research/bondi_dipole/pack_results.sh`](../research/bondi_dipole/pack_results.sh) |

### `bondi-dipole-runaway/`

A positive-active-mass and a negative-active-mass soliton, released at rest in
full 3+1 NR with constraint-solved matter, self-accelerate together; both
same-sector controls are null. Holds the per-cell time series, dressed-star
tables, solve/evolution parameters, code patches, curated frames and movies, and
the derived tables — plus the physics findings, the debugging trail, the matter
model, and the launch reference as standalone documents.

## Campaign packs

Deepest evolving-geodesic path saving per campaign, newest first. Depths are
comparable across rows; scores are not (they are objective-mode specific).

| Pack | Campaign | Search | Evals | Best depth |
|---|---|---|---|---|
| [`qball-trajectory-cmaes-refinement/`](qball-trajectory-cmaes-refinement/) | `qball_traj_fgeo_depth_cmaes_v1` | CMA-ES | 200 | **48.38 %** (47.97 % conservative) |
| [`qball-trajectory-geodesic-depth-search/`](qball-trajectory-geodesic-depth-search/) | `qball_traj_fgeo_depth_v1` | MAP-Elites | 200 | 45.90 % |
| [`qball-trajectory-evolving-geodesic-shortcut-search/`](qball-trajectory-evolving-geodesic-shortcut-search/) | `qball_traj_fgeo_v1` | MAP-Elites | — | — |
