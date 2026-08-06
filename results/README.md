# Publishable simulation results

Light, GitHub-friendly extracts of campaign outputs used by papers.
Full evolutions (frames, plotfiles, gridinit, metric stacks) stay in the
gitignored `/runs` tree on the machine that produced them.

| Paper pack | Article | Regenerator |
|---|---|---|
| [`matter-first-automated-discovery-of-transient-spacetime-shortcuts/`](matter-first-automated-discovery-of-transient-spacetime-shortcuts/) | [`research/neuralspacetime/article/research.tex`](../research/neuralspacetime/article/research.tex) | [`research/neuralspacetime/pack_publishable_results.sh`](../research/neuralspacetime/pack_publishable_results.sh) |

## Campaign packs

Deepest evolving-geodesic path saving per campaign, newest first. Depths are
comparable across rows; scores are not (they are objective-mode specific).

| Pack | Campaign | Search | Evals | Best depth |
|---|---|---|---|---|
| [`qball-trajectory-cmaes-refinement/`](qball-trajectory-cmaes-refinement/) | `qball_traj_fgeo_depth_cmaes_v1` | CMA-ES | 200 | **48.38 %** (47.97 % conservative) |
| [`qball-trajectory-geodesic-depth-search/`](qball-trajectory-geodesic-depth-search/) | `qball_traj_fgeo_depth_v1` | MAP-Elites | 200 | 45.90 % |
| [`qball-trajectory-evolving-geodesic-shortcut-search/`](qball-trajectory-evolving-geodesic-shortcut-search/) | `qball_traj_fgeo_v1` | MAP-Elites | — | — |
