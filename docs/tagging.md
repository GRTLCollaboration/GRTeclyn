# Tagging criteria

During regridding, each tagging criterion marks cells that should be covered by the next refinement level. If several taggers are used, a cell is refined when any one of them tags it. AMReX then groups the tagged cells into grid boxes, subject to its blocking and proper-nesting requirements.

## `ChiTagger`

`ChiTagger` refines regions where the conformal factor \(\chi\) is varying rapidly. Its criterion is based on the magnitude of the second spatial derivatives of \(\chi\), scaled by the grid spacing.

This tagger is used by the `BinaryBH` example.

- `tagging.threshold` sets the same tagging threshold on every level. Smaller values generally tag more cells.
- `tagging.thresholds` optionally supplies a separate threshold for each AMR level and takes precedence over `tagging.threshold`.

## `PunctureTagger`

`PunctureTagger` maintains refinement around each tracked puncture. The tagged radius is proportional to the puncture mass. Around the finest level it is approximately

\[
r_{\mathrm{tag}} =
\texttt{finest\_level\_factor}\,M.
\]

On the next coarser tagging level, this radius is multiplied by `level_separation`. This helps keep the finest refinement boundaries away from the puncture.

This tagger is used by the `BinaryBH` example when puncture tracking is enabled.

- `puncture_tagging.finest_level_factor` controls the size of the finest tagged region relative to the puncture mass.
- `puncture_tagging.level_separation` controls how quickly the tagged region around each puncture grows between successive refinement levels.

Larger values increase the refined volume around the punctures.

## `ExtractionTagger`

`ExtractionTagger` ensures that each spherical extraction surface is covered by its requested refinement level. It tags a sphere slightly larger than the extraction radius, providing a buffer between the extraction surface and the refinement boundary. On successively coarser levels, the tagged radius is increased by powers of `level_separation` to maintain nested refinement regions.

This tagger is used by the `BinaryBH` example when Weyl extraction is enabled.

- `extraction_tagging.level_separation` controls how quickly the tagged region grows between successive refinement levels.

Larger values place coarser refinement boundaries farther from the extraction surface.

## `FixedGridsTagger`

`FixedGridsTagger` creates a fixed hierarchy of nested boxes centred on `geometry.center`. The linear size of the tagged region is halved on each successive level and does not depend on the evolved fields.

This tagger is used by the `KleinGordon` and `ScalarField` examples. It has no tagger-specific parameters; its placement and extent are determined by `geometry.center` and `geometry.prob_extent`.
