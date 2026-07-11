# Stage 4 integration review: external cell-source and spatial context

## Reference eligibility and boundaries

- GSE156326: six human single-cell matrices only, comprising three normal-skin and three hypertrophic-scar samples. Mouse matrices were excluded before QC.
- GSE241124: 16 human spatial matrices, four each from Skin, Wound1, Wound7 and Wound30 across four donors. This is a healthy-volunteer wound-healing reference, not a burn cohort.
- GSE163446 was excluded: it is a burn/tenotomy tendon-injury model, not human burn skin.

## GSE156326 cell-source localization

Fixed Reactome ECM/collagen programs from the GSE178411 bulk analysis are descriptively highest in scar fibroblast compartments, not keratinocyte or immune compartments. In HTS samples, the ECM program score is highest in:

1. `fibroblast_PCOLCE2_IGFBP6`: 0.322; present in all three scar samples.
2. `fibroblast_PI16_ELN`: 0.297; present in all three scar samples.
3. `activated_fibroblast_POSTN_ADAM12`: 0.265; present in all three scar samples.

This supports fibroblast-associated localization of the fixed ECM/collagen program. It is descriptive cell-source evidence, not an independent statistical replication or a causal cell-state claim.

## GSE241124 spatial context

At the donor/sample level, all four conditions have four spatial samples. The descriptive mean ECM program scores are Skin 0.255, Wound1 0.181, Wound7 0.273 and Wound30 0.365. Corresponding collagen-formation scores are 0.327, 0.205, 0.361 and 0.505.

Thus, in an independent human non-burn wound reference, ECM/collagen programs are lower at Wound1 and higher by Wound30. This is compatible with wound-remodeling context, but cannot validate burn-specific or HTS-specific kinetics.

## Technical note

The R `hdf5r` package could not compile against the current Homebrew HDF5 API. Spatial H5 matrices were therefore audited and scored reproducibly with Python `h5py` and `scipy`; all 16 matrices mapped fully to the subset of publicly annotated spots.

## Decision

`gate_status = go_to_figure_and_manuscript_planning`

The evidence chain is now:

`GSE178411 subject-aware bulk stage program → GSE156326 human HTS fibroblast-source localization → GSE241124 human wound spatial context`

The paper remains a cross-sectional, public-data study. It must not claim individual temporal progression, burn-scar prediction, therapeutic target validation or causal mechanism.
