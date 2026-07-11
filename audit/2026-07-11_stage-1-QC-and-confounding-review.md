# Stage 1 QC and confounding review: GSE178411

## Technical QC

- 28,395 Entrez rows; 18,829 retained after low-expression filtering.
- 103 primary cross-sectional specimens: normal skin 24, early wound 22, late wound 29, HTS 28.
- TMM-normalized logCPM distributions have comparable medians and no isolated sample-level distribution failure.
- PCA separates normal skin, wound stages and HTS; this is a biological-stage signal candidate, not proof of causality.

## Critical design issue

Sample IDs with the `E` prefix are imbalanced across groups: HTS 19/28, late wound 4/29, normal skin 1/24, early wound 0/22. The metadata do not label this prefix as a sequencing batch. It may represent a collection/cohort difference and can confound HTS contrasts.

Other observed structural differences:

- HTS median age 19 years; normal skin 30.5 years; late wound 35.5 years.
- 75 source subjects for 103 primary samples; 23 subjects contributed multiple specimens.
- 14 late-wound, 11 early-wound and 2 HTS specimens have a same-subject normal-skin sample.

## Gate decision

`gate_status = go_with_mandatory_sensitivity`

Proceed to expression models only under all three analyses below:

1. Stage-only model with subject correlation.
2. Age/sex-adjusted complete-case sensitivity model with subject correlation.
3. Non-`E` sample-ID restricted sensitivity, especially for HTS contrasts.

A stage-associated module is eligible for the paper only if direction and broad functional theme are robust to these sensitivity models. If HTS results disappear under non-`E` restriction, they cannot be interpreted as robust scar biology.
