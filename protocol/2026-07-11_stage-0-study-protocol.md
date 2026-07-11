# Stage 0 protocol: human burn-wound to hypertrophic-scar transcriptomic transition

## Provisional title

**Transcriptomic programs across acute burn wounds and hypertrophic scars in a large human surgical cohort: a stage-associated analysis with single-cell localization**

This is provisional. It is not a claim that acute samples predict scars in the same patient.

## Primary objective

Identify stage-associated transcriptomic modules that distinguish normal skin, acute burn wounds and hypertrophic scars in GSE178411, and determine whether fixed modules have plausible cell-source localization in independent public burn or scar single-cell data.

## Primary data unit and inference boundary

- Each GSE178411 specimen is one bulk RNA-seq sample.
- Acute wounds and hypertrophic scars are cross-sectional groups unless patient identifiers prove longitudinal linkage.
- The principal inference is group association, not patient-level temporal evolution or scar prediction.

## Planned group definitions

1. Normal skin: uninjured donor skin.
2. Acute burn wound: retain metadata-defined early/late or other substage labels separately.
3. Hypertrophic scar: retain only samples explicitly designated HTS.

No group may be silently merged after results are observed.

## Pre-analysis metadata gate

Before expression testing, verify:

- Matrix is raw counts, not FPKM/TPM.
- Gene identifier namespace and duplicate-gene handling.
- Sample accession to count-matrix column mapping.
- Exact stage labels, sample sites, time-from-burn/wound timing, age/sex and patient linkage if available.
- Whether samples are multiple specimens from the same participant.

If there are repeated samples per patient, use a participant block/random effect or a patient-level aggregate; do not treat them as independent.

## Statistical flow

1. Filter low-count genes using a pre-specified edgeR rule.
2. TMM normalization and library-size/mean-variance QC.
3. MDS/PCA and sample correlation; no automatic exclusion.
4. limma-voom or edgeR quasi-likelihood models for pre-specified contrasts.
5. Primary contrasts: acute wound vs normal skin, HTS vs normal skin, HTS vs acute wound. If acute substages are adequately annotated, early vs late wound is a planned secondary comparison.
6. Program-level inference: Hallmark/GO/Reactome gene-set testing after contrast results; avoid hub-gene selection as a main endpoint.
7. Internal robustness: stage-stratified split or patient-aware resampling only if sample metadata permit. If no independent partition is defensible, label analysis as single-cohort discovery.
8. External scRNA reference mapping only for fixed bulk modules; no reverse gene selection after scRNA inspection.

## Stop/go gates

- **Gate 0:** counts and sample metadata reconstructable. Otherwise stop.
- **Gate 1:** stage grouping not fully explained by an unresolvable technical batch. Otherwise restrict claim or stop.
- **Gate 2:** at least one module/program is consistent across a pre-specified internal robustness check. Otherwise do not proceed to single-cell localization.
- **Gate 3:** cell-source localization is consistent in external scRNA data. Otherwise keep bulk result only and reduce claim.

## Prohibited claims

- No early molecular biomarker predicts HTS development.
- No causal stage trajectory at the individual patient level.
- No drug efficacy, therapeutic recommendation, or clinical prediction model.
- No single-cell cell-level p value used as independent patient-level validation.
