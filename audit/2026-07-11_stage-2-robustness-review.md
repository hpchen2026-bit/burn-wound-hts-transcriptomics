# Stage 2 robustness review: GSE178411 stage models

## Model results

All four stage contrasts show extensive stage-associated differential expression. Because this is surgical tissue with distinct wound/scar stages, gene counts are not the main scientific endpoint and will not be used as a novelty claim.

## Mandatory sensitivity concordance

| Contrast | Stage-only vs age/sex-adjusted logFC correlation | Stage-only vs non-E-restricted logFC correlation | FDR<0.05 shared with non-E restriction |
|---|---:|---:|---:|
| Early wound vs normal skin | 0.9999 | 0.9996 | 15,031 |
| Late wound vs normal skin | 0.9997 | 0.9990 | 14,930 |
| HTS vs normal skin | 0.9994 | 0.9719 | 10,260 |
| HTS vs late wound | 0.9981 | 0.9460 | 7,179 |

All shared significant genes retain the same effect direction in the non-E-restricted sensitivity model.

## Interpretation

The stage-associated bulk signal is robust to age/sex adjustment and to exclusion of the metadata-unlabelled E-prefix sample subset. This does not establish individual-level temporal progression or causality. It supports proceeding from single-gene lists to pre-specified **module-level** analysis.

## Decision

`gate_status = go_to_program_level_analysis`

Next stage will test Hallmark and Reactome/GO program changes across early wound, late wound and HTS, then identify modules that are stable in the two required sensitivity models. No hub-gene, classifier, drug-prediction or therapeutic claim will be generated.
