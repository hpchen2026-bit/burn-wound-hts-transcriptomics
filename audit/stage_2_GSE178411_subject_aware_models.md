# GSE178411 Stage 2 subject-aware models

gate_status: review_required

- Three pre-specified model families were run: stage-only, age/sex adjusted, and non-E-prefix restricted sensitivity.
- Subject correlation was estimated with limma duplicateCorrelation in each model family.
- E-prefix is an unexplained source metadata pattern, not a labelled sequencing batch; non-E restriction is a robustness test, not batch correction.
- Downstream module analysis is allowed only for contrasts with compatible direction and broad signal under the sensitivity models.
