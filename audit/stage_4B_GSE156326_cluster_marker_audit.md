# GSE156326 human Stage 4B clustering and marker audit

gate_status: marker_review_required

- Clusters are exploratory descriptive partitions after minimal QC, normalization and PCA-neighbor graph construction.
- Cluster marker tests are not used for tissue-condition inference because cells are nested within six source samples.
- Final coarse cell labels require marker review and sample-composition review; condition comparison must be donor-aware pseudobulk if performed.
- The only permitted next use is mapping frozen bulk module scores to reviewed cell compartments.
