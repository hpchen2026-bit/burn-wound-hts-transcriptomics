# GSE178411 Stage 1 count-matrix QC

gate_status: review_required

- Counts were read as 28,395 Entrez-ID rows and 108 GEO sample columns; the missing header label for the gene-ID field was handled explicitly.
- Primary QC includes 103 specimens: normal skin, early wound, late wound and HTS. Chronic wound and normal scar remain excluded from primary analysis.
- Low-expression filtering and TMM normalization were performed. No differential-expression test was run.
- PCA and library-size metrics must be reviewed with stage, subject repetition, age, sex, location and days-since-injury metadata before fixing the limma-voom design.
