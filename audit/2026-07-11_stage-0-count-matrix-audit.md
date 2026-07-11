# GSE178411 Stage 0 count-matrix audit

- Matrix header columns: 108 sample IDs.
- Data rows: 28395 gene rows.
- Data field-count distribution: {109: 28395}. Each row has one numeric gene identifier plus 108 counts.
- Header has no explicit gene-ID label; the first data field is treated as the gene identifier, while all 108 header fields are sample IDs.
- First gene identifiers: 100287102, 653635, 102466751, 100302278, 645520.
- All-zero gene rows: 1325.
- Sample metadata rows: 108; exact count-header mapping: 108/108.
- Exact source-stage labels: {'Early Wound': 22, 'Late wound': 29, 'Normal scar': 2, 'Normal skin': 24, 'HTS': 28, 'Chronic wound': 3}.
- Unique source subjects: 75; subjects with >1 specimen: 23.

## Gate 0 status

`go_to_metadata_structure_audit`

The counts matrix is structurally usable and maps to all GEO sample titles. The discrepancy between the GEO series summary groups (normal skin 26, acute wounds 54, HTS 30) and exact sample-title labels (normal skin 24, normal scar 2; early/late/chronic wound 54; HTS 28, normal scar 2) must be explicitly resolved before defining analysis groups. No expression analysis has been performed.
