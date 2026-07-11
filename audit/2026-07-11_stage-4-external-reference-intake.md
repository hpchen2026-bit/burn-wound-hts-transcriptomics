# Stage 4 external-reference intake audit

## Archived references

### GSE156326: human HTS single-cell reference

- Archive retained at `02-data-raw/GEO/GSE156326/original/GSE156326_RAW.tar`.
- Source metadata and SHA256 are retained separately.
- Initial working extraction contains only 18 files for the six eligible human matrices: three normal-skin and three human-scar samples.
- All seven mouse samples were excluded before any analysis.
- Role: identify likely human scar cell sources for fixed bulk ECM/collagen modules. It is not an independent burn cohort and will not validate acute burn-stage effects.

### GSE241124: human acute-wound spatial reference

- Archive retained at `02-data-raw/GEO/GSE241124/original/GSE241124_RAW.tar`.
- Initial working extraction contains 16 filtered H5 matrices; image ZIP files remain preserved in the raw archive and are not needed for initial count/spot audit.
- Source metadata yields 16 samples: four each of Skin, Wound1, Wound7 and Wound30, across four donors.
- Role: localize fixed acute-wound modules in normal skin and early/remodeling wound spatial contexts. This is a healthy-volunteer skin-wounding reference, not a burn cohort.

## Explicit exclusion

GSE163446 remains excluded because it is a tendon injury model after severe burn/tenotomy, not human burn skin.

## Gate status

`gate_status = software_setup_and_reference_QC_pending`

Next: install/test the single-cell and H5 environment, then conduct donor-aware QC and annotation for GSE156326 human samples. Spatial scoring in GSE241124 follows only after the bulk module list is frozen and single-cell QC is complete.
