# Stage 4 external single-cell reference selection audit

## GSE163446: excluded

GEO identifies GSE163446 as single-cell sequencing from a **tendon injury site after severe burn/tenotomy injury with or without sciatic neurectomy**. It is not a human burn-skin reference. It will not be used for cell-source localization in this human burn-wound/HTS project.

## GSE156326: eligible HTS reference

GSE156326 includes processed human normal-skin and human scar single-cell matrices (three samples each), plus separate mouse scar experiments. Only the six human samples are eligible. The mouse samples are excluded.

## Required before Stage 4 completion

1. Archive and SHA256-verify GSE156326 source archive and metadata inside this project.
2. Identify and audit a **human burn-skin** single-cell dataset before adding it to the project. A non-skin, non-human or mechanically different injury model cannot be substituted.
3. Use external scRNA only to assign likely cell sources to bulk modules frozen at Stage 3; it cannot be used to discover or revise the modules.
