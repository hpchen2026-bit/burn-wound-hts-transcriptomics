# Burn wound–HTS transcriptomics

Reproducible analysis code and derived results for a public-data study of stage-associated transcriptomic programs across human burn wounds and hypertrophic scars (HTS).

## Study summary

The primary bulk RNA-sequencing cohort is GEO **GSE178411**. The primary analysis includes normal skin (n=24), early wound (n=22), late wound (n=29), and HTS (n=28) specimens. The analysis uses participant-aware bulk models and prespecified sensitivity analyses, then maps fixed programs to human HTS single-cell data (GSE156326) and a human wound spatial reference (GSE241124).

This is a **cross-sectional secondary analysis**. It does not predict HTS, establish individual temporal trajectories, identify causal mechanisms, or evaluate treatments.

## Repository contents

- `protocol/` — study protocol and group-lock records.
- `code/` — R and Python scripts for QC, models, program analyses, external-reference analysis and figures.
- `tables/` — derived analytical tables used in the manuscript.
- `audit/` — auditable analysis and interpretation boundaries.
- `data_sources/` — GEO source URLs, download manifest and checksums.

## Data access

The repository does **not** redistribute GEO archives. Obtain the source data from:

- GSE178411: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178411
- GSE156326: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156326
- GSE241124: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241124

The `data_sources/` files record the exact source URLs and SHA256 checksums used for the local analysis archive.

## Reproduction environment

Primary bulk and single-cell scripts were run with R 4.6.1 and the `edgeR`, `limma`, `msigdbr`, `Seurat`, `readr`, `dplyr`, and `ggplot2` packages. Spatial H5 audit and scoring scripts use Python 3 with `h5py`, `scipy`, `pandas`, `numpy`, and `matplotlib`.

Run scripts from the repository root in numerical order. Source archives should be placed under the paths described in `data_sources/download_manifest.csv`; scripts do not download or expose any credentials.

## Transparency statement

The repository is a versioned reproducibility archive for an analysis of public, deidentified data. It contains no raw patient-level files, credentials, or private institutional material.

## Citation

A citable release DOI will be added after the manuscript and archive are finalized.
