#!/usr/bin/env Rscript

# Stage 4A: donor-aware QC and coarse compartment annotation for GSE156326 human samples only.
# Mouse samples are never read. This script does not modify immutable raw archives.

suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})
root <- normalizePath(".")
in_dir <- file.path(root, "05-results/intermediate_references/GSE156326_human")
table_dir <- file.path(root, "04-tables")
audit_dir <- file.path(root, "05-results/audit")
fig_dir <- file.path(root, "05-results/figures/stage_4A_GSE156326")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c(
  GSM4729097 = "human_skin_1", GSM4729098 = "human_skin_2", GSM4729099 = "human_skin_3",
  GSM4729100 = "human_scar_1", GSM4729101 = "human_scar_2", GSM4729102 = "human_scar_3"
)
condition <- ifelse(grepl("skin", samples), "normal_skin", "hypertrophic_scar")

read_one <- function(gsm, label, cond) {
  stem <- paste0(gsm, "_", label)
  mtx <- file.path(in_dir, paste0(stem, "_matrix.mtx.gz"))
  feature <- list.files(in_dir, pattern = paste0("^", stem, "_(features|genes)\\.tsv\\.gz$"), full.names = TRUE)
  barcode <- file.path(in_dir, paste0(stem, "_barcodes.tsv.gz"))
  stopifnot(file.exists(mtx), length(feature) == 1, file.exists(barcode))
  x <- ReadMtx(mtx = mtx, features = feature, cells = barcode, feature.column = 2)
  o <- CreateSeuratObject(x, project = gsm, min.cells = 3, min.features = 0)
  o$GSM <- gsm; o$sample_label <- label; o$condition <- cond
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  o
}

objs <- Map(read_one, names(samples), unname(samples), condition)
names(objs) <- names(samples)
pre <- bind_rows(lapply(objs, function(o) tibble(GSM = unique(o$GSM), condition = unique(o$condition), cells_before = ncol(o), median_features_before = median(o$nFeature_RNA), median_counts_before = median(o$nCount_RNA), median_mt_before = median(o$percent.mt))))

# Transparent minimal QC. High-feature cells are flagged for later doublet review rather than removed automatically.
objs <- lapply(objs, function(o) subset(o, subset = nFeature_RNA >= 200 & percent.mt < 25))
post <- bind_rows(lapply(objs, function(o) tibble(GSM = unique(o$GSM), condition = unique(o$condition), cells_after = ncol(o), median_features_after = median(o$nFeature_RNA), median_counts_after = median(o$nCount_RNA), median_mt_after = median(o$percent.mt), high_feature_flag_n = sum(o$nFeature_RNA > 7500))))
qc <- left_join(pre, post, by = c("GSM", "condition"))
write_csv(qc, file.path(table_dir, "GSE156326_human_stage4A_QC_by_sample.csv"))
if (any(qc$cells_after < 100)) stop("A human sample retains fewer than 100 cells after QC; stop for review")

merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs), project = "GSE156326_human")
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)

# Coarse compartments are a transparent marker-score aid; they require marker review before final labels.
markers <- list(
  fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"),
  keratinocyte = c("KRT14", "KRT5", "KRT1", "KRT10", "EPCAM"),
  immune = c("PTPRC", "LYZ", "CD3D", "MS4A1", "NKG7"),
  endothelial = c("PECAM1", "VWF", "KDR", "EMCN"),
  pericyte_smc = c("RGS5", "MCAM", "CSPG4", "ACTA2", "TAGLN")
)
markers <- lapply(markers, function(x) intersect(x, rownames(merged)))
for (marker_name in names(markers)) {
  merged <- AddModuleScore(merged, features = list(markers[[marker_name]]), name = paste0(marker_name, "_score"), search = FALSE)
}
score_cols <- paste0(names(markers), "_score1")
score_mat <- FetchData(merged, vars = score_cols)
merged$coarse_compartment_preliminary <- names(markers)[max.col(as.matrix(score_mat), ties.method = "first")]

p <- DimPlot(merged, reduction = "pca", group.by = "condition") + ggtitle("GSE156326 human cells: PCA after minimal QC")
ggsave(file.path(fig_dir, "FigS_stage4A_PCA_by_condition.png"), p, width = 8, height = 6, dpi = 300)
p2 <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "GSM", pt.size = 0)
ggsave(file.path(fig_dir, "FigS_stage4A_QC_violin.png"), p2, width = 12, height = 7, dpi = 300)

compartment_summary <- merged@meta.data %>% as.data.frame() %>% count(GSM, condition, coarse_compartment_preliminary, name = "cells")
write_csv(compartment_summary, file.path(table_dir, "GSE156326_human_stage4A_preliminary_compartments.csv"))
saveRDS(merged, file.path(root, "05-results/intermediate_GSE156326_human_stage4A_seurat.rds"))
writeLines(c(
  "# GSE156326 human Stage 4A single-cell QC", "", "gate_status: review_required", "",
  "- Only six human matrices were read; all mouse matrices were excluded at input.",
  "- Minimal QC retained cells with >=200 detected features and <25% mitochondrial reads.",
  "- High-feature cells were flagged, not removed; formal doublet assessment remains pending.",
  "- Coarse compartment labels are marker-score aids and must not be treated as final cell annotations without marker review.",
  "- This reference is used for HTS cell-source localization only and does not establish burn-stage dynamics."
), file.path(audit_dir, "stage_4A_GSE156326_human_scRNA_QC.md"))
capture.output(sessionInfo(), file = file.path(audit_dir, "stage_4A_GSE156326_session_info.txt"))
cat("GSE156326 human Stage 4A complete; review required.\n")
