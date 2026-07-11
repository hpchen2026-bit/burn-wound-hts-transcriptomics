#!/usr/bin/env Rscript

# Stage 4B: unsupervised clustering and marker audit for GSE156326 human samples.
# Clusters are descriptive; no cell-level condition p values are used as confirmatory evidence.

suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(readr); library(ggplot2); library(tibble) })
root <- normalizePath(".")
x <- readRDS(file.path(root, "05-results/intermediate_GSE156326_human_stage4A_seurat.rds"))
table_dir <- file.path(root, "04-tables")
audit_dir <- file.path(root, "05-results/audit")
fig_dir <- file.path(root, "05-results/figures/stage_4B_GSE156326")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

x <- FindNeighbors(x, dims = 1:20, verbose = FALSE)
x <- FindClusters(x, resolution = 0.4, verbose = FALSE)
x <- RunUMAP(x, dims = 1:20, seed.use = 20260711, verbose = FALSE)
x <- JoinLayers(x)
markers <- FindAllMarkers(x, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox", verbose = FALSE)
write_csv(markers, file.path(table_dir, "GSE156326_human_stage4B_cluster_markers.csv"))
top_markers <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 15, with_ties = FALSE) %>% ungroup()
write_csv(top_markers, file.path(table_dir, "GSE156326_human_stage4B_top15_cluster_markers.csv"))
composition <- x@meta.data %>% as.data.frame() %>% count(seurat_clusters, GSM, condition, coarse_compartment_preliminary, name = "cells")
write_csv(composition, file.path(table_dir, "GSE156326_human_stage4B_cluster_sample_composition.csv"))
p <- DimPlot(x, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("GSE156326 human cells: descriptive clusters")
ggsave(file.path(fig_dir, "FigS_stage4B_UMAP_clusters.png"), p, width = 8, height = 6, dpi = 300)
p2 <- DimPlot(x, reduction = "umap", group.by = "condition") + ggtitle("GSE156326 human cells: tissue condition overlay")
ggsave(file.path(fig_dir, "FigS_stage4B_UMAP_condition.png"), p2, width = 8, height = 6, dpi = 300)
saveRDS(x, file.path(root, "05-results/intermediate_GSE156326_human_stage4B_clustered_seurat.rds"))
writeLines(c(
  "# GSE156326 human Stage 4B clustering and marker audit", "", "gate_status: marker_review_required", "",
  "- Clusters are exploratory descriptive partitions after minimal QC, normalization and PCA-neighbor graph construction.",
  "- Cluster marker tests are not used for tissue-condition inference because cells are nested within six source samples.",
  "- Final coarse cell labels require marker review and sample-composition review; condition comparison must be donor-aware pseudobulk if performed.",
  "- The only permitted next use is mapping frozen bulk module scores to reviewed cell compartments."
), file.path(audit_dir, "stage_4B_GSE156326_cluster_marker_audit.md"))
cat("GSE156326 Stage 4B clustering complete; marker review required.\n")
