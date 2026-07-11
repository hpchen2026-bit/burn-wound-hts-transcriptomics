#!/usr/bin/env Rscript

# Stage 4A diagnostic: inspect donor/sample structure before any condition-level interpretation.
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(readr); library(tibble) })
root <- normalizePath(".")
x <- readRDS(file.path(root, "05-results/intermediate_GSE156326_human_stage4A_seurat.rds"))
fig_dir <- file.path(root, "05-results/figures/stage_4A_GSE156326")
table_dir <- file.path(root, "04-tables")
pcs <- Embeddings(x, "pca")[, 1:2, drop = FALSE] %>% as.data.frame() %>% rownames_to_column("cell") %>%
  mutate(GSM = x$GSM, condition = x$condition)
centroids <- pcs %>% group_by(GSM, condition) %>% summarise(PC_1 = mean(PC_1), PC_2 = mean(PC_2), cells = n(), .groups = "drop")
write_csv(centroids, file.path(table_dir, "GSE156326_human_stage4A_PCA_centroids.csv"))
p <- ggplot(centroids, aes(PC_1, PC_2, colour = condition, label = GSM, size = cells)) + geom_point() + geom_text(vjust = -0.8) + theme_classic(base_size = 13) + labs(title = "GSE156326: sample centroids in PCA space")
ggsave(file.path(fig_dir, "FigS_stage4A_PCA_sample_centroids.png"), p, width = 8, height = 6, dpi = 300)
p2 <- DimPlot(x, reduction = "pca", group.by = "GSM") + ggtitle("GSE156326 human cells: PCA by source sample")
ggsave(file.path(fig_dir, "FigS_stage4A_PCA_by_sample.png"), p2, width = 9, height = 7, dpi = 300)
cat("GSE156326 sample diagnostics complete.\n")
