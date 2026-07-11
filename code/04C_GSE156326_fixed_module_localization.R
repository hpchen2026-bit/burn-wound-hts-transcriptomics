#!/usr/bin/env Rscript

# Stage 4C: map pre-frozen bulk programs to reviewed GSE156326 human cell compartments.
# Scores are descriptive cell-source localization. There is no cell-level tissue-condition hypothesis test.

suppressPackageStartupMessages({ library(Seurat); library(msigdbr); library(dplyr); library(readr); library(ggplot2); library(tibble) })
root <- normalizePath(".")
x <- readRDS(file.path(root, "05-results/intermediate_GSE156326_human_stage4B_clustered_seurat.rds"))
table_dir <- file.path(root, "04-tables")
audit_dir <- file.path(root, "05-results/audit")
fig_dir <- file.path(root, "05-results/figures/stage_4C_GSE156326")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Labels are based on the saved cluster-marker audit. Ambiguous fibroblast clusters remain broad fibroblast subtypes.
labels <- c(
  `0` = "activated_fibroblast_POSTN_ADAM12", `1` = "fibroblast_CXCL12_APOE", `2` = "fibroblast_CXCL14_CFD",
  `3` = "basal_keratinocyte", `4` = "T_cell", `5` = "vascular_endothelial_ACKR1", `6` = "myeloid_macrophage",
  `7` = "vascular_endothelial_PECAM1", `8` = "pericyte_smc", `9` = "fibroblast_PCOLCE2_IGFBP6",
  `10` = "fibroblast_PI16_ELN", `11` = "fibroblast_TWIST1", `12` = "fibroblast_CILP",
  `13` = "lymphatic_endothelial", `14` = "melanocyte", `15` = "endothelial_IGF2", `16` = "schwann_like"
)
x$reviewed_cell_compartment <- unname(labels[as.character(x$seurat_clusters)])

reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>% filter(!is.na(gene_symbol))
modules <- c(
  ECM_organization = "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  collagen_formation = "REACTOME_COLLAGEN_FORMATION",
  collagen_biosynthesis = "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES",
  neutrophil_degranulation = "REACTOME_NEUTROPHIL_DEGRANULATION",
  innate_immune_system = "REACTOME_INNATE_IMMUNE_SYSTEM",
  cell_cycle_mitotic = "REACTOME_CELL_CYCLE_MITOTIC"
)
for (module_name in names(modules)) {
  genes <- unique(reactome$gene_symbol[reactome$gs_name == modules[[module_name]]])
  genes <- intersect(genes, rownames(x))
  if (length(genes) < 5) stop("Insufficient matched genes for module: ", module_name)
  x <- AddModuleScore(x, features = list(genes), name = paste0(module_name, "_score"), search = FALSE)
}
score_cols <- paste0(names(modules), "_score1")
meta <- x@meta.data %>% as.data.frame() %>% rownames_to_column("cell")
summary <- meta %>% group_by(GSM, condition, reviewed_cell_compartment) %>%
  summarise(cells = n(), across(all_of(score_cols), mean), .groups = "drop")
write_csv(summary, file.path(table_dir, "GSE156326_stage4C_module_score_by_sample_compartment.csv"))

# A second summary prevents cell count from being misread as donor replication.
condition_summary <- summary %>% group_by(condition, reviewed_cell_compartment) %>%
  summarise(samples = n_distinct(GSM), total_cells = sum(cells), across(all_of(score_cols), mean), .groups = "drop")
write_csv(condition_summary, file.path(table_dir, "GSE156326_stage4C_module_score_by_condition_compartment_descriptive.csv"))

plot_df <- condition_summary %>% select(condition, reviewed_cell_compartment, ECM_organization_score1, samples, total_cells) %>% filter(total_cells >= 30)
p <- ggplot(plot_df, aes(reorder(reviewed_cell_compartment, ECM_organization_score1), ECM_organization_score1, colour = condition, size = total_cells)) +
  geom_point(position = position_dodge(width = 0.45)) + coord_flip() + theme_classic(base_size = 12) + labs(title = "GSE156326: ECM program score by reviewed compartment", x = NULL, y = "AddModuleScore (descriptive)", size = "cells")
ggsave(file.path(fig_dir, "FigS_stage4C_ECM_score_by_compartment.png"), p, width = 10, height = 7, dpi = 300)

saveRDS(x, file.path(root, "05-results/intermediate_GSE156326_human_stage4C_module_scored.rds"))
writeLines(c(
  "# GSE156326 Stage 4C fixed bulk-module localization", "", "gate_status: review_required", "",
  "- Reactome program definitions were fixed from Stage 3 before this single-cell scoring step.",
  "- Scores are descriptive and aggregated by source sample and reviewed compartment.",
  "- No cell-level condition p value, differential expression or causal mechanism claim is reported.",
  "- GSE156326 supports HTS cell-source localization only; it cannot validate acute burn-wound stage dynamics."
), file.path(audit_dir, "stage_4C_GSE156326_fixed_module_localization.md"))
cat("GSE156326 Stage 4C module localization complete; review required.\n")
