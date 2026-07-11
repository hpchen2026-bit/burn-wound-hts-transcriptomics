#!/usr/bin/env Rscript

# Stage 3: Hallmark program-level analysis of fixed stage contrasts.
# Programs are tested in stage-only and non-E-restricted models; stable programs require compatible direction and FDR in both.

suppressPackageStartupMessages({
  library(limma)
  library(msigdbr)
  library(readr)
  library(dplyr)
  library(tibble)
})

root <- normalizePath(".")
table_dir <- file.path(root, "04-tables")
audit_dir <- file.path(root, "05-results/audit")

hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  filter(!is.na(ncbi_gene)) %>%
  transmute(gs_name, entrez_id = as.character(ncbi_gene)) %>% distinct()
sets <- split(hallmark$entrez_id, hallmark$gs_name)
sets <- lapply(sets, unique)
contrasts <- c("early_wound_vs_normal_skin", "late_wound_vs_normal_skin", "hypertrophic_scar_vs_normal_skin", "hypertrophic_scar_vs_late_wound")

run_camera <- function(mode, contrast) {
  tab <- read_csv(file.path(table_dir, paste0("GSE178411_", mode, "_", contrast, "_limma_voom.csv")), show_col_types = FALSE) %>%
    filter(!is.na(entrez_id), !is.na(t)) %>% mutate(entrez_id = as.character(entrez_id)) %>% distinct(entrez_id, .keep_all = TRUE)
  stat <- tab$t; names(stat) <- tab$entrez_id
  ans <- cameraPR(statistic = stat, index = sets, inter.gene.cor = 0.01, sort = TRUE) %>%
    rownames_to_column("gs_name") %>% as_tibble() %>%
    mutate(mode = mode, contrast = contrast)
  write_csv(ans, file.path(table_dir, paste0("GSE178411_", mode, "_", contrast, "_Hallmark_cameraPR.csv")))
  ans
}

all_res <- bind_rows(lapply(contrasts, function(cn) bind_rows(run_camera("stage_only", cn), run_camera("nonE_restricted", cn))))
write_csv(all_res, file.path(table_dir, "GSE178411_stage3_Hallmark_all_models.csv"))

stable <- all_res %>% select(gs_name, mode, contrast, Direction, PValue, FDR, NGenes) %>%
  tidyr::pivot_wider(names_from = mode, values_from = c(Direction, PValue, FDR, NGenes)) %>%
  mutate(stable_FDR05 = FDR_stage_only < 0.05 & FDR_nonE_restricted < 0.05 & Direction_stage_only == Direction_nonE_restricted) %>%
  arrange(contrast, desc(stable_FDR05), FDR_stage_only)
write_csv(stable, file.path(table_dir, "GSE178411_stage3_Hallmark_stable_programs.csv"))

summary <- stable %>% group_by(contrast) %>% summarise(stable_programs_FDR05 = sum(stable_FDR05), .groups = "drop")
write_csv(summary, file.path(table_dir, "GSE178411_stage3_Hallmark_stable_summary.csv"))
writeLines(c(
  "# GSE178411 Stage 3 Hallmark program analysis", "", "gate_status: review_required", "",
  "- Hallmark sets were tested from fixed stage-contrast t statistics using cameraPR.",
  "- Stable program status requires FDR < 0.05 and identical direction in both stage-only and non-E-restricted models.",
  "- This provides pathway-level robustness, not individual-level temporal causality or a therapeutic mechanism."
), file.path(audit_dir, "stage_3_GSE178411_Hallmark_programs.md"))
cat("Stage 3 Hallmark analysis complete; review required.\n")
