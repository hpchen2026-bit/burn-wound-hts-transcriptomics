#!/usr/bin/env Rscript

# Stage 2: subject-aware stage models and mandatory confounding sensitivities.
# This script reports associations only; no temporal, causal or predictive claim is created.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(readr)
  library(dplyr)
  library(tibble)
})

root <- normalizePath(".")
obj <- readRDS(file.path(root, "05-results/intermediate_GSE178411_stage1_DGEList.rds"))
base_y <- obj$DGEList
base_manifest <- obj$manifest
base_y$samples$subject <- base_manifest$subject
base_y$samples$matrix_column_id <- base_manifest$matrix_column_id
base_y$samples$analysis_group <- base_manifest$analysis_group
base_y$samples$age <- suppressWarnings(as.numeric(base_manifest$age))
base_y$samples$sex <- base_manifest$sex

table_dir <- file.path(root, "04-tables")
audit_dir <- file.path(root, "05-results/audit")

make_contrast_matrix <- function(design) {
  cn <- colnames(design)
  groups <- c("normal_skin", "early_wound", "late_wound", "hypertrophic_scar")
  want <- list(
    early_wound_vs_normal_skin = c(early_wound = 1, normal_skin = -1),
    late_wound_vs_normal_skin = c(late_wound = 1, normal_skin = -1),
    hypertrophic_scar_vs_normal_skin = c(hypertrophic_scar = 1, normal_skin = -1),
    hypertrophic_scar_vs_late_wound = c(hypertrophic_scar = 1, late_wound = -1)
  )
  mat <- sapply(want, function(spec) {
    v <- rep(0, length(cn)); names(v) <- cn
    for (g in names(spec)) v[paste0("group", g)] <- spec[[g]]
    v
  })
  colnames(mat) <- names(want)
  mat
}

run_model <- function(y, manifest, mode) {
  grp <- factor(manifest$analysis_group, levels = c("normal_skin", "early_wound", "late_wound", "hypertrophic_scar"))
  if (any(is.na(grp))) stop("Unknown group in ", mode)
  if (mode == "stage_only") {
    design <- model.matrix(~ 0 + grp)
  } else if (mode == "age_sex_adjusted") {
    design <- model.matrix(~ 0 + grp + age + sex, data = manifest)
  } else if (mode == "nonE_restricted") {
    design <- model.matrix(~ 0 + grp)
  } else stop("Unknown mode")
  colnames(design) <- sub("grp", "group", colnames(design), fixed = TRUE)

  # Subject block is retained even though many subjects are singletons; consensus correlation uses informative repeated blocks.
  corfit <- duplicateCorrelation(cpm(y, log = TRUE, prior.count = 2), design, block = manifest$subject)
  v <- voom(y, design, plot = FALSE, block = manifest$subject, correlation = corfit$consensus)
  fit <- lmFit(v, design, block = manifest$subject, correlation = corfit$consensus)
  cm <- make_contrast_matrix(design)
  fit2 <- eBayes(contrasts.fit(fit, cm), robust = TRUE, trend = TRUE)

  summary_rows <- list()
  for (contrast in colnames(cm)) {
    tt <- topTable(fit2, coef = contrast, number = Inf, sort.by = "P") %>% rownames_to_column("entrez_id")
    write_csv(tt, file.path(table_dir, paste0("GSE178411_", mode, "_", contrast, "_limma_voom.csv")))
    summary_rows[[contrast]] <- tibble(mode = mode, contrast = contrast, samples = ncol(y), genes = nrow(tt), FDR05 = sum(tt$adj.P.Val < 0.05), up_FDR05 = sum(tt$adj.P.Val < 0.05 & tt$logFC > 0), down_FDR05 = sum(tt$adj.P.Val < 0.05 & tt$logFC < 0), consensus_subject_correlation = corfit$consensus)
  }
  bind_rows(summary_rows)
}

# stage-only full primary set
m1 <- base_manifest
idx1 <- match(m1$matrix_column_id, colnames(base_y)); y1 <- base_y[, idx1]
s1 <- run_model(y1, m1, "stage_only")

# age/sex complete cases; one late-wound sample with unavailable age is excluded only from this sensitivity model
m2 <- base_manifest %>% filter(!is.na(suppressWarnings(as.numeric(age))), sex %in% c("male", "female")) %>% mutate(age = as.numeric(age), sex = factor(sex))
idx2 <- match(m2$matrix_column_id, colnames(base_y)); y2 <- base_y[, idx2]
s2 <- run_model(y2, m2, "age_sex_adjusted")

# E-prefix restricted sensitivity: tests whether HTS signal depends on metadata-unlabelled E-coded specimen subset
m3 <- base_manifest %>% filter(!startsWith(matrix_column_id, "E"))
idx3 <- match(m3$matrix_column_id, colnames(base_y)); y3 <- base_y[, idx3]
s3 <- run_model(y3, m3, "nonE_restricted")

summary <- bind_rows(s1, s2, s3)
write_csv(summary, file.path(table_dir, "GSE178411_stage2_model_summary.csv"))
writeLines(c(
  "# GSE178411 Stage 2 subject-aware models", "", "gate_status: review_required", "",
  "- Three pre-specified model families were run: stage-only, age/sex adjusted, and non-E-prefix restricted sensitivity.",
  "- Subject correlation was estimated with limma duplicateCorrelation in each model family.",
  "- E-prefix is an unexplained source metadata pattern, not a labelled sequencing batch; non-E restriction is a robustness test, not batch correction.",
  "- Downstream module analysis is allowed only for contrasts with compatible direction and broad signal under the sensitivity models."
), file.path(audit_dir, "stage_2_GSE178411_subject_aware_models.md"))
capture.output(sessionInfo(), file = file.path(audit_dir, "stage_2_GSE178411_session_info.txt"))
cat("GSE178411 Stage 2 models complete; review required.\n")
