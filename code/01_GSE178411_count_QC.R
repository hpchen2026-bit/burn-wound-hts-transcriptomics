#!/usr/bin/env Rscript

# Stage 1: count-matrix QC and design-structure audit for GSE178411.
# No differential-expression contrast is tested in this script.

suppressPackageStartupMessages({
  library(edgeR)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

root <- normalizePath(".")
counts_path <- file.path(root, "02-data-raw/GEO/original/GSE178411_counts.txt.gz")
manifest_path <- file.path(root, "04-tables/GSE178411_sample_manifest_v0.2.csv")
table_dir <- file.path(root, "04-tables")
audit_dir <- file.path(root, "05-results/audit")
fig_dir <- file.path(root, "05-results/figures/stage_1_GSE178411_QC")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# GEO counts header lacks the first gene-ID column label.
con <- gzfile(counts_path, "rt")
header <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
raw <- read.delim(con, header = FALSE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
close(con)
if (ncol(raw) != length(header) + 1L) stop("Unexpected count matrix dimensions")
colnames(raw) <- c("entrez_id", header)

manifest <- read_csv(manifest_path, show_col_types = FALSE) %>%
  filter(primary_analysis_role == "primary_cross_sectional") %>%
  mutate(analysis_group = factor(analysis_group, levels = c("normal_skin", "early_wound", "late_wound", "hypertrophic_scar"))) %>%
  arrange(match(matrix_column_id, header))
stopifnot(nrow(manifest) == 103L, all(manifest$matrix_column_verified == TRUE), all(manifest$matrix_column_id %in% header))

count_mat <- as.matrix(raw[, manifest$matrix_column_id, drop = FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- as.character(raw$entrez_id)
if (anyDuplicated(rownames(count_mat))) stop("Duplicate Entrez identifiers require explicit collapsing before analysis")

# Filter only after sample groups are locked.
y <- DGEList(counts = count_mat, samples = as.data.frame(manifest))
keep <- filterByExpr(y, group = manifest$analysis_group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
logcpm <- cpm(y, log = TRUE, prior.count = 2)

# QC metrics and PCA/MDS are diagnostic only.
pc <- prcomp(t(logcpm), center = TRUE, scale. = FALSE)
pc_df <- as.data.frame(pc$x[, 1:2]) %>% rownames_to_column("matrix_column_id") %>%
  left_join(manifest %>% select(matrix_column_id, GSM, analysis_group, subject, repeat_subject_sample_count, age, sex, location), by = "matrix_column_id")
var_explained <- summary(pc)$importance[2, 1:2] * 100
qc <- manifest %>%
  transmute(matrix_column_id, GSM, analysis_group, subject, repeat_subject_sample_count, age, sex, location,
            library_size = y$samples$lib.size, TMM_norm_factor = y$samples$norm.factors,
            logCPM_median = apply(logcpm, 2, median), logCPM_IQR = apply(logcpm, 2, IQR),
            PC1 = pc_df$PC1[match(matrix_column_id, pc_df$matrix_column_id)],
            PC2 = pc_df$PC2[match(matrix_column_id, pc_df$matrix_column_id)])
write_csv(qc, file.path(table_dir, "GSE178411_stage1_QC_metrics.csv"))
write_csv(manifest, file.path(table_dir, "GSE178411_primary_cross_sectional_manifest.csv"))
write_csv(tibble(total_entrez_rows = nrow(raw), retained_after_filter = nrow(y), samples = ncol(y), zero_rows_before_filter = sum(rowSums(count_mat) == 0)), file.path(table_dir, "GSE178411_stage1_filter_summary.csv"))

p <- ggplot(pc_df, aes(PC1, PC2, colour = analysis_group, shape = as.factor(repeat_subject_sample_count))) +
  geom_point(size = 3.2) +
  labs(title = "GSE178411 primary stage groups: logCPM PCA", x = sprintf("PC1 (%.1f%%)", var_explained[1]), y = sprintf("PC2 (%.1f%%)", var_explained[2]), colour = "stage", shape = "specimens per subject") +
  theme_classic(base_size = 13)
ggsave(file.path(fig_dir, "FigS_stage1_logCPM_PCA.png"), p, width = 8.5, height = 6.5, dpi = 300)

png(file.path(fig_dir, "FigS_stage1_logCPM_boxplot.png"), width = 3200, height = 1600, res = 250)
boxplot(logcpm, las = 2, col = as.integer(manifest$analysis_group), ylab = "log2 CPM", main = "GSE178411 primary stage groups: TMM-normalized distributions")
dev.off()

writeLines(c(
  "# GSE178411 Stage 1 count-matrix QC",
  "", "gate_status: review_required", "",
  "- Counts were read as 28,395 Entrez-ID rows and 108 GEO sample columns; the missing header label for the gene-ID field was handled explicitly.",
  "- Primary QC includes 103 specimens: normal skin, early wound, late wound and HTS. Chronic wound and normal scar remain excluded from primary analysis.",
  "- Low-expression filtering and TMM normalization were performed. No differential-expression test was run.",
  "- PCA and library-size metrics must be reviewed with stage, subject repetition, age, sex, location and days-since-injury metadata before fixing the limma-voom design."
), file.path(audit_dir, "stage_1_GSE178411_count_QC.md"))
capture.output(sessionInfo(), file = file.path(audit_dir, "stage_1_GSE178411_session_info.txt"))
saveRDS(list(DGEList = y, manifest = manifest, logCPM = logcpm), file.path(root, "05-results/intermediate_GSE178411_stage1_DGEList.rds"))
cat("GSE178411 Stage 1 QC complete; gate_status=review_required.\n")
