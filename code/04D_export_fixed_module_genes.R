#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(msigdbr); library(dplyr); library(readr) })
root <- normalizePath(".")
reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>% filter(!is.na(gene_symbol))
modules <- c(
  ECM_organization = "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  collagen_formation = "REACTOME_COLLAGEN_FORMATION",
  collagen_biosynthesis = "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES",
  neutrophil_degranulation = "REACTOME_NEUTROPHIL_DEGRANULATION",
  innate_immune_system = "REACTOME_INNATE_IMMUNE_SYSTEM",
  cell_cycle_mitotic = "REACTOME_CELL_CYCLE_MITOTIC"
)
out <- bind_rows(lapply(names(modules), function(nm) reactome %>% filter(gs_name == modules[[nm]]) %>% transmute(module = nm, gene_symbol) %>% distinct()))
write_csv(out, file.path(root, "04-tables/stage4_fixed_reactome_module_genes.csv"))
cat("Fixed module genes exported.\n")
