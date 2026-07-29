#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Seurat); library(Scissor) })
source("scripts/utils.R")
args <- parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- load_config(args$config %||% "config/analysis_config.R")
required <- c("single_cell", "bulk", "phenotype")
if (!all(required %in% names(args))) {
  stop("Supply --single_cell, --bulk, and --phenotype.", call. = FALSE)
}
sce <- readRDS(args$single_cell)
bulk <- as.matrix(read.csv(args$bulk, row.names = 1, check.names = FALSE))
phenotype <- read.csv(args$phenotype, row.names = 1, check.names = FALSE)
if (!identical(colnames(bulk), rownames(phenotype))) {
  stop("Bulk sample columns must match phenotype row names in the same order.")
}
family <- args$family %||% "gaussian"
cutoff <- as.numeric(args$cutoff %||% 0.03)
result <- Scissor(bulk, sce, phenotype, family = family, cutoff = cutoff)
out <- args$output %||% file.path(cfg$paths$output_dir, "application", "scissor_result.rds")
ensure_dir(dirname(out))
saveRDS(result, out)
