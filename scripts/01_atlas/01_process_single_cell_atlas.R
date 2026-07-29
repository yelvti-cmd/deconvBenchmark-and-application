#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

source("scripts/utils.R")
args <- parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- load_config(args$config %||% "config/analysis_config.R")
set.seed(cfg$seed)

input <- args$input %||% cfg$paths$single_cell_rds
output <- args$output %||% file.path(cfg$paths$checkpoint_dir, "single_cell_atlas_processed.rds")
require_file(input, "Single-cell Seurat object")
ensure_dir(dirname(output))

obj <- readRDS(input)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 500 & nFeature_RNA < 20000 & percent.mt < 10)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj))

integration_vars <- intersect(
  c(cfg$metadata$source, cfg$metadata$platform, cfg$metadata$sample_id),
  colnames(obj[[]])
)
if (!length(integration_vars)) {
  stop("No configured integration variables are present in the metadata.", call. = FALSE)
}
obj <- RunHarmony(obj, group.by.vars = integration_vars, reduction = "pca")
pcs <- seq_len(min(30L, ncol(Embeddings(obj, "harmony"))))
obj <- FindNeighbors(obj, reduction = "harmony", dims = pcs)
obj <- FindClusters(obj, resolution = 0.2)
obj <- RunUMAP(obj, reduction = "harmony", dims = pcs)
saveRDS(obj, output)
