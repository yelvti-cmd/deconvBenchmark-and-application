#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
source("scripts/utils.R")
args <- parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- load_config(args$config %||% "config/analysis_config.R")
normalization <- args$normalization %||% "none"
reference_name <- args$reference %||% "all_samples"

set.seed(cfg$seed)
generation_file <- file.path(
  cfg$paths$output_dir,
  paste0(reference_name, "_", normalization, "_generation_object.rds")
)
proportion_file <- file.path(
  cfg$paths$output_dir,
  "real_sample_proportions.rds"
)
require_file(generation_file, "Generation-cell object")
require_file(proportion_file, "Real sample-proportion matrix")

object <- readRDS(generation_file)
real_proportions <- readRDS(proportion_file)
cell_type_col <- cfg$metadata$cell_type
object$cell_type_for_simulation <- object@meta.data[[cell_type_col]]
split_objects <- SplitObject(object, split.by = "cell_type_for_simulation")
cell_types <- names(split_objects)
counts_by_type <- lapply(split_objects, get_counts)
genes <- rownames(counts_by_type[[1L]])

draw_proportions <- function(distribution) {
  n <- length(cell_types)
  if (distribution == "uniform") {
    p <- runif(n, 0.01, 0.99)
  } else if (startsWith(distribution, "normal_")) {
    sd_value <- as.numeric(sub("normal_", "", distribution))
    p <- abs(rnorm(n, 0, sd_value))
  } else if (distribution == "bimodal") {
    n_low <- ceiling(n / 2)
    p <- sample(c(
      abs(rnorm(n_low, 0, 0.1)),
      abs(rnorm(n - n_low, 1, 0.1))
    ))
  } else if (distribution == "realistic") {
    aligned <- matrix(0, nrow(real_proportions), n)
    colnames(aligned) <- cell_types
    common <- intersect(colnames(real_proportions), cell_types)
    aligned[, common] <- real_proportions[, common, drop = FALSE]
    p <- aligned[sample(seq_len(nrow(aligned)), 1L), ]
    jitter <- cfg$pseudobulk$realistic_jitter
    p <- p * runif(n, 1 - jitter, 1 + jitter)
  } else {
    stop("Unsupported distribution: ", distribution)
  }
  p <- pmax(as.numeric(p), 1e-9)
  setNames(p / sum(p), cell_types)
}

simulate_one <- function(proportions) {
  target <- round(cfg$pseudobulk$cells_per_sample * proportions)
  target[target < 1L] <- 1L
  target[which.max(target)] <- target[which.max(target)] +
    cfg$pseudobulk$cells_per_sample - sum(target)

  mixture <- numeric(length(genes))
  truth <- matrix(0, length(genes), length(cell_types))
  rownames(truth) <- genes
  colnames(truth) <- cell_types

  for (cell_type in cell_types) {
    matrix_i <- counts_by_type[[cell_type]]
    selected <- sample(
      seq_len(ncol(matrix_i)),
      target[[cell_type]],
      replace = target[[cell_type]] > ncol(matrix_i)
    )
    sampled <- matrix_i[, selected, drop = FALSE]
    mixture <- mixture + Matrix::rowSums(sampled)
    truth[, cell_type] <- Matrix::rowMeans(sampled)
  }
  list(
    mixture = mixture,
    proportions = target / sum(target),
    expression = truth
  )
}

distributions <- cfg$pseudobulk$distributions
n_per_distribution <- cfg$pseudobulk$samples_per_distribution
n_total <- length(distributions) * n_per_distribution
mixtures <- matrix(0, length(genes), n_total, dimnames = list(genes, NULL))
proportions <- matrix(0, length(cell_types), n_total, dimnames = list(cell_types, NULL))
expression <- array(
  0,
  dim = c(length(genes), length(cell_types), n_total),
  dimnames = list(genes, cell_types, NULL)
)

index <- 1L
sample_names <- character(n_total)
for (distribution in distributions) {
  for (replicate_id in seq_len(n_per_distribution)) {
    result <- simulate_one(draw_proportions(distribution))
    sample_name <- paste(distribution, replicate_id, sep = "_")
    mixtures[, index] <- result$mixture
    proportions[, index] <- result$proportions
    expression[, , index] <- result$expression
    sample_names[[index]] <- sample_name
    index <- index + 1L
  }
}

colnames(mixtures) <- sample_names
colnames(proportions) <- sample_names
dimnames(expression)[[3L]] <- sample_names
prefix <- file.path(
  cfg$paths$output_dir,
  paste(reference_name, normalization, sep = "_")
)
saveRDS(mixtures, paste0(prefix, "_pseudobulk.rds"))
saveRDS(proportions, paste0(prefix, "_true_proportions.rds"))
saveRDS(expression, paste0(prefix, "_true_expression.rds"))
