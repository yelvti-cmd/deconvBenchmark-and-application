parse_cli <- function(args) {
  values <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must be supplied as --key value pairs.", call. = FALSE)
    }
    values[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  values
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) y else x
}

load_config <- function(path) {
  if (!file.exists(path)) {
    stop("Configuration file not found: ", path, call. = FALSE)
  }
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  if (!exists("config", envir = env, inherits = FALSE)) {
    stop("The configuration file must define an object named 'config'.", call. = FALSE)
  }
  env$config
}

require_file <- function(path, description = "Input file") {
  if (!file.exists(path)) {
    stop(description, " not found: ", path, call. = FALSE)
  }
  invisible(path)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

get_counts <- function(object, assay = "RNA") {
  assay_object <- object[[assay]]
  if (inherits(assay_object, "Assay5")) {
    SeuratObject::LayerData(object, assay = assay, layer = "counts")
  } else {
    SeuratObject::GetAssayData(object, assay = assay, slot = "counts")
  }
}
