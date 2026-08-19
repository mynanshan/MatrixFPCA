#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]))
} else {
  normalizePath(file.path("application", "adni", "summarize_continuum_sensitivity.R"))
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
setwd(project_root)

options <- list(
  output_dir = file.path("application", "adni", "results"),
  reference_pcs = 30L
)
for (argument in commandArgs(trailingOnly = TRUE)) {
  if (!grepl("^--[^=]+=", argument)) {
    stop("Options must use --name=value syntax: ", argument)
  }
  pair <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
  key <- gsub("-", "_", pair[[1L]], fixed = TRUE)
  if (!key %in% names(options)) stop("Unknown option: --", pair[[1L]])
  value <- paste(pair[-1L], collapse = "=")
  options[[key]] <- if (is.integer(options[[key]])) as.integer(value) else value
}

output_dir <- normalizePath(options$output_dir)
diagnostic_files <- Sys.glob(file.path(
  output_dir, "latent_coordinate_diagnostics_pcs*.csv"
))
if (length(diagnostic_files) < 2L) {
  stop("At least two suffixed continuum fits are required for sensitivity analysis.")
}

extract_dimension <- function(path) {
  as.integer(sub(".*_pcs([0-9]+)\\.csv$", "\\1", basename(path)))
}
dimensions <- vapply(diagnostic_files, extract_dimension, integer(1L))
ordering <- order(dimensions)
dimensions <- dimensions[ordering]
diagnostic_files <- diagnostic_files[ordering]

reference_dimension <- options$reference_pcs
if (!reference_dimension %in% dimensions) {
  stop("Reference dimension pcs", reference_dimension, " is not available.")
}
read_scores <- function(dimension) {
  path <- file.path(
    output_dir, sprintf("fpc_scores_embeddings_pcs%d.csv", dimension)
  )
  if (!file.exists(path)) stop("Missing score file: ", path)
  fread(path, select = c("PTID", "latent_coordinate"))
}
reference_scores <- read_scores(reference_dimension)
setnames(reference_scores, "latent_coordinate", "reference_coordinate")

sensitivity <- rbindlist(lapply(seq_along(dimensions), function(index) {
  dimension <- dimensions[[index]]
  diagnostic <- fread(diagnostic_files[[index]])
  scores <- read_scores(dimension)
  comparison <- merge(scores, reference_scores, by = "PTID", sort = FALSE)
  data.table(
    score_pcs = dimension,
    subjects = nrow(scores),
    retained_product_variance_fraction =
      diagnostic$input_retained_product_variance_fraction[[1L]],
    principal_curve_variance_explained =
      diagnostic$principal_curve_variance_explained[[1L]],
    linear_axis_variance_explained =
      diagnostic$linear_axis_variance_explained[[1L]],
    nonlinear_improvement = diagnostic$nonlinear_improvement[[1L]],
    coordinate_vs_age_spearman = diagnostic$coordinate_vs_age_spearman[[1L]],
    reference_pcs = reference_dimension,
    matched_reference_subjects = nrow(comparison),
    coordinate_vs_reference_spearman = cor(
      comparison$latent_coordinate, comparison$reference_coordinate,
      method = "spearman"
    ),
    coordinate_vs_reference_pearson = cor(
      comparison$latent_coordinate, comparison$reference_coordinate,
      method = "pearson"
    ),
    mean_absolute_coordinate_difference = mean(abs(
      comparison$latent_coordinate - comparison$reference_coordinate
    )),
    maximum_absolute_coordinate_difference = max(abs(
      comparison$latent_coordinate - comparison$reference_coordinate
    ))
  )
}))

output_path <- file.path(output_dir, "continuum_dimension_sensitivity.csv")
fwrite(sensitivity, output_path)
print(sensitivity)
cat("Wrote: ", output_path, "\n", sep = "")
