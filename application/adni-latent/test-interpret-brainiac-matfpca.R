#!/usr/bin/env Rscript

# Fast source-level regression checks for the interpretation analysis. The full
# Full data dependencies are integration-tested by the reproducible main run.

path <- file.path(
  "application", "adni-latent", "interpret_brainiac_matfpca.R"
)
stopifnot(file.exists(path))
text <- readLines(path, warn = FALSE)
stopifnot(
  any(grepl("image_reconstruction_performed = FALSE", text, fixed = TRUE)),
  any(grepl("score_reliability <- function", text, fixed = TRUE)),
  any(grepl("bootstrap_coordinate <- function", text, fixed = TRUE)),
  any(grepl("cross_layer_score_concordance.csv", text, fixed = TRUE)),
  any(grepl("localized_eigenmode_candidates.csv", text, fixed = TRUE)),
  any(grepl("coordinate_vs_full_abs_spearman", text, fixed = TRUE)),
  any(grepl("observed_mri_continuum_", text, fixed = TRUE)),
  any(grepl("fpc_score_embedding_comparison_by_continuum_", text, fixed = TRUE)),
  any(grepl("raw_vs_brainiac_observed_mri_continua.png", text, fixed = TRUE)),
  any(grepl("minimum_subjects = 250L", text, fixed = TRUE))
)
cat("BrainIAC MatrixFPCA interpretation source checks passed.\n")
