#!/usr/bin/env Rscript

# Install only software dependencies. Competitor source code remains in this
# directory and result files remain under experiment/results-methods/.
arguments <- commandArgs(trailingOnly = TRUE)
repository <- "https://cloud.r-project.org"
library_path <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
if (!is.na(library_path) && nzchar(library_path)) {
  dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(library_path, .libPaths()))
}

packages <- c(
  "fda", "rTensor", "data.table", "Matrix", "argparse",
  "fdapace", "face", "MFPCA", "mfaces", "mGSFPCA", "fdarep",
  "future", "future.apply", "RhpcBLASctl"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (!length(missing)) {
  cat("All benchmark R dependencies are already installed.\n")
  quit(status = 0L)
}
cat("Installing: ", paste(missing, collapse = ", "), "\n", sep = "")
install.packages(missing, repos = repository)
still_missing <- missing[!vapply(
  missing, requireNamespace, logical(1L), quietly = TRUE
)]
if (length(still_missing)) {
  stop("Installation did not provide: ", paste(still_missing, collapse = ", "))
}
cat("Dependency installation completed.\n")
