#!/usr/bin/env Rscript

# Keep result ingestion and tabulation in the established analysis script,
# then replace its ggplot figures with the base-graphics versions.
args <- commandArgs(trailingOnly = TRUE)
legacy_script <- file.path("experiment", "analyze_simu-simple-ggplot.R")
base_plot_script <- file.path("experiment", "plot_simu-simple-base.R")

status <- system2("Rscript", c(legacy_script, args))
if (!identical(status, 0L)) {
  stop("Simulation result aggregation failed with status ", status, ".")
}

output_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path("experiment", "summary")
status <- system2("Rscript", c(base_plot_script, output_dir, output_dir))
if (!identical(status, 0L)) {
  stop("Base-graphics figure generation failed with status ", status, ".")
}
