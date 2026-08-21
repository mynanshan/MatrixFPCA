#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]))
} else {
  normalizePath(file.path("application", "adni-latent", "compare_brainiac_layers.R"))
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
setwd(project_root)

defaults <- c(
  block03 = file.path(
    "application", "adni-latent", "results", "brainiac_patch_full_layer03"
  ),
  block06 = file.path(
    "application", "adni-latent", "results", "brainiac_patch_full_layer06"
  ),
  final = file.path(
    "application", "adni-latent", "results", "brainiac_patch_full"
  )
)
arguments <- commandArgs(trailingOnly = TRUE)
output_dir <- file.path(
  "application", "adni-latent", "results", "brainiac_layer_comparison_full"
)
for (argument in arguments) {
  if (grepl("^--output-dir=", argument)) {
    output_dir <- sub("^--output-dir=", "", argument)
  } else {
    stop("Only --output-dir=PATH is supported.")
  }
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  "effective_ranks.csv", "sampled_separability_diagnostic.csv",
  "latent_distance_summary.csv",
  "marginal_eigenspectra.csv", "run_summary.json"
)
missing <- unlist(lapply(defaults, function(directory) {
  paths <- file.path(directory, required_files)
  paths[!file.exists(paths)]
}), use.names = FALSE)
if (length(missing)) {
  stop("Missing layer-analysis outputs: ", paste(missing, collapse = ", "))
}
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required.")

layer_depth <- c(block03 = 3L, block06 = 6L, final = 12L)
layer_label <- c(block03 = "Block 3", block06 = "Block 6", final = "Final normalized")

rows <- lapply(names(defaults), function(layer) {
  directory <- defaults[[layer]]
  effective <- fread(file.path(directory, "effective_ranks.csv"))
  separability <- fread(file.path(directory, "sampled_separability_diagnostic.csv"))
  distances <- fread(file.path(directory, "latent_distance_summary.csv"))
  spectra <- fread(file.path(directory, "marginal_eigenspectra.csv"))
  summary <- jsonlite::read_json(file.path(directory, "run_summary.json"), simplifyVector = TRUE)
  within <- distances[
    pair_type == "within_subject" & representation == "patch_matrix"
  ]
  between <- distances[
    pair_type == "between_subject" & representation == "patch_matrix"
  ]
  token_rank <- effective[mode == "token"]$retained_rank
  feature_rank <- effective[mode == "compressed_feature"]$retained_rank
  time_rank <- effective[mode == "time"]$retained_rank
  data.table(
    layer = layer,
    layer_label = layer_label[[layer]],
    layer_depth = layer_depth[[layer]],
    subjects = summary$subjects,
    visits = summary$visits,
    feature_pca_retained_pve = summary$feature_pca_retained_pve,
    token_effective_rank = effective[mode == "token"]$effective_rank,
    token_retained_rank = token_rank,
    token_pve_at_retained_rank = spectra[
      mode == "token" & component == token_rank
    ]$cumulative_pve,
    feature_effective_rank = effective[mode == "compressed_feature"]$effective_rank,
    feature_retained_rank = feature_rank,
    feature_pve_at_retained_rank = spectra[
      mode == "feature" & component == feature_rank
    ]$cumulative_pve,
    time_effective_rank = effective[mode == "time"]$effective_rank,
    time_retained_rank = time_rank,
    retained_product_components = token_rank * feature_rank * time_rank,
    separability_offdiag_correlation =
      separability$empirical_model_offdiag_correlation,
    separability_relative_frobenius_error = separability$relative_frobenius_error,
    within_patch_median_distance = within$median_distance,
    between_patch_median_distance = between$median_distance,
    within_between_median_ratio = within$median_distance / between$median_distance,
    lag_distance_spearman = within$lag_distance_spearman
  )
})
comparison <- rbindlist(rows)[order(layer_depth)]
fwrite(comparison, file.path(output_dir, "brainiac_layer_comparison.csv"))

png(file.path(output_dir, "brainiac_layer_comparison.png"),
    width = 1400, height = 1200, res = 150)
par(mfrow = c(2L, 2L), mar = c(6.5, 4.6, 3.0, 1.2))
labels <- sub(" normalized$", "", comparison$layer_label)
token_positions <- barplot(
  comparison$token_effective_rank, names.arg = labels, las = 2,
  col = "#4393C3", ylab = "Effective rank", main = "Spatial-token covariance",
  ylim = c(0, max(comparison$token_effective_rank) * 1.16)
)
text(
  token_positions, comparison$token_effective_rank,
  labels = paste0("retained r=", comparison$token_retained_rank),
  pos = 3, cex = 0.78
)
barplot(
  log10(comparison$retained_product_components),
  names.arg = labels, las = 2, col = "#1B7837",
  ylab = "log10 retained product modes", main = "MatrixFPCA complexity",
  ylim = c(0, max(log10(comparison$retained_product_components)) * 1.16)
)
barplot(
  rbind(comparison$separability_offdiag_correlation,
        comparison$separability_relative_frobenius_error),
  beside = TRUE, names.arg = labels, las = 2,
  col = c("#762A83", "#F1B6DA"), ylim = c(0, 1.25),
  ylab = "Metric", main = "Sampled Kronecker diagnostic"
)
legend(
  "topright", legend = c("covariance correlation", "relative error"),
  fill = c("#762A83", "#F1B6DA"), bty = "n", cex = 0.7
)
barplot(
  comparison$within_between_median_ratio, names.arg = labels, las = 2,
  col = "#9970AB", ylab = "Within / between median distance",
  main = "Subject specificity (lower is better)"
)
dev.off()

cat("Wrote ", normalizePath(output_dir), "\n", sep = "")
