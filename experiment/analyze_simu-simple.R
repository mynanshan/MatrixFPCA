#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1L) args[[1L]] else
  file.path("experiment", "results")
output_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path(results_dir, "summary")

files <- list.files(
  results_dir,
  pattern = "^simu_.*\\.csv$",
  full.names = TRUE
)
if (!length(files)) {
  stop("No simulation CSV files found in ", results_dir)
}

read_result <- function(path) {
  out <- tryCatch(
    fread(path),
    error = function(e) {
      warning("Skipping unreadable file ", path, ": ", e$message)
      NULL
    }
  )
  if (!is.null(out)) out[, source_file := basename(path)]
  out
}

res <- rbindlist(lapply(files, read_result), use.names = TRUE, fill = TRUE)

# Compatibility with the preceding result schema.
if (!"dependence_param" %in% names(res)) {
  if ("amplitude_var" %in% names(res)) {
    setnames(res, "amplitude_var", "dependence_param")
  } else {
    res[, dependence_param := 0]
  }
}
if (!"dependence" %in% names(res)) res[, dependence := "independent"]
if (!"distr" %in% names(res)) res[, distr := "unknown"]

required <- c(
  "seed", "weight", "n", "p", "q", "sig", "distr", "dependence",
  "dependence_param", "s1", "s2", "nt",
  "relsseU", "relsseV", "reliseK", "relsiseC"
)
missing <- setdiff(required, names(res))
if (length(missing)) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

numeric_columns <- c(
  "seed", "n", "p", "q", "sig", "dependence_param", "s1", "s2", "nt",
  "relsseU", "relsseV", "reliseK", "relsiseC"
)
res[, (numeric_columns) := lapply(.SD, as.numeric), .SDcols = numeric_columns]
res <- res[weight %in% c("identity", "optimal")]

# Stored errors are logarithms; convert them back before averaging.
metric_columns <- c("relsseU", "relsseV", "reliseK", "relsiseC")
for (column in metric_columns) {
  res[, paste0(column, "_raw") := exp(get(column))]
}

setting_columns <- c(
  "n", "p", "q", "sig", "distr", "dependence", "dependence_param",
  "s1", "s2", "nt", "weight"
)

summary_long <- res[
  ,
  .(
    n_rep = .N,
    relsseU_mean = mean(relsseU_raw),
    relsseU_sd = sd(relsseU_raw),
    relsseV_mean = mean(relsseV_raw),
    relsseV_sd = sd(relsseV_raw),
    reliseK_mean = mean(reliseK_raw),
    reliseK_sd = sd(reliseK_raw),
    relsiseC_mean = mean(relsiseC_raw),
    relsiseC_sd = sd(relsiseC_raw)
  ),
  by = setting_columns
]
setorder(
  summary_long,
  dependence, distr, n, p, sig, dependence_param, weight
)

summary_wide <- dcast(
  summary_long,
  n + p + q + sig + distr + dependence + dependence_param + s1 + s2 + nt ~
    weight,
  value.var = c(
    "n_rep",
    "relsseU_mean", "relsseU_sd",
    "relsseV_mean", "relsseV_sd",
    "reliseK_mean", "reliseK_sd",
    "relsiseC_mean", "relsiseC_sd"
  )
)

# Paired replicate-level comparison. Ratios below one favor optimal weighting.
paired <- dcast(
  res,
  seed + n + p + q + sig + distr + dependence + dependence_param +
    s1 + s2 + nt ~ weight,
  value.var = paste0(metric_columns, "_raw")
)
for (metric in metric_columns) {
  paired[
    ,
    paste0(metric, "_optimal_ratio") :=
      get(paste0(metric, "_raw_optimal")) /
      get(paste0(metric, "_raw_identity"))
  ]
}

plot_data <- melt(
  summary_long,
  id.vars = setting_columns,
  measure.vars = c(
    "relsseU_mean", "relsseV_mean", "reliseK_mean", "relsiseC_mean"
  ),
  variable.name = "metric",
  value.name = "mean_error"
)
plot_data[
  ,
  metric := factor(
    metric,
    levels = c(
      "relsseU_mean", "relsseV_mean", "reliseK_mean", "relsiseC_mean"
    ),
    labels = c("U", "V", "K", "Full covariance")
  )
]

theme_paper <- theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

dimension_data <- plot_data[
  n == 50 & sig == 0.5 & distr == "normal" &
    p %in% c(15, 30, 60)
]
dimension_plot <- ggplot(
  dimension_data,
  aes(p, mean_error, color = weight, shape = weight)
) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_grid(metric ~ dependence, scales = "free_y") +
  labs(x = "Spatial dimension (p = q)", y = "Mean relative error") +
  theme_paper

sample_size_data <- plot_data[
  p == 30 & sig == 0.5 & distr == "normal" &
    dependence == "independent"
]
sample_size_plot <- ggplot(
  sample_size_data,
  aes(n, mean_error, color = weight, shape = weight)
) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~metric, scales = "free_y") +
  labs(x = "Number of subjects", y = "Mean relative error") +
  theme_paper

noise_data <- plot_data[
  n == 50 & p == 30 & distr == "normal" &
    dependence == "independent"
]
noise_plot <- ggplot(
  noise_data,
  aes(sig, mean_error, color = weight, shape = weight)
) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~metric, scales = "free_y") +
  labs(x = "Noise standard deviation", y = "Mean relative error") +
  theme_paper

distribution_data <- plot_data[
  n == 50 & p == 30 & sig == 0.5 & dependence == "independent"
]
distribution_plot <- ggplot(
  distribution_data,
  aes(distr, mean_error, color = weight, shape = weight, group = weight)
) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~metric, scales = "free_y") +
  labs(x = "FPC-score distribution", y = "Mean relative error") +
  theme_paper

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(res, file.path(output_dir, "simulation_raw_combined.csv"))
fwrite(summary_long, file.path(output_dir, "simulation_summary_long.csv"))
fwrite(summary_wide, file.path(output_dir, "simulation_summary_wide.csv"))
fwrite(paired, file.path(output_dir, "simulation_paired_ratios.csv"))

save_plot_if_data <- function(data, plot, filename, width, height) {
  if (!nrow(data)) {
    message("Skipping ", filename, ": required experiment series is absent.")
    return(invisible(FALSE))
  }
  ggsave(
    file.path(output_dir, filename),
    plot,
    width = width,
    height = height
  )
  invisible(TRUE)
}

save_plot_if_data(
  dimension_data, dimension_plot, "error_vs_dimension.pdf", 11, 8
)
save_plot_if_data(
  sample_size_data, sample_size_plot, "error_vs_sample_size.pdf", 8, 5
)
save_plot_if_data(
  noise_data, noise_plot, "error_vs_noise.pdf", 8, 5
)
save_plot_if_data(
  distribution_data, distribution_plot, "error_vs_distribution.pdf", 8, 5
)

cat("Read", length(files), "files and", nrow(res), "method rows.\n")
cat("Summaries and figures written to", output_dir, "\n")
