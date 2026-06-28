#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1L) args[[1L]] else
  file.path("experiment", "results")
output_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path("experiment", "summary")
expected_repetitions <- if (length(args) >= 3L) as.integer(args[[3L]]) else NA_integer_
manifest_dirs <- c(
  file.path(results_dir, "manifests"),
  file.path("logs", "manifests")
)
dimension_grid <- c(20, 30, 50, 100, 150, 200)

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
  "sseU", "sseV", "iseK", "relsseU", "relsseV", "reliseK", "relsiseC"
)
missing <- setdiff(required, names(res))
if (length(missing)) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

numeric_columns <- c(
  "seed", "n", "p", "q", "sig", "dependence_param", "s1", "s2", "nt",
  "sseU", "sseV", "iseK", "relsseU", "relsseV", "reliseK", "relsiseC"
)
res[, (numeric_columns) := lapply(.SD, as.numeric), .SDcols = numeric_columns]
# Read results created before the method was renamed without requiring a rerun.
res[weight == "optimal", weight := "gaussian_refined"]
res <- res[weight %in% c("identity", "gaussian_refined")]

read_expected_from_manifests <- function(manifest_dirs) {
  manifests <- unlist(lapply(
    manifest_dirs,
    list.files,
    pattern = "^simu_simple_[0-9]{8}_[0-9]{6}(_[A-Za-z0-9_-]+)?\\.tsv$",
    full.names = TRUE
  ))
  if (!length(manifests)) return(NULL)

  manifest_info <- data.table(
    path = manifests,
    timestamp = sub(
      "^simu_simple_([0-9]{8}_[0-9]{6})(_[A-Za-z0-9_-]+)?\\.tsv$",
      "\\1",
      basename(manifests)
    )
  )
  latest_timestamp <- max(manifest_info$timestamp)
  latest_manifests <- manifest_info[timestamp == latest_timestamp, path]
  expected <- rbindlist(lapply(latest_manifests, fread), use.names = TRUE)
  setnames(
    expected,
    c(
      "replicate", "sigma", "distribution", "dep_param"
    ),
    c(
      "seed", "sig", "distr", "dependence_param"
    ),
    skip_absent = TRUE
  )
  expected[, q := p]
  expected[, `:=`(s1 = 2, s2 = 4, nt = 51)]
  expected[
    ,
    .(
      config_id, seed, n, p, q, sig, distr, dependence, dependence_param,
      s1, s2, nt
    )
  ]
}

make_expected_configs <- function(repetitions) {
  configs <- data.table(
    n = numeric(),
    p = numeric(),
    q = numeric(),
    sig = numeric(),
    distr = character(),
    dependence = character(),
    dependence_param = numeric(),
    s1 = numeric(),
    s2 = numeric(),
    nt = numeric()
  )

  add_config <- function(n, p, sig, distr, dependence, dependence_param) {
    candidate <- data.table(
      n = n, p = p, q = p, sig = sig, distr = distr,
      dependence = dependence, dependence_param = dependence_param,
      s1 = 2, s2 = 4, nt = 51
    )
    configs <<- unique(rbind(configs, candidate, use.names = TRUE))
  }

  for (n in c(25, 50, 100)) add_config(n, 30, 0.5, "normal", "independent", 0)
  distribution_grid <- c("normal", "uniform", "beta32", "t5")
  for (p in dimension_grid) {
    for (distr in distribution_grid) {
      add_config(50, p, 0.5, distr, "independent", 0)
    }
  }
  for (sig in c(0, 0.5, 1)) add_config(50, 30, sig, "normal", "independent", 0)
  for (p in dimension_grid) {
    add_config(50, p, 0.5, "normal", "common_amplitude", 1)
    add_config(50, p, 0.5, "normal", "common_signed", 0)
    add_config(50, p, 0.5, "normal", "heavy_amplitude", 1.5)
  }

  configs[, config_id := sprintf("%03d", .I)]
  configs[
    CJ(config_id = config_id, seed = seq_len(repetitions)),
    on = "config_id",
    allow.cartesian = TRUE
  ][
    ,
    .(
      config_id, seed, n, p, q, sig, distr, dependence, dependence_param,
      s1, s2, nt
    )
  ]
}

run_columns <- c(
  "seed", "n", "p", "q", "sig", "distr", "dependence", "dependence_param",
  "s1", "s2", "nt"
)
setting_columns_no_weight <- setdiff(run_columns, "seed")

expected_runs <- read_expected_from_manifests(manifest_dirs)
if (is.null(expected_runs)) {
  if (is.na(expected_repetitions)) {
    expected_repetitions <- max(res$seed, na.rm = TRUE)
  }
  if (!is.finite(expected_repetitions) || expected_repetitions < 1L) {
    stop("Could not determine expected repetitions.")
  }
  expected_runs <- make_expected_configs(expected_repetitions)
} else {
  expected_repetitions <- max(expected_runs$seed, na.rm = TRUE)
}
actual_runs <- unique(res[, c(run_columns, "source_file"), with = FALSE])
actual_run_keys <- unique(actual_runs[, ..run_columns])
missing_runs <- expected_runs[!actual_run_keys, on = run_columns]

expected_method_rows <- expected_runs[, {
  row_identity <- copy(.SD)
  row_identity[, weight := "identity"]
  row_refined <- copy(.SD)
  row_refined[, weight := "gaussian_refined"]
  rbind(row_identity, row_refined)
}, by = config_id]
expected_method_rows[, config_id := NULL]
actual_method_rows <- unique(res[, c(run_columns, "weight"), with = FALSE])
missing_method_rows <- expected_method_rows[
  !actual_method_rows,
  on = c(run_columns, "weight")
]

completion_runs <- copy(expected_runs)
completion_runs[, present := FALSE]
completion_runs[actual_run_keys, present := TRUE, on = run_columns]
completion_by_setting <- completion_runs[
  ,
  .(
    expected_runs = .N,
    completed_runs = sum(present),
    missing_runs = sum(!present),
    missing_seeds = paste(seed[!present], collapse = ",")
  ),
  by = c("config_id", setting_columns_no_weight)
]
completion_by_setting[, completion_rate := completed_runs / expected_runs]
setorder(completion_by_setting, config_id)

# Stored errors are logarithms; convert them back before averaging.
metric_columns <- c("relsseU", "relsseV", "reliseK", "relsiseC")
absolute_metric_columns <- c("sseU", "sseV", "iseK")
for (column in c(metric_columns, absolute_metric_columns)) {
  res[, paste0(column, "_raw") := exp(get(column))]
}
res[, sseU_pixel_raw := sseU_raw / p^2]
res[, sseV_pixel_raw := sseV_raw / q^2]

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
    sseU_pixel_mean = mean(sseU_pixel_raw),
    sseU_pixel_sd = sd(sseU_pixel_raw),
    relsseV_mean = mean(relsseV_raw),
    relsseV_sd = sd(relsseV_raw),
    sseV_pixel_mean = mean(sseV_pixel_raw),
    sseV_pixel_sd = sd(sseV_pixel_raw),
    reliseK_mean = mean(reliseK_raw),
    reliseK_sd = sd(reliseK_raw),
    iseK_mean = mean(iseK_raw),
    iseK_sd = sd(iseK_raw),
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
    "relsseU_mean", "relsseU_sd", "sseU_pixel_mean", "sseU_pixel_sd",
    "relsseV_mean", "relsseV_sd", "sseV_pixel_mean", "sseV_pixel_sd",
    "reliseK_mean", "reliseK_sd", "iseK_mean", "iseK_sd",
    "relsiseC_mean", "relsiseC_sd"
  )
)

# Paired replicate-level comparison. Ratios below one favor Gaussian-refined
# weighting.
paired <- dcast(
  res,
  seed + n + p + q + sig + distr + dependence + dependence_param +
    s1 + s2 + nt ~ weight,
  value.var = paste0(metric_columns, "_raw")
)
for (metric in metric_columns) {
  paired[
    ,
    paste0(metric, "_gaussian_refined_ratio") :=
      get(paste0(metric, "_raw_gaussian_refined")) /
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
plot_data[
  ,
  `:=`(
    weight_label = factor(
      weight,
      levels = c("identity", "gaussian_refined"),
      labels = c("Identity", "Gaussian-refined")
    ),
    distr_label = factor(
      distr,
      levels = c("normal", "uniform", "beta32", "t5"),
      labels = c("Gaussian", "Uniform", "Beta(3,2)", "t(5)")
    )
  )
]

theme_paper <- theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

dimension_metric_data <- melt(
  summary_long,
  id.vars = c(setting_columns, "n_rep"),
  measure.vars = list(
    mean_error = c(
      "sseU_pixel_mean", "sseV_pixel_mean", "reliseK_mean", "relsiseC_mean"
    ),
    sd_error = c(
      "sseU_pixel_sd", "sseV_pixel_sd", "reliseK_sd", "relsiseC_sd"
    )
  ),
  variable.name = "metric_index"
)
dimension_metric_data[
  ,
  metric := factor(
    metric_index,
    levels = 1:4,
    labels = c(
      "U pixel MSE", "V pixel MSE", "K relative ISE", "Full covariance"
    )
  )
]
dimension_metric_data[
  ,
  `:=`(
    weight_label = factor(
      weight,
      levels = c("identity", "gaussian_refined"),
      labels = c("Identity", "Gaussian-refined")
    ),
    distr_label = factor(
      distr,
      levels = c("normal", "uniform", "beta32", "t5"),
      labels = c("Gaussian", "Uniform", "Beta(3,2)", "t(5)")
    ),
    mcse = sd_error / sqrt(n_rep)
  )
]
dimension_metric_data[
  ,
  `:=`(
    error_lower = mean_error * exp(-1.96 * mcse / mean_error),
    error_upper = mean_error * exp(1.96 * mcse / mean_error)
  )
]

dimension_data <- dimension_metric_data[
  n == 50 & sig == 0.5 & dependence == "independent" &
    p %in% dimension_grid
]
dimension_plot <- ggplot(
  dimension_data,
  aes(p, mean_error, color = weight_label, shape = weight_label)
) +
  geom_ribbon(
    aes(
      ymin = error_lower, ymax = error_upper,
      fill = weight_label, group = weight_label
    ),
    alpha = 0.12,
    color = NA
  ) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_grid(metric ~ distr_label, scales = "free_y") +
  labs(
    x = "Spatial dimension (p = q)", y = "Mean error",
    color = "Weighting", shape = "Weighting", fill = "Weighting"
  ) +
  theme_paper

dependent_dimension_data <- dimension_metric_data[
  n == 50 & sig == 0.5 & distr == "normal" &
    dependence != "independent" & p %in% dimension_grid
]
dependent_dimension_plot <- ggplot(
  dependent_dimension_data,
  aes(p, mean_error, color = weight_label, shape = weight_label)
) +
  geom_ribbon(
    aes(
      ymin = error_lower, ymax = error_upper,
      fill = weight_label, group = weight_label
    ),
    alpha = 0.12,
    color = NA
  ) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_grid(metric ~ dependence, scales = "free_y") +
  labs(
    x = "Spatial dimension (p = q)", y = "Mean error",
    color = "Weighting", shape = "Weighting", fill = "Weighting"
  ) +
  theme_paper

sample_size_data <- plot_data[
  p == 30 & sig == 0.5 & distr == "normal" &
    dependence == "independent"
]
sample_size_plot <- ggplot(
  sample_size_data,
  aes(n, mean_error, color = weight_label, shape = weight_label)
) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~metric, scales = "free_y") +
  labs(
    x = "Number of subjects", y = "Mean relative error",
    color = "Weighting", shape = "Weighting"
  ) +
  theme_paper

noise_data <- plot_data[
  n == 50 & p == 30 & distr == "normal" &
    dependence == "independent"
]
noise_plot <- ggplot(
  noise_data,
  aes(sig, mean_error, color = weight_label, shape = weight_label)
) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~metric, scales = "free_y") +
  labs(
    x = "Noise standard deviation", y = "Mean relative error",
    color = "Weighting", shape = "Weighting"
  ) +
  theme_paper

distribution_data <- plot_data[
  n == 50 & p == 30 & sig == 0.5 & dependence == "independent"
]
distribution_plot <- ggplot(
  distribution_data,
  aes(
    distr_label, mean_error, color = weight_label, shape = weight_label,
    group = weight_label
  )
) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~metric, scales = "free_y") +
  labs(
    x = "FPC-score distribution", y = "Mean relative error",
    color = "Weighting", shape = "Weighting"
  ) +
  theme_paper

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(res, file.path(output_dir, "simulation_raw_combined.csv"))
fwrite(summary_long, file.path(output_dir, "simulation_summary_long.csv"))
fwrite(summary_wide, file.path(output_dir, "simulation_summary_wide.csv"))
fwrite(paired, file.path(output_dir, "simulation_paired_ratios.csv"))
fwrite(expected_runs, file.path(output_dir, "simulation_expected_runs.csv"))
fwrite(completion_by_setting, file.path(output_dir, "simulation_completion_by_setting.csv"))
fwrite(missing_runs, file.path(output_dir, "simulation_missing_runs.csv"))
fwrite(missing_method_rows, file.path(output_dir, "simulation_missing_method_rows.csv"))

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
  dependent_dimension_data, dependent_dimension_plot,
  "error_vs_dimension_dependent.pdf", 11, 8
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
cat(
  "Expected", nrow(expected_runs), "seed/config runs across",
  uniqueN(expected_runs$config_id), "settings and", expected_repetitions,
  "seeds.\n"
)
cat(
  "Detected", nrow(missing_runs), "missing seed/config runs and",
  nrow(missing_method_rows), "missing weight-level rows.\n"
)
cat("Summaries and figures written to", output_dir, "\n")
