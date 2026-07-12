#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1L) args[[1L]] else
  file.path("experiment", "results")
output_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path("experiment", "summary")
expected_repetitions <- if (length(args) >= 3L) as.integer(args[[3L]]) else NA_integer_
requested_plot_nbasis <- if (length(args) >= 4L) args[[4L]] else NA_character_
requested_plot_weight_folds <- if (length(args) >= 5L) {
  as.integer(args[[5L]])
} else {
  NA_integer_
}
requested_plot_eigensystem <- if (length(args) >= 6L) args[[6L]] else NA_character_
manifest_dirs <- c(
  file.path(results_dir, "manifests"),
  file.path("logs", "manifests")
)
dimension_grid <- c(20, 50, 100, 150, 200)

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
if (!"rank_row" %in% names(res)) res[, rank_row := p]
if (!"rank_col" %in% names(res)) res[, rank_col := q]
if (!"rank_time" %in% names(res)) res[, rank_time := 5]
if (!"eigensystem" %in% names(res)) res[, eigensystem := "legacy"]
if (!"eigvec_set" %in% names(res)) res[, eigvec_set := 2]
if (!"eigval_set" %in% names(res)) res[, eigval_set := 1]
if (!"nbasis" %in% names(res)) res[, nbasis := "theory"]
if (!"weight_folds" %in% names(res)) res[, weight_folds := 1]
res[, nbasis := as.character(nbasis)]
optional_metric_columns <- c(
  "iseKc", "reliseKc",
  "sseGammaU", "sseGammaV", "relsseGammaU", "relsseGammaV"
)
for (column in optional_metric_columns) {
  if (!column %in% names(res)) res[, (column) := NA_real_]
}

required <- c(
  "seed", "weight", "n", "p", "q", "sig", "distr", "dependence",
  "dependence_param", "s1", "s2", "nt", "nbasis", "weight_folds",
  "rank_row", "rank_col",
  "rank_time", "eigensystem", "eigvec_set", "eigval_set",
  "sseU", "sseV", "iseK", "relsseU", "relsseV", "reliseK", "relsiseC"
)
missing <- setdiff(required, names(res))
if (length(missing)) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

numeric_columns <- c(
  "seed", "n", "p", "q", "sig", "dependence_param", "s1", "s2", "nt",
  "weight_folds", "rank_row", "rank_col", "rank_time", "eigvec_set", "eigval_set",
  "sseU", "sseV", "iseK", "relsseU", "relsseV", "reliseK", "relsiseC",
  optional_metric_columns
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
  if (!"eigensystem" %in% names(expected)) expected[, eigensystem := "legacy"]
  if (!"eigvec_set" %in% names(expected)) expected[, eigvec_set := 2]
  if (!"eigval_set" %in% names(expected)) expected[, eigval_set := 1]
  if (!"rank_time" %in% names(expected)) expected[, rank_time := 5]
  if (!"rank_row" %in% names(expected)) expected[, rank_row := p]
  if (!"rank_col" %in% names(expected)) expected[, rank_col := q]
  if (!"nbasis" %in% names(expected)) expected[, nbasis := "theory"]
  if (!"weight_folds" %in% names(expected)) expected[, weight_folds := 1]
  expected[, nbasis := as.character(nbasis)]
  expected[
    ,
    .(
      config_id, seed, n, p, q, sig, distr, dependence, dependence_param,
      s1, s2, nt, nbasis, weight_folds,
      rank_row, rank_col, rank_time, eigensystem,
      eigvec_set, eigval_set
    )
  ]
}

make_expected_configs <- function(
  repetitions,
  default_nbasis,
  default_weight_folds,
  default_eigensystems
) {
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
    nt = numeric(),
    nbasis = character(),
    weight_folds = numeric(),
    rank_row = numeric(),
    rank_col = numeric(),
    rank_time = numeric(),
    eigensystem = character(),
    eigvec_set = numeric(),
    eigval_set = numeric()
  )

  add_config <- function(n, p, sig, distr, dependence, dependence_param) {
    for (row in seq_len(nrow(default_eigensystems))) {
      eig <- default_eigensystems[row]
      candidate <- data.table(
        n = n, p = p, q = p, sig = sig, distr = distr,
        dependence = dependence, dependence_param = dependence_param,
        s1 = 2, s2 = 4, nt = 51,
        nbasis = default_nbasis,
        weight_folds = default_weight_folds,
        rank_row = p, rank_col = p, rank_time = eig$rank_time,
        eigensystem = eig$eigensystem,
        eigvec_set = eig$eigvec_set,
        eigval_set = eig$eigval_set
      )
      configs <<- unique(rbind(configs, candidate, use.names = TRUE))
    }
  }

  for (n in c(25, 50, 75, 100)) add_config(n, 30, 0.1, "normal", "independent", 0)
  distribution_grid <- c("normal", "uniform", "beta32", "t5")
  for (p in dimension_grid) {
    for (distr in distribution_grid) {
      add_config(50, p, 0.1, distr, "independent", 0)
    }
  }
  for (sig in c(0.1, 0.5, 1)) add_config(50, 30, sig, "normal", "independent", 0)
  for (p in dimension_grid) {
    add_config(50, p, 0.1, "normal", "common_amplitude", 1)
    add_config(50, p, 0.1, "normal", "common_signed", 0)
    add_config(50, p, 0.1, "normal", "heavy_amplitude", 1.5)
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
      s1, s2, nt, nbasis, weight_folds,
      rank_row, rank_col, rank_time, eigensystem,
      eigvec_set, eigval_set
    )
  ]
}

run_columns <- c(
  "seed", "n", "p", "q", "sig", "distr", "dependence", "dependence_param",
  "s1", "s2", "nt", "nbasis", "weight_folds",
  "rank_row", "rank_col", "rank_time",
  "eigensystem", "eigvec_set", "eigval_set"
)
setting_columns_no_weight <- setdiff(run_columns, "seed")
available_nbasis <- sort(unique(res$nbasis))
default_expected_nbasis <- if ("theory" %in% available_nbasis) "theory" else
  available_nbasis[[1L]]
available_weight_folds <- sort(unique(res$weight_folds))
default_expected_weight_folds <- if (10 %in% available_weight_folds) 10 else
  available_weight_folds[[1L]]
default_expected_eigensystems <- unique(res[
  ,
  .(eigensystem, eigvec_set, eigval_set, rank_time)
])
if (!nrow(default_expected_eigensystems)) {
  default_expected_eigensystems <- data.table(
    eigensystem = "ar1",
    eigvec_set = 2,
    eigval_set = 1,
    rank_time = 5
  )
}

expected_runs <- read_expected_from_manifests(manifest_dirs)
if (is.null(expected_runs)) {
  if (is.na(expected_repetitions)) {
    expected_repetitions <- max(res$seed, na.rm = TRUE)
  }
  if (!is.finite(expected_repetitions) || expected_repetitions < 1L) {
    stop("Could not determine expected repetitions.")
  }
  expected_runs <- make_expected_configs(
    expected_repetitions,
    default_nbasis = default_expected_nbasis,
    default_weight_folds = default_expected_weight_folds,
    default_eigensystems = default_expected_eigensystems
  )
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
metric_columns <- c("relsseU", "relsseV", "reliseK", "reliseKc", "relsiseC")
diagnostic_metric_columns <- c("relsseGammaU", "relsseGammaV")
absolute_metric_columns <- c("sseU", "sseV", "iseK", "iseKc", "sseGammaU", "sseGammaV")
for (column in c(metric_columns, diagnostic_metric_columns, absolute_metric_columns)) {
  res[, paste0(column, "_raw") := exp(get(column))]
}
res[, sseU_pixel_raw := sseU_raw / p^2]
res[, sseV_pixel_raw := sseV_raw / q^2]

mean_metric <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}
sd_metric <- function(x) {
  if (sum(!is.na(x)) < 2L) return(NA_real_)
  sd(x, na.rm = TRUE)
}

setting_columns <- c(
  "n", "p", "q", "sig", "distr", "dependence", "dependence_param",
  "s1", "s2", "nt", "nbasis", "weight_folds",
  "rank_row", "rank_col", "rank_time",
  "eigensystem", "eigvec_set", "eigval_set", "weight"
)

summary_long <- res[
  ,
  .(
    n_rep = .N,
    relsseU_mean = mean_metric(relsseU_raw),
    relsseU_sd = sd_metric(relsseU_raw),
    sseU_pixel_mean = mean_metric(sseU_pixel_raw),
    sseU_pixel_sd = sd_metric(sseU_pixel_raw),
    relsseV_mean = mean_metric(relsseV_raw),
    relsseV_sd = sd_metric(relsseV_raw),
    sseV_pixel_mean = mean_metric(sseV_pixel_raw),
    sseV_pixel_sd = sd_metric(sseV_pixel_raw),
    reliseK_mean = mean_metric(reliseK_raw),
    reliseK_sd = sd_metric(reliseK_raw),
    iseK_mean = mean_metric(iseK_raw),
    iseK_sd = sd_metric(iseK_raw),
    reliseKc_mean = mean_metric(reliseKc_raw),
    reliseKc_sd = sd_metric(reliseKc_raw),
    iseKc_mean = mean_metric(iseKc_raw),
    iseKc_sd = sd_metric(iseKc_raw),
    relsiseC_mean = mean_metric(relsiseC_raw),
    relsiseC_sd = sd_metric(relsiseC_raw),
    relsseGammaU_mean = mean_metric(relsseGammaU_raw),
    relsseGammaU_sd = sd_metric(relsseGammaU_raw),
    relsseGammaV_mean = mean_metric(relsseGammaV_raw),
    relsseGammaV_sd = sd_metric(relsseGammaV_raw)
  ),
  by = setting_columns
]
setorder(
  summary_long,
  dependence, distr, n, p, sig, dependence_param, weight
)

summary_wide <- dcast(
  summary_long,
  n + p + q + sig + distr + dependence + dependence_param + s1 + s2 + nt +
    nbasis + weight_folds + rank_row + rank_col + rank_time + eigensystem +
    eigvec_set + eigval_set ~
    weight,
  value.var = c(
    "n_rep",
    "relsseU_mean", "relsseU_sd", "sseU_pixel_mean", "sseU_pixel_sd",
    "relsseV_mean", "relsseV_sd", "sseV_pixel_mean", "sseV_pixel_sd",
    "reliseK_mean", "reliseK_sd", "iseK_mean", "iseK_sd",
    "reliseKc_mean", "reliseKc_sd", "iseKc_mean", "iseKc_sd",
    "relsiseC_mean", "relsiseC_sd",
    "relsseGammaU_mean", "relsseGammaU_sd",
    "relsseGammaV_mean", "relsseGammaV_sd"
  )
)

# Paired replicate-level comparison. Ratios below one favor Gaussian-refined
# weighting.
paired <- dcast(
  res,
  seed + n + p + q + sig + distr + dependence + dependence_param +
    s1 + s2 + nt + nbasis + weight_folds +
    rank_row + rank_col + rank_time + eigensystem +
    eigvec_set + eigval_set ~ weight,
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

theme_paper <- theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 9, hjust = 0.5)
  )

plot_nbasis <- requested_plot_nbasis
if (is.na(plot_nbasis)) {
  plot_nbasis <- if ("theory" %in% available_nbasis) "theory" else
    available_nbasis[[1L]]
}
if (!plot_nbasis %in% available_nbasis) {
  warning(
    "Requested nbasis = ", plot_nbasis, " is absent; using ",
    available_nbasis[[1L]], "."
  )
  plot_nbasis <- available_nbasis[[1L]]
}
plot_weight_folds <- requested_plot_weight_folds
if (is.na(plot_weight_folds)) {
  plot_weight_folds <- if (10 %in% available_weight_folds) 10 else
    available_weight_folds[[1L]]
}
if (!plot_weight_folds %in% available_weight_folds) {
  warning(
    "Requested weight_folds = ", plot_weight_folds, " is absent; using ",
    available_weight_folds[[1L]], "."
  )
  plot_weight_folds <- available_weight_folds[[1L]]
}
available_eigensystems <- sort(unique(res$eigensystem))
plot_eigensystem <- requested_plot_eigensystem
if (is.na(plot_eigensystem)) {
  plot_eigensystem <- if ("ar1" %in% available_eigensystems) "ar1" else
    available_eigensystems[[1L]]
}
if (!plot_eigensystem %in% available_eigensystems) {
  warning(
    "Requested eigensystem = ", plot_eigensystem, " is absent; using ",
    available_eigensystems[[1L]], "."
  )
  plot_eigensystem <- available_eigensystems[[1L]]
}
raw_plot_data <- melt(
  res,
  id.vars = c("seed", setting_columns),
  measure.vars = c(
    "relsseU_raw", "relsseV_raw", "reliseK_raw", "relsiseC_raw"
  ),
  variable.name = "metric",
  value.name = "error"
)
raw_plot_data[
  ,
  metric := factor(
    metric,
    levels = c(
      "relsseU_raw", "relsseV_raw", "reliseK_raw", "relsiseC_raw"
    ),
    labels = c(
      "U relative Frobenius", "V relative Frobenius",
      "K relative ISE", "Full covariance"
    )
  )
]
raw_plot_data[
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

# Uniform scores are retained in the saved summaries but omitted from figures:
# their behavior is very similar to the more asymmetric Beta(3,2) case.
raw_plot_data <- raw_plot_data[distr != "uniform"]

mean_plot_data <- raw_plot_data[
  ,
  .(
    mean_error = mean(error, na.rm = TRUE),
    se_error = sd(error, na.rm = TRUE) / sqrt(sum(!is.na(error)))
  ),
  by = c(setdiff(setting_columns, "weight"), "metric", "weight_label", "distr_label")
]
mean_plot_data[
  ,
  `:=`(
    error_lower = pmax(.Machine$double.xmin, mean_error - 1.96 * se_error),
    error_upper = mean_error + 1.96 * se_error
  )
]

difference_id_columns <- c(
  "seed", setdiff(setting_columns, "weight"), "metric", "distr_label"
)
difference_formula <- stats::as.formula(paste(
  paste(difference_id_columns, collapse = " + "), "~ weight"
))
difference_plot_data <- dcast(
  raw_plot_data,
  difference_formula,
  value.var = "error"
)
difference_plot_data[, error_difference := identity - gaussian_refined]

mean_trend_panel <- function(data, x_var, panel_title, x_label, band = TRUE) {
  plot <- ggplot(
    data,
    aes(
      x = .data[[x_var]], y = mean_error,
      color = weight_label, fill = weight_label, group = weight_label
    )
  )
  if (band) {
    plot <- plot + geom_ribbon(
      aes(ymin = error_lower, ymax = error_upper),
      alpha = 0.14, color = NA
    )
  } else {
    plot <- plot + geom_errorbar(
      aes(ymin = error_lower, ymax = error_upper),
      position = position_dodge(width = 0.35), width = 0.12,
      linewidth = 0.45
    )
  }
  plot +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.7) +
    scale_y_log10() +
    facet_wrap(~metric, ncol = 1, scales = "free_y", strip.position = "right") +
    labs(title = panel_title, x = x_label, y = NULL, color = "Weighting", fill = "Weighting") +
    theme_paper
}

difference_box_panel <- function(data, x_var, panel_title, x_label) {
  ggplot(
    data,
    aes(
      x = factor(.data[[x_var]]), y = error_difference,
      group = factor(.data[[x_var]])
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey45") +
    geom_boxplot(
      width = 0.48, fill = "grey82", color = "grey25",
      linewidth = 0.4, outlier.alpha = 0.35, outlier.size = 0.65
    ) +
    facet_wrap(~metric, ncol = 1, scales = "free_y", strip.position = "right") +
    labs(title = panel_title, x = x_label, y = NULL) +
    theme_paper +
    theme(legend.position = "none")
}

condition_pair <- function(
  mean_data, difference_data, condition_var, condition_value,
  condition_label, x_var, x_label, band = TRUE
) {
  mean_subset <- mean_data[get(condition_var) == condition_value]
  difference_subset <- difference_data[get(condition_var) == condition_value]
  mean_panel <- mean_trend_panel(
    mean_subset, x_var, paste0(condition_label, "\nMean error"), x_label, band
  )
  difference_panel <- difference_box_panel(
    difference_subset, x_var,
    paste0(condition_label, "\nIdentity - Gaussian-refined"), x_label
  )
  patchwork::wrap_plots(
    list(mean_panel, difference_panel), nrow = 1, widths = c(1.35, 1)
  )
}

assemble_condition_pairs <- function(
  mean_data, difference_data, condition_var, condition_values,
  condition_labels, x_var, x_label, ncol = length(condition_values),
  band = TRUE
) {
  pairs <- Map(
    function(value, label) condition_pair(
      mean_data, difference_data, condition_var, value, label,
      x_var, x_label, band
    ),
    condition_values, condition_labels
  )
  patchwork::wrap_plots(pairs, ncol = ncol, guides = "collect") &
    theme(legend.position = "top")
}
raw_plot_data[
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

dimension_data <- mean_plot_data[
  nbasis == plot_nbasis & weight_folds == plot_weight_folds &
    eigensystem == plot_eigensystem &
    n == 50 & sig == 0.1 & dependence == "independent" &
    p %in% dimension_grid
]
dimension_difference_data <- difference_plot_data[
  nbasis == plot_nbasis & weight_folds == plot_weight_folds &
    eigensystem == plot_eigensystem &
    n == 50 & sig == 0.1 & dependence == "independent" &
    p %in% dimension_grid
]
dimension_plot <- assemble_condition_pairs(
  dimension_data, dimension_difference_data,
  "distr", c("normal", "beta32", "t5"),
  c("Gaussian", "Beta(3,2)", "t(5)"),
  "p", "Spatial dimension (p = q)"
)

dependent_dimension_data <- mean_plot_data[
  nbasis == plot_nbasis & weight_folds == plot_weight_folds &
    eigensystem == plot_eigensystem &
    n == 50 & sig == 0.1 & distr == "normal" &
    dependence != "independent" & p %in% dimension_grid
]
dependent_dimension_difference_data <- difference_plot_data[
  nbasis == plot_nbasis & weight_folds == plot_weight_folds &
    eigensystem == plot_eigensystem &
    n == 50 & sig == 0.1 & distr == "normal" &
    dependence != "independent" & p %in% dimension_grid
]
dependent_dimension_plot <- assemble_condition_pairs(
  dependent_dimension_data, dependent_dimension_difference_data,
  "dependence",
  c("common_amplitude", "common_signed", "heavy_amplitude"),
  c("Common amplitude", "Common signed", "Heavy amplitude"),
  "p", "Spatial dimension (p = q)"
)

sample_size_data <- mean_plot_data[
  nbasis == plot_nbasis & weight_folds == plot_weight_folds &
    eigensystem == plot_eigensystem &
    p == 30 & sig == 0.1 & distr == "normal" &
    dependence == "independent"
]
sample_size_difference_data <- difference_plot_data[
  nbasis == plot_nbasis & weight_folds == plot_weight_folds &
    eigensystem == plot_eigensystem &
    p == 30 & sig == 0.1 & distr == "normal" & dependence == "independent"
]
metric_values <- levels(raw_plot_data$metric)
sample_size_plot <- assemble_condition_pairs(
  sample_size_data, sample_size_difference_data,
  "metric", metric_values, metric_values,
  "n", "Number of subjects", ncol = 2
)

noise_data <- mean_plot_data[
  nbasis == plot_nbasis & weight_folds == plot_weight_folds &
    eigensystem == plot_eigensystem &
    n == 50 & p == 30 & distr == "normal" &
    dependence == "independent"
]
noise_difference_data <- difference_plot_data[
  nbasis == plot_nbasis & weight_folds == plot_weight_folds &
    eigensystem == plot_eigensystem &
    n == 50 & p == 30 & distr == "normal" & dependence == "independent"
]
noise_plot <- assemble_condition_pairs(
  noise_data, noise_difference_data,
  "metric", metric_values, metric_values,
  "sig", "Noise standard deviation", ncol = 2
)

distribution_data <- mean_plot_data[
  nbasis == plot_nbasis & weight_folds == plot_weight_folds &
    eigensystem == plot_eigensystem &
    n == 50 & p == 50 & sig == 0.1 & dependence == "independent"
]
distribution_difference_data <- difference_plot_data[
  nbasis == plot_nbasis & weight_folds == plot_weight_folds &
    eigensystem == plot_eigensystem &
    n == 50 & p == 50 & sig == 0.1 & dependence == "independent"
]
distribution_plot <- assemble_condition_pairs(
  distribution_data, distribution_difference_data,
  "metric", metric_values, metric_values,
  "distr_label", "FPC-score distribution", ncol = 2, band = FALSE
)

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
  dimension_data, dimension_plot,
  "error_vs_dimension_n50_sig0p1_independent.pdf", 15, 8
)
save_plot_if_data(
  dependent_dimension_data, dependent_dimension_plot,
  "error_vs_dimension_dependent_n50_sig0p1_normal.pdf", 15, 8
)
save_plot_if_data(
  sample_size_data, sample_size_plot,
  "error_vs_sample_size_p30_sig0p1_normal_independent.pdf", 11, 7
)
save_plot_if_data(
  noise_data, noise_plot,
  "error_vs_noise_n50_p30_normal_independent.pdf", 11, 7
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
