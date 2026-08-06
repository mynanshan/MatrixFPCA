#!/usr/bin/env Rscript

library(data.table)

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
  "sseGammaU", "sseGammaV", "relsseGammaU", "relsseGammaV",
  "mseX", "relmseX"
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
metric_columns <- c(
  "relsseU", "relsseV", "reliseK", "reliseKc", "relsiseC", "relmseX"
)
diagnostic_metric_columns <- c("relsseGammaU", "relsseGammaV")
absolute_metric_columns <- c(
  "sseU", "sseV", "iseK", "iseKc", "sseGammaU", "sseGammaV", "mseX"
)
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
    relmseX_mean = mean_metric(relmseX_raw),
    relmseX_sd = sd_metric(relmseX_raw),
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
    "relmseX_mean", "relmseX_sd",
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

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(res, file.path(output_dir, "simulation_raw_combined.csv"))
fwrite(summary_long, file.path(output_dir, "simulation_summary_long.csv"))
fwrite(summary_wide, file.path(output_dir, "simulation_summary_wide.csv"))
fwrite(paired, file.path(output_dir, "simulation_paired_ratios.csv"))
fwrite(expected_runs, file.path(output_dir, "simulation_expected_runs.csv"))
fwrite(completion_by_setting, file.path(output_dir, "simulation_completion_by_setting.csv"))
fwrite(missing_runs, file.path(output_dir, "simulation_missing_runs.csv"))
fwrite(missing_method_rows, file.path(output_dir, "simulation_missing_method_rows.csv"))
n_method_rows <- nrow(res)
metric_columns <- c("relsseU", "relsseV", "reliseK", "relsiseC", "relmseX")
raw_metric_columns <- paste0(metric_columns, "_raw")
for (column in metric_columns) {
  raw_column <- paste0(column, "_raw")
  if (!raw_column %in% names(res)) res[, (raw_column) := exp(get(column))]
}
res <- res[weight %in% c("identity", "gaussian_refined")]
res[, nbasis := as.character(nbasis)]
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
res <- res[
  nbasis == plot_nbasis & weight_folds == plot_weight_folds &
    eigensystem == plot_eigensystem
]
res[, weight_label := factor(
  weight,
  levels = c("identity", "gaussian_refined"),
  labels = c("Identity", "Refined")
)]

id_columns <- c(
  "seed", "n", "p", "q", "sig", "distr", "dependence",
  "dependence_param", "s1", "s2", "nt", "nbasis", "weight_folds",
  "rank_row", "rank_col", "rank_time", "eigensystem", "eigvec_set",
  "eigval_set", "weight", "weight_label"
)
long <- melt(
  res,
  id.vars = id_columns,
  measure.vars = raw_metric_columns,
  variable.name = "metric_key",
  value.name = "error"
)
metric_labels <- c(
  relsseU_raw = "U rel. SSE",
  relsseV_raw = "V rel. SSE",
  reliseK_raw = "K rel. ISE",
  relsiseC_raw = "Full cov. rel. SISE",
  # The common time-grid integration weight cancels in this relative error.
  relmseX_raw = "X(t) rel. SISE"
)
covariance_metric_labels <- unname(metric_labels[names(metric_labels) != "relmseX_raw"])
dimension_metric_labels <- covariance_metric_labels
# The preceding five-row dimension plots were 8.4 inches high.  This four-row
# height preserves their 1.18632-inch inner panel height under the same margins.
dimension_plot_height <- 7.01568
long[, metric := factor(
  metric_key,
  levels = names(metric_labels),
  labels = unname(metric_labels)
)]

# Uniform remains in the CSV summaries but is intentionally omitted from plots.
long <- long[distr != "uniform"]
setting_columns <- setdiff(id_columns, c("seed", "weight", "weight_label"))
summary_data <- long[
  ,
  .(
    median_error = median(error, na.rm = TRUE),
    lower = quantile(error, probs = 0.025, na.rm = TRUE),
    upper = quantile(error, probs = 0.975, na.rm = TRUE)
  ),
  by = c(setting_columns, "metric", "weight_label")
]

difference_id <- c("seed", setting_columns, "metric")
difference_formula <- as.formula(paste(
  paste(difference_id, collapse = " + "), "~ weight"
))
difference_data <- dcast(long, difference_formula, value.var = "error")
difference_data[, difference := identity - gaussian_refined]

method_colors <- c(Identity = "#F8766D", Refined = "#00BFC4")
method_fills <- grDevices::adjustcolor(method_colors, alpha.f = 0.16)
names(method_fills) <- names(method_colors)

pretty_log_ticks <- function(limits) {
  powers <- floor(log10(limits[[1L]])):ceiling(log10(limits[[2L]]))
  ticks <- as.vector(outer(c(1, 3), 10^powers))
  ticks[ticks >= limits[[1L]] & ticks <= limits[[2L]]]
}

mean_panel_limits <- function(data) {
  positive <- c(data$lower, data$upper)
  positive <- positive[is.finite(positive) & positive > 0]
  limits <- range(positive)
  log_limits <- log10(limits)
  10^(log_limits + c(-0.06, 0.06) * diff(log_limits))
}

draw_mean_panel <- function(
  data, x_var, show_x = FALSE, show_y = TRUE, ylim = NULL
) {
  x <- sort(unique(data[[x_var]]))
  if (is.null(ylim)) ylim <- mean_panel_limits(data)
  x_limits <- range(x)
  x_padding <- if (length(x) > 1L) 0.5 * min(diff(x)) else 0.5
  x_limits <- x_limits + c(-x_padding, x_padding)
  plot(
    NA, xlim = x_limits, ylim = ylim, log = "y",
    axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i"
  )
  grid(col = "grey90", lty = 1)
  for (method in levels(data$weight_label)) {
    values <- data[weight_label == method][order(get(x_var))]
    polygon(
      c(values[[x_var]], rev(values[[x_var]])),
      c(values$lower, rev(values$upper)),
      border = NA, col = method_fills[[method]]
    )
    lines(
      values[[x_var]], values$median_error,
      col = method_colors[[method]], lwd = 1.6
    )
    points(
      values[[x_var]], values$median_error,
      col = method_colors[[method]], pch = 16, cex = 0.9
    )
  }
  box(col = "grey30")
  if (show_x) {
    axis(
      1, at = x, labels = format(x, trim = TRUE),
      cex.axis = 0.92, gap.axis = -1
    )
  }
  if (show_y) {
    ticks <- pretty_log_ticks(ylim)
    axis(2, at = ticks, labels = format(ticks, scientific = TRUE), las = 1, cex.axis = 0.88)
  }
}

draw_difference_panel <- function(data, x_var, show_x = FALSE, show_y = TRUE) {
  x <- sort(unique(data[[x_var]]))
  groups <- lapply(x, function(value) data[get(x_var) == value, difference])
  boxplot_stats <- boxplot(groups, plot = FALSE)
  ylim <- range(c(0, boxplot_stats$stats), finite = TRUE)
  padding <- max(diff(ylim) * 0.06, .Machine$double.eps)
  ylim <- ylim + c(-padding, padding)
  plot(
    NA, xlim = c(0.5, length(x) + 0.5), ylim = ylim,
    axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i"
  )
  grid(col = "grey90", lty = 1)
  abline(h = 0, col = "#2C7FB8", lwd = 1.3)
  boxplot(
    groups, at = seq_along(x), add = TRUE, axes = FALSE,
    boxwex = 0.55, col = "grey82", border = "grey25",
    medcol = "grey10", medlwd = 1.3,
    whiskcol = "grey35", staplecol = "grey35",
    outline = FALSE
  )
  box(col = "grey30")
  if (show_x) {
    axis(
      1, at = seq_along(x), labels = format(x, trim = TRUE),
      cex.axis = 0.92, gap.axis = -1
    )
  }
  if (show_y) axis(2, las = 1, cex.axis = 0.88)
}

plot_region_device_coordinates <- function() {
  figure_region <- par("fig")
  plot_region <- par("plt")
  outer_region <- par("omd")
  inner_coordinates <- c(
    figure_region[[1L]] + plot_region[[1L]] * diff(figure_region[1L:2L]),
    figure_region[[1L]] + plot_region[[2L]] * diff(figure_region[1L:2L]),
    figure_region[[3L]] + plot_region[[3L]] * diff(figure_region[3L:4L]),
    figure_region[[3L]] + plot_region[[4L]] * diff(figure_region[3L:4L])
  )
  c(
    outer_region[[1L]] + inner_coordinates[1L:2L] *
      diff(outer_region[1L:2L]),
    outer_region[[3L]] + inner_coordinates[3L:4L] *
      diff(outer_region[3L:4L])
  )
}

add_outer_annotations <- function(
  x_label, row_labels = NULL, row_y = NULL,
  x_label_y = 0.018, row_label_x = 0.992
) {
  par(
    fig = c(0, 1, 0, 1), mar = rep(0, 4), oma = rep(0, 4),
    new = TRUE
  )
  plot.new()
  plot.window(
    xlim = c(0, 1), ylim = c(0, 1),
    xaxs = "i", yaxs = "i"
  )
  legend(
    "top", legend = names(method_colors), col = method_colors,
    lty = 1, lwd = 2, pch = 16, pt.cex = 0.9,
    horiz = TRUE, bty = "n", cex = 1.05, inset = 0.03
  )
  text(0.5, x_label_y, x_label, cex = 1.05)
  if (!is.null(row_labels)) {
    if (is.null(row_y)) {
      panel_top <- 0.88
      panel_bottom <- 0.09
      row_y <- panel_top -
        (seq_along(row_labels) - 0.5) * (panel_top - panel_bottom) /
          length(row_labels)
    }
    text(
      row_label_x, row_y, row_labels,
      srt = -90, adj = 0.5, cex = 1.02
    )
  }
}

equal_plot_row_heights <- function(
  device_height, outer_margins, bottom_margins, top_margins
) {
  line_height <- par("csi") * par("mex")
  inner_height <- device_height -
    sum(outer_margins[c(1L, 3L)]) * line_height
  margin_heights <- (bottom_margins + top_margins) * line_height
  plot_height <- (inner_height - sum(margin_heights)) /
    length(margin_heights)
  if (plot_height <= 0) {
    stop("The plotting device is too short for the requested panel margins.")
  }
  plot_height + margin_heights
}

plot_condition_grid <- function(
  means, differences, condition_var, condition_values, condition_labels,
  x_var, x_label, filename, metrics = dimension_metric_labels,
  width = 11.5, height = dimension_plot_height
) {
  nrow <- length(metrics)
  ncol <- 2L * length(condition_values)
  outer_margins <- c(3.2, 2.8, 4.2, 4.4)
  bottom_margins <- c(rep(0.9, nrow - 1L), 2.8)
  top_margins <- c(2.5, rep(0.6, nrow - 1L))
  pdf(file.path(output_dir, filename), width = width, height = height)
  on.exit(dev.off())
  par(oma = outer_margins)
  panel_layout <- matrix(seq_len(nrow * ncol), nrow = nrow, byrow = TRUE)
  # Establish the automatic layout text scaling before converting margin
  # lines to inches; dense layouts reduce csi and therefore margin heights.
  layout(panel_layout)
  row_heights <- equal_plot_row_heights(
    height, outer_margins, bottom_margins, top_margins
  )
  layout(panel_layout, heights = row_heights)
  row_centers <- numeric(nrow)
  chart_bottom <- NA_real_
  chart_right <- NA_real_
  for (row in seq_along(metrics)) {
    row_mean_data <- means[
      metric == metrics[[row]] & get(condition_var) %in% condition_values
    ]
    shared_mean_ylim <- mean_panel_limits(row_mean_data)
    for (condition_index in seq_along(condition_values)) {
      condition_value <- condition_values[[condition_index]]
      condition_label <- condition_labels[[condition_index]]
      mean_subset <- means[
        metric == metrics[[row]] & get(condition_var) == condition_value
      ]
      difference_subset <- differences[
        metric == metrics[[row]] & get(condition_var) == condition_value
      ]
      bottom <- row == nrow
      par(mar = c(if (bottom) 2.8 else 0.9, 3.2, if (row == 1L) 2.5 else 0.6, 0.5))
      draw_mean_panel(
        mean_subset, x_var, show_x = bottom, show_y = TRUE,
        ylim = shared_mean_ylim
      )
      if (condition_index == 1L) {
        panel_coordinates <- plot_region_device_coordinates()
        row_centers[[row]] <- mean(panel_coordinates[3L:4L])
        if (bottom) chart_bottom <- panel_coordinates[[3L]]
      }
      if (row == 1L) title(main = paste(condition_label, "Median error", sep = "\n"), cex.main = 0.96)

      par(mar = c(if (bottom) 2.8 else 0.9, 3.2, if (row == 1L) 2.5 else 0.6, 0.5))
      draw_difference_panel(difference_subset, x_var, show_x = bottom, show_y = TRUE)
      if (condition_index == length(condition_values)) {
        chart_right <- plot_region_device_coordinates()[[2L]]
      }
      if (row == 1L) title(main = paste(condition_label, "Identity - Refined", sep = "\n"), cex.main = 0.96)
    }
  }
  line_height <- par("csi") * par("mex")
  add_outer_annotations(
    x_label, metrics, row_y = row_centers,
    x_label_y = chart_bottom - 3.1 * line_height / par("din")[[2L]],
    row_label_x = chart_right + line_height / par("din")[[1L]]
  )
}

plot_metric_grid <- function(
  means, differences, x_var, x_label, filename, width = 8.6, height = 6.4,
  categorical_x = FALSE, metrics = covariance_metric_labels
) {
  nrow <- 2L
  outer_margins <- c(3.2, 1.8, 4.2, 1.2)
  bottom_margins <- c(0.9, 2.8)
  top_margins <- rep(2.6, nrow)
  pdf(file.path(output_dir, filename), width = width, height = height)
  on.exit(dev.off())
  par(oma = outer_margins)
  panel_layout <- matrix(
    seq_len(2L * length(metrics)), nrow = nrow, byrow = TRUE
  )
  layout(panel_layout)
  row_heights <- equal_plot_row_heights(
    height, outer_margins, bottom_margins, top_margins
  )
  layout(panel_layout, heights = row_heights)
  for (index in seq_along(metrics)) {
    row <- if (index <= 2L) 1L else 2L
    bottom <- row == 2L
    metric_name <- metrics[[index]]
    mean_subset <- means[metric == metric_name]
    difference_subset <- differences[metric == metric_name]
    if (categorical_x) {
      levels_x <- if (is.factor(mean_subset[[x_var]])) {
        levels(droplevels(mean_subset[[x_var]]))
      } else {
        unique(as.character(mean_subset[[x_var]]))
      }
      mean_subset[, plot_x := match(as.character(get(x_var)), levels_x)]
      difference_subset[, plot_x := match(as.character(get(x_var)), levels_x)]
      panel_x_var <- "plot_x"
    } else {
      levels_x <- NULL
      panel_x_var <- x_var
    }

    par(mar = c(if (bottom) 2.8 else 0.9, 3.2, 2.6, 0.5))
    draw_mean_panel(
      mean_subset, panel_x_var,
      show_x = bottom && !categorical_x, show_y = TRUE
    )
    title(main = paste(metric_name, "Median error", sep = "\n"), cex.main = 0.96)
    if (bottom && categorical_x) {
      axis(
        1, at = seq_along(levels_x), labels = levels_x,
        cex.axis = 0.9, gap.axis = -1
      )
    }

    par(mar = c(if (bottom) 2.8 else 0.9, 3.2, 2.6, 0.5))
    draw_difference_panel(
      difference_subset, panel_x_var,
      show_x = bottom && !categorical_x, show_y = TRUE
    )
    title(main = "Identity - Refined", cex.main = 0.96)
    if (bottom && categorical_x) {
      axis(
        1, at = seq_along(levels_x), labels = levels_x,
        cex.axis = 0.9, gap.axis = -1
      )
    }
  }
  add_outer_annotations(x_label)
}

dimension_grid <- c(20, 50, 100, 150, 200)
dimension_means <- summary_data[
  n == 50 & sig == 0.1 & dependence == "independent" & p %in% dimension_grid
]
dimension_differences <- difference_data[
  n == 50 & sig == 0.1 & dependence == "independent" & p %in% dimension_grid
]
plot_condition_grid(
  dimension_means, dimension_differences,
  "distr", c("normal", "beta32", "t5"), c("Gaussian", "Beta(3,2)", "t(5)"),
  "p", "Spatial dimension (p = q)",
  "error_vs_dimension_n50_sig0p1_independent.pdf"
)

dependent_means <- summary_data[
  n == 50 & sig == 0.1 & distr == "normal" & dependence != "independent" &
    p %in% dimension_grid
]
dependent_differences <- difference_data[
  n == 50 & sig == 0.1 & distr == "normal" & dependence != "independent" &
    p %in% dimension_grid
]
plot_condition_grid(
  dependent_means, dependent_differences,
  "dependence",
  c("common_amplitude", "common_signed", "heavy_amplitude"),
  c("Common amplitude", "Common signed", "Heavy amplitude"),
  "p", "Spatial dimension (p = q)",
  "error_vs_dimension_dependent_n50_sig0p1_normal.pdf"
)

sample_means <- summary_data[
  p == 30 & sig == 0.1 & distr == "normal" & dependence == "independent"
]
sample_differences <- difference_data[
  p == 30 & sig == 0.1 & distr == "normal" & dependence == "independent"
]
plot_metric_grid(
  sample_means, sample_differences, "n", "Number of subjects",
  "error_vs_sample_size_p30_sig0p1_normal_independent.pdf"
)

noise_means <- summary_data[
  n == 50 & p == 30 & distr == "normal" & dependence == "independent"
]
noise_differences <- difference_data[
  n == 50 & p == 30 & distr == "normal" & dependence == "independent"
]
plot_metric_grid(
  noise_means, noise_differences, "sig", "Noise standard deviation",
  "error_vs_noise_n50_p30_normal_independent.pdf"
)

cat("Base-graphics figures written to", output_dir, "\n")
cat("Read", length(files), "files and", n_method_rows, "method rows.\n")
cat(
  "Expected", nrow(expected_runs), "seed/config runs across",
  uniqueN(expected_runs$config_id), "settings and", expected_repetitions,
  "seeds.\n"
)
cat(
  "Detected", nrow(missing_runs), "missing seed/config runs and",
  nrow(missing_method_rows), "missing weight-level rows.\n"
)
