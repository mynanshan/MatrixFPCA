#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
summary_dir <- if (length(args) >= 1L) args[[1L]] else
  file.path("experiment", "summary")
output_dir <- if (length(args) >= 2L) args[[2L]] else summary_dir
raw_path <- file.path(summary_dir, "simulation_raw_combined.csv")
if (!file.exists(raw_path)) stop("Missing combined simulation results: ", raw_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

res <- fread(raw_path)
metric_columns <- c("relsseU", "relsseV", "reliseK", "relsiseC")
raw_metric_columns <- paste0(metric_columns, "_raw")
for (column in metric_columns) {
  raw_column <- paste0(column, "_raw")
  if (!raw_column %in% names(res)) res[, (raw_column) := exp(get(column))]
}
res <- res[weight %in% c("identity", "gaussian_refined")]
res[, nbasis := as.character(nbasis)]
plot_nbasis <- if ("theory" %in% res$nbasis) "theory" else sort(unique(res$nbasis))[[1L]]
plot_weight_folds <- if (10 %in% res$weight_folds) 10 else sort(unique(res$weight_folds))[[1L]]
plot_eigensystem <- if ("ar1" %in% res$eigensystem) "ar1" else sort(unique(res$eigensystem))[[1L]]
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
  relsseU_raw = "U relative Frobenius",
  relsseV_raw = "V relative Frobenius",
  reliseK_raw = "K relative ISE",
  relsiseC_raw = "Full covariance"
)
long[, metric := factor(
  metric_key,
  levels = names(metric_labels),
  labels = unname(metric_labels)
)]

# Uniform remains in the CSV summaries but is intentionally omitted from plots.
long <- long[distr != "uniform"]
setting_columns <- setdiff(id_columns, c("seed", "weight", "weight_label"))
mean_data <- long[
  ,
  .(
    mean_error = mean(error, na.rm = TRUE),
    se_error = sd(error, na.rm = TRUE) / sqrt(sum(!is.na(error)))
  ),
  by = c(setting_columns, "metric", "weight_label")
]
mean_data[, `:=`(
  lower = mean_error * exp(-1.96 * se_error / mean_error),
  upper = mean_error * exp(1.96 * se_error / mean_error)
)]

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
      values[[x_var]], values$mean_error,
      col = method_colors[[method]], lwd = 1.6
    )
    points(
      values[[x_var]], values$mean_error,
      col = method_colors[[method]], pch = 16, cex = 0.9
    )
  }
  box(col = "grey30")
  if (show_x) axis(1, at = x, labels = format(x, trim = TRUE), cex.axis = 0.92)
  if (show_y) {
    ticks <- pretty_log_ticks(ylim)
    axis(2, at = ticks, labels = format(ticks, scientific = TRUE), las = 1, cex.axis = 0.88)
  }
}

draw_difference_panel <- function(data, x_var, show_x = FALSE, show_y = TRUE) {
  x <- sort(unique(data[[x_var]]))
  groups <- lapply(x, function(value) data[get(x_var) == value, difference])
  values <- unlist(groups, use.names = FALSE)
  ylim <- range(c(0, values), finite = TRUE)
  padding <- max(diff(ylim) * 0.06, .Machine$double.eps)
  ylim <- ylim + c(-padding, padding)
  plot(
    NA, xlim = c(0.5, length(x) + 0.5), ylim = ylim,
    axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i"
  )
  grid(col = "grey90", lty = 1)
  abline(h = 0, col = "grey45", lwd = 0.8)
  boxplot(
    groups, at = seq_along(x), add = TRUE, axes = FALSE,
    boxwex = 0.55, col = "grey82", border = "grey25",
    medcol = "grey10", medlwd = 1.3,
    whiskcol = "grey35", staplecol = "grey35",
    outcol = "grey45", outpch = 16, outcex = 0.5
  )
  box(col = "grey30")
  if (show_x) axis(1, at = seq_along(x), labels = format(x, trim = TRUE), cex.axis = 0.92)
  if (show_y) axis(2, las = 1, cex.axis = 0.88)
}

add_outer_annotations <- function(x_label, row_labels = NULL) {
  par(
    fig = c(0, 1, 0, 1), mar = rep(0, 4), oma = rep(0, 4),
    new = TRUE
  )
  plot.new()
  legend(
    "top", legend = names(method_colors), col = method_colors,
    lty = 1, lwd = 2, pch = 16, pt.cex = 0.9,
    horiz = TRUE, bty = "n", cex = 1.05, inset = 0.012
  )
  text(0.5, 0.018, x_label, cex = 1.05)
  if (!is.null(row_labels)) {
    panel_top <- 0.88
    panel_bottom <- 0.09
    row_y <- panel_top -
      (seq_along(row_labels) - 0.5) * (panel_top - panel_bottom) /
        length(row_labels)
    text(
      0.992, row_y, row_labels,
      srt = -90, adj = 0.5, cex = 1.02
    )
  }
}

plot_condition_grid <- function(
  means, differences, condition_var, condition_values, condition_labels,
  x_var, x_label, filename, width = 11.5, height = 7.2
) {
  metrics <- levels(long$metric)
  nrow <- length(metrics)
  ncol <- 2L * length(condition_values)
  pdf(file.path(output_dir, filename), width = width, height = height)
  on.exit(dev.off())
  layout(matrix(seq_len(nrow * ncol), nrow = nrow, byrow = TRUE))
  par(oma = c(3.2, 3.8, 4.2, 4.4))
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
      if (row == 1L) title(main = paste(condition_label, "Mean error", sep = "\n"), cex.main = 0.96)

      par(mar = c(if (bottom) 2.8 else 0.9, 3.2, if (row == 1L) 2.5 else 0.6, 0.5))
      draw_difference_panel(difference_subset, x_var, show_x = bottom, show_y = TRUE)
      if (row == 1L) title(main = paste(condition_label, "Identity - Refined", sep = "\n"), cex.main = 0.96)
    }
  }
  add_outer_annotations(x_label, metrics)
}

plot_metric_grid <- function(
  means, differences, x_var, x_label, filename, width = 8.6, height = 6.4,
  categorical_x = FALSE
) {
  metrics <- levels(long$metric)
  pdf(file.path(output_dir, filename), width = width, height = height)
  on.exit(dev.off())
  layout(matrix(seq_len(8L), nrow = 2L, byrow = TRUE))
  par(oma = c(3.2, 3.8, 4.2, 1.2))
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
    title(main = paste(metric_name, "Mean error", sep = "\n"), cex.main = 0.96)
    if (bottom && categorical_x) axis(1, at = seq_along(levels_x), labels = levels_x, cex.axis = 0.9)

    par(mar = c(if (bottom) 2.8 else 0.9, 3.2, 2.6, 0.5))
    draw_difference_panel(
      difference_subset, panel_x_var,
      show_x = bottom && !categorical_x, show_y = TRUE
    )
    title(main = "Identity - Refined", cex.main = 0.96)
    if (bottom && categorical_x) axis(1, at = seq_along(levels_x), labels = levels_x, cex.axis = 0.9)
  }
  add_outer_annotations(x_label)
}

dimension_grid <- c(20, 50, 100, 150, 200)
dimension_means <- mean_data[
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

dependent_means <- mean_data[
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

sample_means <- mean_data[
  p == 30 & sig == 0.1 & distr == "normal" & dependence == "independent"
]
sample_differences <- difference_data[
  p == 30 & sig == 0.1 & distr == "normal" & dependence == "independent"
]
plot_metric_grid(
  sample_means, sample_differences, "n", "Number of subjects",
  "error_vs_sample_size_p30_sig0p1_normal_independent.pdf"
)

noise_means <- mean_data[
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
