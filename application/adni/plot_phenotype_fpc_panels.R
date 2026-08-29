#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]))
} else {
  normalizePath(file.path("application", "adni", "plot_phenotype_fpc_panels.R"))
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
setwd(project_root)
source(file.path("R", "matfpca.R"))

result_dir <- file.path("application", "adni", "results")
model <- readRDS(file.path(result_dir, "adni_matfpca_model_pcs30.rds"))
associations <- fread(file.path(
  result_dir, "fpc_association_all_components.csv"
))

trapezoid_weights <- function(grid) {
  c(
    (grid[[2L]] - grid[[1L]]) / 2,
    (grid[3:length(grid)] - grid[1:(length(grid) - 2L)]) / 2,
    (grid[[length(grid)]] - grid[[length(grid) - 1L]]) / 2
  )
}

# Reconstruct the same two-mode integral eigensystem used by the association
# analysis, so plotted T2 rows are the diagnosed non-constant temporal mode.
time_weights <- trapezoid_weights(model$t_grid_years)
sqrt_weights <- sqrt(time_weights)
weighted_covariance <- sweep(
  sweep(model$covariance$K, 1L, sqrt_weights, "*"),
  2L, sqrt_weights, "*"
)
time_eigen <- eigen(
  0.5 * (weighted_covariance + t(weighted_covariance)), symmetric = TRUE
)
time_values <- pmax(time_eigen$values, 0)
time_functions <- sweep(time_eigen$vectors, 1L, sqrt_weights, "/")
full_eigen <- eigen.matfcov(model$covariance)

select_components <- function(effect_column, t1_count = 3L, t2_count = 3L) {
  # T2 scores need enough longitudinal information to be interpretable.
  # A reliability threshold prevents tiny, heavily shrunk components from
  # winning solely because of an unstable sample association.
  t1 <- associations[time_component == 1L][
    order(-get(effect_column))
  ][seq_len(t1_count)]
  t2_candidates <- associations[
    time_component == 2L & mean_score_reliability >= 0.20
  ][order(-get(effect_column))]
  t2 <- t2_candidates[seq_len(min(t2_count, nrow(t2_candidates)))]
  rbindlist(list(
    cbind(selection_block = "stable T1", t1),
    cbind(selection_block = "dynamic T2", t2)
  ), use.names = TRUE)
}

age_components <- select_components(
  "age_partial_r2_given_group_and_schedule"
)
conversion_components <- select_components(
  "mci_converter_partial_r2_given_age_and_schedule"
)
age_components[, `:=`(
  phenomenon = "older age",
  effect_r2 = age_partial_r2_given_group_and_schedule,
  effect_beta = age_schedule_adjusted_standardized_beta,
  effect_fdr = age_schedule_adjusted_fdr_within_time
)]
conversion_components[, `:=`(
  phenomenon = "MCI-to-AD conversion",
  effect_r2 = mci_converter_partial_r2_given_age_and_schedule,
  effect_beta = mci_converter_age_schedule_adjusted_standardized_beta,
  effect_fdr = mci_converter_age_schedule_adjusted_fdr_within_time
)]
selected_components <- rbindlist(list(
  age_components, conversion_components
), use.names = TRUE, fill = TRUE)
selected_components[, selection_rank := seq_len(.N), by = .(
  phenomenon, selection_block
)]
fwrite(
  selected_components[, .(
    phenomenon, selection_block, selection_rank, component,
    row_component, column_component, time_component,
    product_rank_within_full_spatial_x2, theoretical_variance_fraction,
    mean_score_reliability, effect_r2, effect_beta, effect_fdr
  )],
  file.path(result_dir, "phenotype_fpc_panel_selection.csv")
)

target_times <- 0:4
time_indices <- vapply(
  target_times,
  function(value) which.min(abs(model$t_grid_years - value)),
  integer(1L)
)
display_times <- model$t_grid_years[time_indices]
brain_mask <- model$brain_mask
diverging_palette <- colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(255L)

format_q <- function(value) {
  if (!is.finite(value)) return("NA")
  if (value < 0.001) format(value, scientific = TRUE, digits = 2) else
    sprintf("%.3f", value)
}

draw_label_panel <- function(row, phenomenon_label) {
  plot.new()
  text(
    0.98, 0.58,
    sprintf("R%d × C%d × T%d", row$row_component,
      row$column_component, row$time_component),
    adj = 1, cex = 0.90, font = 2
  )
  text(
    0.98, 0.40,
    sprintf(
      "%s  |  FVE %.4f%%\nreliability %.2f  |  partial R² %.1f%%\nq = %s",
      row$selection_block,
      100 * row$theoretical_variance_fraction,
      row$mean_score_reliability,
      100 * row$effect_r2,
      format_q(row$effect_fdr)
    ),
    adj = 1, cex = 0.66, col = "grey20"
  )
}

draw_spatial_panel <- function(values, limit, title = "") {
  displayed <- values
  displayed[!brain_mask] <- NA_real_
  image(
    seq_len(ncol(displayed)), seq_len(nrow(displayed)),
    t(displayed[nrow(displayed):1L, , drop = FALSE]),
    col = diverging_palette, zlim = c(-limit, limit),
    axes = FALSE, useRaster = TRUE, xlab = "", ylab = "", main = title
  )
  contour(
    seq_len(ncol(brain_mask)), seq_len(nrow(brain_mask)),
    t(brain_mask[nrow(brain_mask):1L, , drop = FALSE] * 1),
    levels = 0.5, add = TRUE, drawlabels = FALSE,
    col = "grey15", lwd = 0.65
  )
  box(col = "grey15", lwd = 0.6)
}

plot_component_grid <- function(rows, output_file, title, direction_label) {
  panel_columns <- length(time_indices) + 1L
  png(
    output_file,
    width = 1960L,
    height = 285L * nrow(rows) + 150L,
    res = 150L,
    bg = "white"
  )
  layout(
    matrix(seq_len(nrow(rows) * panel_columns), nrow = nrow(rows), byrow = TRUE),
    widths = c(1.75, rep(1, length(time_indices)))
  )
  par(
    mar = c(0.6, 0.6, 2.0, 0.6),
    oma = c(2.0, 0.5, 2.8, 0.5),
    cex.main = 0.82,
    font.main = 1
  )

  for (i in seq_len(nrow(rows))) {
    row <- rows[i]
    draw_label_panel(row, direction_label)
    spatial_function <- outer(
      full_eigen$funcs$Phi[, row$row_component],
      full_eigen$funcs$Psi[, row$column_component]
    )
    oriented_sign <- ifelse(row$effect_beta >= 0, 1, -1)
    one_sd_panels <- lapply(time_indices, function(time_index) {
      oriented_sign * sqrt(row$eigenvalue) * spatial_function *
        time_functions[time_index, row$time_component]
    })
    row_limit <- max(abs(unlist(lapply(
      one_sd_panels, function(panel) panel[brain_mask]
    ))))
    if (!is.finite(row_limit) || row_limit <= 0) row_limit <- 1

    for (j in seq_along(one_sd_panels)) {
      panel_title <- if (i == 1L) {
        if (j == 1L) "Baseline" else sprintf("Year %.0f", target_times[[j]])
      } else ""
      draw_spatial_panel(one_sd_panels[[j]], row_limit, panel_title)
    }
  }
  mtext(title, side = 3, outer = TRUE, line = 1.15, cex = 1.15, font = 2)
  mtext(
    sprintf(
      "One-SD product-FPC deviations; red/blue show positive/negative spatial loadings in the %s direction. Each row has its own symmetric color scale.",
      direction_label
    ),
    side = 1, outer = TRUE, line = 0.55, cex = 0.67, col = "grey25"
  )
  dev.off()
}

plot_component_grid(
  age_components,
  file.path(result_dir, "age_associated_fpc_time_panels.png"),
  "Product FPCs most associated with baseline age",
  "older-age"
)
plot_component_grid(
  conversion_components,
  file.path(result_dir, "conversion_associated_fpc_time_panels.png"),
  "Product FPCs most associated with MCI-to-AD conversion",
  "converter"
)

message("Wrote phenotype-oriented FPC panel plots and selection table to ", result_dir)
