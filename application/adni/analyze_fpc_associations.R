#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(rTensor)
})

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]))
} else {
  normalizePath(file.path("application", "adni", "analyze_fpc_associations.R"))
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
setwd(project_root)
source(file.path("R", "matfpca.R"))

result_dir <- file.path("application", "adni", "results")
fit_file <- file.path(result_dir, "cache", "matfpca_fit.rds")
model_file <- file.path(result_dir, "adni_matfpca_model_pcs30.rds")
manifest_file <- file.path(result_dir, "analysis_manifest.csv")

stopifnot(file.exists(fit_file), file.exists(model_file), file.exists(manifest_file))
fit <- readRDS(fit_file)
model <- readRDS(model_file)
manifest <- fread(manifest_file, na.strings = c("", "NA"))

subject_ids <- model$subject_ids
ranks <- fit$ranks
stopifnot(identical(as.integer(ranks), c(34L, 24L, 1L)))

# Recreate the centered, training-visit data used for score estimation. This
# reuses the pooled image cache and the already fitted mean/covariance models;
# no covariance quantity is re-estimated here.
image_cache_files <- Sys.glob(file.path(
  result_dir, "cache",
  sprintf("images_*_pool%02d.rds", model$arguments$pool)
))
image_cache <- NULL
for (candidate_file in image_cache_files) {
  candidate <- readRDS(candidate_file)
  if (identical(as.integer(candidate$image_ids), as.integer(manifest$image_id))) {
    image_cache <- candidate
    break
  }
}
if (is.null(image_cache)) stop("No image cache matches analysis_manifest.csv.")

mean_design <- function(mean_model, times) {
  cbind(
    intercept = 1,
    splines::ns(
      times,
      knots = mean_model$knots,
      intercept = FALSE,
      Boundary.knots = mean_model$Boundary.knots
    )
  )
}

image_matrix <- image_cache$image_matrix
image_matrix[, !as.vector(model$brain_mask)] <- 0
mean_matrix <- mean_design(model$mean_model, manifest$time_year) %*%
  model$mean_model$coefficients
centered_images <- image_matrix - mean_matrix
rm(image_matrix, mean_matrix)

train_rows <- which(manifest$split == "train")
train_by_subject <- split(
  train_rows,
  factor(manifest$PTID[train_rows], levels = subject_ids)
)
X_list <- lapply(train_by_subject, function(rows) {
  subject_array <- array(
    NA_real_,
    dim = c(model$image_dim, length(rows))
  )
  for (j in seq_along(rows)) {
    subject_array[, , j] <- matrix(
      centered_images[rows[[j]], ], nrow = model$image_dim[[1L]]
    )
  }
  rTensor::as.tensor(subject_array)
})
t_idx_list <- lapply(
  train_by_subject,
  function(rows) as.integer(manifest$time_index[rows])
)
mat_data <- create_matrix_funcdata_irreg(
  X_list = X_list,
  t_grid = model$t_grid_years,
  t_idx_list = t_idx_list
)
rm(centered_images, X_list)

# Use the same trapezoidal integral eigensystem as the published temporal
# diagnostic. Unlike the automatic retained fit, this explicitly includes the
# material non-constant second temporal eigenfunction.
trapezoid_weights <- function(grid) {
  c(
    (grid[[2L]] - grid[[1L]]) / 2,
    (grid[3:length(grid)] - grid[1:(length(grid) - 2L)]) / 2,
    (grid[[length(grid)]] - grid[[length(grid) - 1L]]) / 2
  )
}
temporal_weights <- trapezoid_weights(model$t_grid_years)
sqrt_weights <- sqrt(temporal_weights)
weighted_covariance <- sweep(
  sweep(model$covariance$K, 1L, sqrt_weights, "*"),
  2L, sqrt_weights, "*"
)
temporal_eigen <- eigen(
  0.5 * (weighted_covariance + t(weighted_covariance)), symmetric = TRUE
)
temporal_values <- pmax(temporal_eigen$values, 0)
temporal_functions <- sweep(temporal_eigen$vectors, 1L, sqrt_weights, "/")

full_eigen <- eigen.matfcov(model$covariance)
diagnostic_eigen <- create_matrix_eigobj(
  Phi = full_eigen$funcs$Phi[, seq_len(ranks[["row"]]), drop = FALSE],
  Psi = full_eigen$funcs$Psi[, seq_len(ranks[["column"]]), drop = FALSE],
  Theta = temporal_functions[, 1:2, drop = FALSE],
  lam_phi = full_eigen$values$lam_phi[seq_len(ranks[["row"]])],
  lam_psi = full_eigen$values$lam_psi[seq_len(ranks[["column"]])],
  lam_theta = temporal_values[1:2],
  t_grid = model$t_grid_years
)

message("Estimating diagnostic scores at ranks 34 x 24 x 2...")
diagnostic_scores <- get_matfpca_scores(
  mat_data,
  diagnostic_eigen,
  sig2 = model$covariance$sig2,
  K_grid = model$covariance$K
)
score_array <- diagnostic_scores@data
score_matrix <- matrix(score_array, nrow = length(subject_ids))
retained_score_matrix <- matrix(fit$scores@data, nrow = length(subject_ids))
retained_basis_score_correlations <- vapply(
  seq_len(ncol(retained_score_matrix)),
  function(j) cor(score_matrix[, j], retained_score_matrix[, j]),
  numeric(1L)
)

# Average prior-variance reduction from the visit schedule. This is a model-
# based score reliability/estimability measure and does not use age or group.
spatial_values <- diagnostic_eigen$values$lam_phi %o%
  diagnostic_eigen$values$lam_psi
reliability_sum <- array(0, dim = c(dim(spatial_values), 2L))
for (i in seq_along(t_idx_list)) {
  tid <- t_idx_list[[i]]
  Ki <- 0.5 * (
    model$covariance$K[tid, tid, drop = FALSE] +
      t(model$covariance$K[tid, tid, drop = FALSE])
  )
  ee <- eigen(Ki, symmetric = TRUE)
  d <- pmax(ee$values, 0)
  projected_theta <- crossprod(
    ee$vectors,
    diagnostic_eigen$funcs$Theta[tid, 1:2, drop = FALSE]
  )
  for (r in seq_len(nrow(spatial_values))) {
    for (c in seq_len(ncol(spatial_values))) {
      denominator <- d + max(1e-8, model$covariance$sig2 / spatial_values[r, c])
      reliability_sum[r, c, ] <- reliability_sum[r, c, ] +
        colSums(sweep(projected_theta^2, 1L, denominator, "/")) *
          diagnostic_eigen$values$lam_theta[1:2]
    }
  }
}
mean_score_reliability <- reliability_sum / length(t_idx_list)

component_table <- as.data.table(arrayInd(
  seq_len(prod(dim(score_array)[-1L])), dim(score_array)[-1L]
))
setnames(component_table, c("row_component", "column_component", "time_component"))
component_table[, score_column :=
  row_component +
    (column_component - 1L) * ranks[["row"]] +
    (time_component - 1L) * ranks[["row"]] * ranks[["column"]]
]
component_table[, eigenvalue :=
  diagnostic_eigen$values$lam_phi[row_component] *
    diagnostic_eigen$values$lam_psi[column_component] *
    diagnostic_eigen$values$lam_theta[time_component]
]
setorder(component_table, -eigenvalue)
component_table[, product_rank_within_34x24x2 := .I]
component_table[, theoretical_variance_fraction := eigenvalue / sum(eigenvalue)]
component_table[, theoretical_cumulative_variance_fraction :=
  cumsum(theoretical_variance_fraction)
]

subject_table <- unique(manifest[, .(
  PTID, baseline_age, outcome_type, baseline_dx
)])
if (anyDuplicated(subject_table$PTID)) stop("Inconsistent subject metadata.")
subject_table <- subject_table[match(subject_ids, PTID)]
visit_schedule <- manifest[split == "train", .(
  training_visits = .N,
  training_time_span = diff(range(time_year))
), by = PTID]
subject_table <- visit_schedule[subject_table, on = "PTID"]
subject_table[, patient_group := fcase(
  outcome_type == "prevalent_ad_at_baseline", "AD at baseline",
  outcome_type == "incident_ad_conversion", "MCI-to-AD converter",
  outcome_type == "right_censored_non_ad", "MCI non-converter",
  default = NA_character_
)]
stopifnot(!anyNA(subject_table$patient_group))
subject_table[, patient_group := factor(
  patient_group,
  levels = c("MCI non-converter", "MCI-to-AD converter", "AD at baseline")
)]

age <- subject_table$baseline_age
group <- subject_table$patient_group
age_z <- as.numeric(scale(age))
training_visits_z <- as.numeric(scale(subject_table$training_visits))
training_time_span_z <- as.numeric(scale(subject_table$training_time_span))
X_age <- cbind(1, age_z)
X_group <- model.matrix(~ group)
X_full <- model.matrix(~ age_z + group)
X_age_schedule <- model.matrix(
  ~ age_z + training_visits_z + training_time_span_z
)
X_group_schedule <- model.matrix(
  ~ group + training_visits_z + training_time_span_z
)
X_full_schedule <- model.matrix(
  ~ age_z + group + training_visits_z + training_time_span_z
)
n <- length(age)

residual_sse <- function(X, Y) {
  residual <- qr.resid(qr(X), Y)
  colSums(residual^2)
}

Y <- score_matrix[, component_table$score_column, drop = FALSE]
score_sd <- apply(Y, 2L, sd)
score_mean <- colMeans(Y)
Yz <- sweep(sweep(Y, 2L, score_mean, "-"), 2L, score_sd, "/")
sse_null <- colSums(Yz^2)
sse_age <- residual_sse(X_age, Yz)
sse_group <- residual_sse(X_group, Yz)
sse_full <- residual_sse(X_full, Yz)
coef_full <- qr.coef(qr(X_full), Yz)
sse_age_schedule <- residual_sse(X_age_schedule, Yz)
sse_group_schedule <- residual_sse(X_group_schedule, Yz)
sse_full_schedule <- residual_sse(X_full_schedule, Yz)
coef_full_schedule <- qr.coef(qr(X_full_schedule), Yz)

component_table[, `:=`(
  empirical_score_sd = score_sd,
  empirical_score_variance_to_eigenvalue = score_sd^2 / eigenvalue,
  mean_score_reliability = as.vector(mean_score_reliability)[score_column],
  age_partial_r2_given_group = pmax((sse_group - sse_full) / sse_group, 0),
  age_adjusted_standardized_beta = coef_full["age_z", ],
  group_eta2_unadjusted = pmax((sse_null - sse_group) / sse_null, 0),
  group_partial_r2_given_age = pmax((sse_age - sse_full) / sse_age, 0)
)]

df_full <- n - ncol(X_full)
age_f <- ((sse_group - sse_full) / 1) / (sse_full / df_full)
group_f <- ((sse_age - sse_full) / 2) / (sse_full / df_full)
component_table[, `:=`(
  age_adjusted_p = pf(age_f, 1, df_full, lower.tail = FALSE),
  group_age_adjusted_p = pf(group_f, 2, df_full, lower.tail = FALSE)
)]

df_full_schedule <- n - ncol(X_full_schedule)
age_schedule_f <-
  ((sse_group_schedule - sse_full_schedule) / 1) /
    (sse_full_schedule / df_full_schedule)
group_schedule_f <-
  ((sse_age_schedule - sse_full_schedule) / 2) /
    (sse_full_schedule / df_full_schedule)
component_table[, `:=`(
  age_partial_r2_given_group_and_schedule = pmax(
    (sse_group_schedule - sse_full_schedule) / sse_group_schedule, 0
  ),
  age_schedule_adjusted_standardized_beta =
    coef_full_schedule["age_z", ],
  age_schedule_adjusted_p =
    pf(age_schedule_f, 1, df_full_schedule, lower.tail = FALSE),
  group_partial_r2_given_age_and_schedule = pmax(
    (sse_age_schedule - sse_full_schedule) / sse_age_schedule, 0
  ),
  group_age_schedule_adjusted_p =
    pf(group_schedule_f, 2, df_full_schedule, lower.tail = FALSE)
)]

spearman <- lapply(seq_len(ncol(Y)), function(j) {
  test <- suppressWarnings(cor.test(age, Y[, j], method = "spearman", exact = FALSE))
  c(rho = unname(test$estimate), p = test$p.value)
})
spearman <- do.call(rbind, spearman)
component_table[, `:=`(
  age_spearman_rho = spearman[, "rho"],
  age_spearman_p = spearman[, "p"]
)]

# A clinically focused binary contrast among subjects who entered with MCI.
mci <- group != "AD at baseline"
mci_converter <- group[mci] == "MCI-to-AD converter"
mci_effect <- lapply(seq_len(ncol(Y)), function(j) {
  values <- Y[mci, j]
  x1 <- values[mci_converter]
  x0 <- values[!mci_converter]
  pooled_sd <- sqrt(
    ((length(x1) - 1) * var(x1) + (length(x0) - 1) * var(x0)) /
      (length(x1) + length(x0) - 2)
  )
  test <- suppressWarnings(wilcox.test(x1, x0, exact = FALSE))
  c(
    standardized_mean_difference = (mean(x1) - mean(x0)) / pooled_sd,
    p = test$p.value
  )
})
mci_effect <- do.call(rbind, mci_effect)
component_table[, `:=`(
  mci_converter_standardized_mean_difference =
    mci_effect[, "standardized_mean_difference"],
  mci_converter_wilcoxon_p = mci_effect[, "p"]
)]

X_mci_base <- model.matrix(
  ~ age_z[mci] + training_visits_z[mci] + training_time_span_z[mci]
)
X_mci_full <- model.matrix(
  ~ age_z[mci] + mci_converter + training_visits_z[mci] +
    training_time_span_z[mci]
)
sse_mci_base <- residual_sse(X_mci_base, Yz[mci, , drop = FALSE])
sse_mci_full <- residual_sse(X_mci_full, Yz[mci, , drop = FALSE])
coef_mci_full <- qr.coef(qr(X_mci_full), Yz[mci, , drop = FALSE])
df_mci_full <- sum(mci) - ncol(X_mci_full)
mci_converter_f <- ((sse_mci_base - sse_mci_full) / 1) /
  (sse_mci_full / df_mci_full)
component_table[, `:=`(
  mci_converter_partial_r2_given_age_and_schedule = pmax(
    (sse_mci_base - sse_mci_full) / sse_mci_base, 0
  ),
  mci_converter_age_schedule_adjusted_standardized_beta =
    coef_mci_full["mci_converterTRUE", ],
  mci_converter_age_schedule_adjusted_p =
    pf(mci_converter_f, 1, df_mci_full, lower.tail = FALSE)
)]

group_means <- rowsum(Yz, group, reorder = FALSE) /
  as.numeric(table(group))
component_table[, `:=`(
  mean_z_mci_nonconverter = group_means["MCI non-converter", ],
  mean_z_mci_converter = group_means["MCI-to-AD converter", ],
  mean_z_ad_baseline = group_means["AD at baseline", ]
)]

component_table[, `:=`(
  age_spearman_fdr_all = p.adjust(age_spearman_p, "BH"),
  age_adjusted_fdr_all = p.adjust(age_adjusted_p, "BH"),
  age_schedule_adjusted_fdr_all = p.adjust(age_schedule_adjusted_p, "BH"),
  group_age_adjusted_fdr_all = p.adjust(group_age_adjusted_p, "BH"),
  group_age_schedule_adjusted_fdr_all =
    p.adjust(group_age_schedule_adjusted_p, "BH"),
  mci_converter_fdr_all = p.adjust(mci_converter_wilcoxon_p, "BH"),
  mci_converter_age_schedule_adjusted_fdr_all =
    p.adjust(mci_converter_age_schedule_adjusted_p, "BH")
)]
component_table[, `:=`(
  age_spearman_fdr_within_time = p.adjust(age_spearman_p, "BH"),
  age_schedule_adjusted_fdr_within_time =
    p.adjust(age_schedule_adjusted_p, "BH"),
  group_age_adjusted_fdr_within_time = p.adjust(group_age_adjusted_p, "BH"),
  group_age_schedule_adjusted_fdr_within_time =
    p.adjust(group_age_schedule_adjusted_p, "BH"),
  mci_converter_fdr_within_time = p.adjust(mci_converter_wilcoxon_p, "BH"),
  mci_converter_age_schedule_adjusted_fdr_within_time =
    p.adjust(mci_converter_age_schedule_adjusted_p, "BH")
), by = time_component]

component_table[, component := sprintf(
  "R%d_C%d_T%d", row_component, column_component, time_component
)]

# Quantify whether each separable spatial loading is diffuse or concentrated.
# These are coordinate summaries, not anatomical labels.
mask <- model$brain_mask
mask_pixels <- sum(mask)
spatial_localization <- unique(component_table[, .(
  row_component, column_component
)])
localization_rows <- lapply(seq_len(nrow(spatial_localization)), function(i) {
  r <- spatial_localization$row_component[[i]]
  c <- spatial_localization$column_component[[i]]
  loading <- diagnostic_eigen$funcs$Phi[, r] %o%
    diagnostic_eigen$funcs$Psi[, c]
  energy <- loading^2
  energy[!mask] <- 0
  energy <- energy / sum(energy)
  ordered_energy <- sort(energy[mask], decreasing = TRUE)
  coordinates <- which(mask, arr.ind = TRUE)
  peak <- which.max(energy)
  peak_coordinate <- arrayInd(peak, dim(energy))
  data.table(
    row_component = r,
    column_component = c,
    spatial_energy50_pixel_fraction =
      which(cumsum(ordered_energy) >= 0.5)[[1L]] / mask_pixels,
    spatial_effective_pixel_fraction =
      (1 / sum(ordered_energy^2)) / mask_pixels,
    energy_centroid_row = sum(coordinates[, 1L] * energy[mask]),
    energy_centroid_column = sum(coordinates[, 2L] * energy[mask]),
    peak_loading_row = peak_coordinate[[1L]],
    peak_loading_column = peak_coordinate[[2L]]
  )
})
spatial_localization <- rbindlist(localization_rows)
component_table <- spatial_localization[component_table,
  on = .(row_component, column_component)]
setcolorder(component_table, c(
  "component", "row_component", "column_component", "time_component",
  "product_rank_within_34x24x2", "eigenvalue",
  "theoretical_variance_fraction", "theoretical_cumulative_variance_fraction",
  "empirical_score_sd", "empirical_score_variance_to_eigenvalue",
  "mean_score_reliability",
  setdiff(names(component_table), c(
    "component", "row_component", "column_component", "time_component",
    "product_rank_within_34x24x2", "eigenvalue",
    "theoretical_variance_fraction", "theoretical_cumulative_variance_fraction",
    "empirical_score_sd", "empirical_score_variance_to_eigenvalue",
    "mean_score_reliability"
  ))
))

fwrite(component_table, file.path(result_dir, "fpc_association_all_components.csv"))
fwrite(
  component_table[time_component == 2L][order(-eigenvalue)],
  file.path(result_dir, "fpc_association_temporal2.csv")
)

top_line <- function(table, columns, count = 5L) {
  selected <- table[seq_len(min(count, nrow(table))), ..columns]
  paste(capture.output(print(selected)), collapse = "\n")
}

summary_lines <- c(
  sprintf("Retained marginal ranks: row=%d, column=%d, time=%d", ranks[[1L]], ranks[[2L]], ranks[[3L]]),
  sprintf("Diagnostic temporal-2 marginal variance fraction: %.6f", temporal_values[[2L]] / sum(temporal_values)),
  sprintf(
    "Absolute score correlation between retained T1 and diagnostic-basis T1: median=%.6f; minimum=%.6f",
    median(abs(retained_basis_score_correlations)),
    min(abs(retained_basis_score_correlations))
  ),
  sprintf("Subjects: %d", n),
  paste("Patient groups:", paste(names(table(group)), as.integer(table(group)), sep = "=", collapse = "; ")),
  "",
  "Top age associations by absolute Spearman rho:",
  top_line(
    component_table[order(-abs(age_spearman_rho))],
    c("component", "age_spearman_rho", "age_spearman_fdr_all", "age_partial_r2_given_group_and_schedule", "age_schedule_adjusted_fdr_all")
  ),
  "",
  "Top age-adjusted three-group associations by partial R2:",
  top_line(
    component_table[order(-group_partial_r2_given_age_and_schedule)],
    c("component", "group_partial_r2_given_age_and_schedule", "group_age_schedule_adjusted_fdr_all", "mean_z_mci_nonconverter", "mean_z_mci_converter", "mean_z_ad_baseline")
  ),
  "",
  "Top temporal-2 age associations by absolute Spearman rho:",
  top_line(
    component_table[time_component == 2L][order(-abs(age_spearman_rho))],
    c("component", "product_rank_within_34x24x2", "theoretical_variance_fraction", "mean_score_reliability", "age_spearman_rho", "age_partial_r2_given_group_and_schedule", "age_schedule_adjusted_fdr_within_time")
  ),
  "",
  "Top temporal-2 age-adjusted group associations by partial R2:",
  top_line(
    component_table[time_component == 2L][order(-group_partial_r2_given_age_and_schedule)],
    c("component", "product_rank_within_34x24x2", "theoretical_variance_fraction", "mean_score_reliability", "group_partial_r2_given_age_and_schedule", "group_age_schedule_adjusted_fdr_within_time", "mean_z_mci_nonconverter", "mean_z_mci_converter", "mean_z_ad_baseline")
  ),
  "",
  "Top temporal-2 converter/non-converter associations after adjusting for age and visit schedule:",
  top_line(
    component_table[time_component == 2L][order(-mci_converter_partial_r2_given_age_and_schedule)],
    c("component", "product_rank_within_34x24x2", "theoretical_variance_fraction", "mean_score_reliability", "mci_converter_standardized_mean_difference", "mci_converter_partial_r2_given_age_and_schedule", "mci_converter_age_schedule_adjusted_fdr_within_time")
  )
)
writeLines(summary_lines, file.path(result_dir, "fpc_association_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
