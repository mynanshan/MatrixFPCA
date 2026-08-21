#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(rTensor)
})
if (!nzchar(Sys.getenv("TZ"))) Sys.setenv(TZ = "America/Vancouver")

required_packages <- c("mgcv", "pracma", "foreach", "splines", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]))
} else {
  normalizePath(file.path("application", "adni-latent", "adni_latent_matfpca.R"))
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
setwd(project_root)
source(file.path("R", "matfpca.R"))
source(file.path("application", "adni-latent", "adni_latent_utils.R"))

defaults <- list(
  metadata = file.path(
    "application", "adni-latent", "latents", "brainiac_patch_full", "metadata.csv"
  ),
  output_dir = file.path(
    "application", "adni-latent", "results", "brainiac_patch_full"
  ),
  seed = 20260819L,
  holdout_fraction = 0.20,
  max_month = 48,
  prediction_step_months = 6,
  feature_rank = 64L,
  mean_df = 5L,
  weight_scheme = "identity",
  weight_folds = 3L,
  opw_iter = 1L,
  nbasis = "theory",
  mode_pve = 0.90,
  row_rank_cap = 24L,
  column_rank_cap = 24L,
  time_rank_cap = 4L,
  top_scores = 20L,
  diagnostic_coordinates = 192L,
  reconstructed_subjects = 3L,
  retrieval_candidates = 250L,
  reconstruction_analysis = FALSE,
  cache = TRUE
)

usage <- function() {
  cat(
    "Usage: Rscript application/adni-latent/adni_latent_matfpca.R [--option=value ...]\n\n",
    "Important options:\n",
    "  --metadata=PATH                 BrainIAC extraction metadata.csv\n",
    "  --output-dir=PATH               Analysis output directory\n",
    "  --feature-rank=64               Training-only feature PCA dimension\n",
    "  --holdout-fraction=0.20         Subject-level validation fraction; one interior target each\n",
    "  --weight-scheme=identity        identity or gaussian_refined\n",
    "  --mode-pve=0.90                 Marginal rank target before rank caps\n",
    "  --row-rank-cap=24               Maximum retained token-mode rank\n",
    "  --column-rank-cap=24            Maximum retained feature-mode rank\n",
    "  --time-rank-cap=4               Maximum retained temporal rank\n",
    "  --prediction-step-months=6      Complete-grid reconstruction spacing\n",
    "  --reconstruction-analysis=false Skip latent/image reconstruction outputs\n",
    "  --cache=true                    Reuse centered tensors and fitted model\n",
    sep = ""
  )
}

coerce_option <- function(value, template) {
  if (is.logical(template)) {
    lowered <- tolower(value)
    if (!lowered %in% c("true", "false", "1", "0", "yes", "no")) {
      stop("Logical option must be true or false, not: ", value)
    }
    return(lowered %in% c("true", "1", "yes"))
  }
  if (is.integer(template)) return(as.integer(value))
  if (is.numeric(template)) return(as.numeric(value))
  value
}

parse_options <- function(default_values) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (any(arguments %in% c("--help", "-h"))) {
    usage()
    quit(save = "no", status = 0)
  }
  result <- default_values
  for (argument in arguments) {
    if (!grepl("^--[^=]+=", argument)) {
      stop("Options must use --name=value syntax: ", argument)
    }
    pair <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", pair[[1L]], fixed = TRUE)
    if (!key %in% names(result)) stop("Unknown option: --", pair[[1L]])
    result[[key]] <- coerce_option(
      paste(pair[-1L], collapse = "="), default_values[[key]]
    )
  }
  result
}

args <- parse_options(defaults)
if (!args$weight_scheme %in% c("identity", "gaussian_refined")) {
  stop("--weight-scheme must be identity or gaussian_refined.")
}
if (!is.finite(args$holdout_fraction) || args$holdout_fraction < 0 ||
    args$holdout_fraction >= 1) {
  stop("--holdout-fraction must lie in [0, 1).")
}
positive_integer_options <- c(
  "feature_rank", "mean_df", "weight_folds", "opw_iter", "row_rank_cap",
  "column_rank_cap", "time_rank_cap", "top_scores", "diagnostic_coordinates",
  "reconstructed_subjects", "retrieval_candidates"
)
if (any(vapply(args[positive_integer_options], function(value) {
  !is.finite(value) || value < 1L
}, logical(1)))) stop("Rank, count, and iteration options must be positive.")
if (args$mean_df < 3L || args$mode_pve <= 0 || args$mode_pve > 1 ||
    args$max_month <= 0 || args$prediction_step_months <= 0) {
  stop("Invalid mean, PVE, or temporal-grid setting.")
}

log_message <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
  flush.console()
}

args$metadata <- normalizePath(args$metadata)
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
args$output_dir <- normalizePath(args$output_dir)
run_summary_path <- file.path(args$output_dir, "run_summary.json")
jsonlite::write_json(
  list(
    status = "in_progress",
    started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    metadata = args$metadata,
    configuration = args
  ),
  run_summary_path,
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null"
)
cache_dir <- file.path(args$output_dir, "cache")
tensor_dir <- file.path(cache_dir, "centered_subjects")
dir.create(tensor_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- fread(args$metadata, na.strings = c("", "NA"))
# Accept the explicit names emitted by extract_brainiac_latents.py while also
# keeping compatibility with early pilot manifests used during development.
column_aliases <- c(
  time_month = "model_time_month",
  viscode = "mri_VISCODE2",
  latent_cols = "latent_columns"
)
for (canonical_name in names(column_aliases)) {
  source_name <- column_aliases[[canonical_name]]
  if (!canonical_name %in% names(metadata) && source_name %in% names(metadata)) {
    metadata[, (canonical_name) := get(source_name)]
  }
}
required_columns <- c(
  "subject_id", "visit_id", "image_id", "time_month", "viscode",
  "latent_path", "nifti_path", "latent_rows", "latent_cols"
)
missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns)) {
  stop("Extraction metadata is missing: ", paste(missing_columns, collapse = ", "))
}
metadata[, subject_id := as.character(subject_id)]
metadata[, latent_path := vapply(
  latent_path, resolve_latent_path, character(1L), metadata_path = args$metadata
)]
metadata[, nifti_path := vapply(
  nifti_path,
  function(path) {
    if (grepl("^/", path)) path else
      normalizePath(file.path(dirname(args$metadata), path), mustWork = FALSE)
  },
  character(1L)
)]
if (any(!file.exists(metadata$latent_path))) stop("One or more latent files are missing.")
if (any(!file.exists(metadata$nifti_path))) stop("One or more source NIfTI files are missing.")
if (uniqueN(metadata[, .(latent_rows, latent_cols)]) != 1L) {
  stop("All latent matrices must have the same shape.")
}
if (any(!is.finite(metadata$time_month)) || any(metadata$time_month > args$max_month)) {
  stop("Metadata contains invalid time or a visit after --max-month.")
}

# Use the same day-resolution time convention as the image-space application.
metadata[, time_tick := as.integer(round(time_month * 365.25 / 12))]
metadata[tolower(viscode) == "bl", time_tick := 0L]
metadata[, time_year := time_tick / 365.25]
if (metadata[, anyDuplicated(paste(subject_id, time_tick, sep = "|"))]) {
  stop("A subject has multiple visits on the same snapped analysis day.")
}
setorder(metadata, subject_id, time_tick, image_id)
metadata[, analysis_row := .I]
visit_counts <- metadata[, .N, by = subject_id]
if (any(visit_counts$N < 3L)) {
  stop("Interior-visit interpolation requires at least three visits per subject.")
}

# Split at subject level before estimating any population quantity.  For each
# validation subject, only the non-target visits are later used to estimate
# that subject's scores; none of its visits enter feature PCA, the mean, or the
# population covariance fit.
subject_ids <- sort(unique(metadata$subject_id))
holdout_count <- if (args$holdout_fraction == 0) 0L else
  max(1L, as.integer(round(args$holdout_fraction * length(subject_ids))))
set.seed(args$seed)
validation_subjects <- if (holdout_count) {
  sample(subject_ids, min(holdout_count, length(subject_ids)))
} else character()
population_subject_ids <- setdiff(subject_ids, validation_subjects)
if (length(population_subject_ids) < 2L) {
  stop("The subject-level split leaves fewer than two population-training subjects.")
}
metadata[, cohort_split := ifelse(
  subject_id %in% validation_subjects, "validation", "population_train"
)]
metadata[, split := ifelse(
  cohort_split == "validation", "validation_conditioning", "population_train"
)]
set.seed(args$seed + 1L)
for (subject in validation_subjects) {
  rows <- which(metadata$subject_id == subject)
  interior <- rows[seq.int(2L, length(rows) - 1L)]
  metadata[sample(interior, 1L), split := "validation_test"]
}
if (length(validation_subjects) && any(
  metadata[
    cohort_split == "validation" & split == "validation_conditioning",
    .N, by = subject_id
  ]$N < 2L
)) {
  stop("A validation subject has fewer than two conditioning visits after holdout.")
}

prediction_ticks <- as.integer(round(
  seq(0, args$max_month, by = args$prediction_step_months) * 365.25 / 12
))
prediction_ticks <- sort(unique(c(
  prediction_ticks, 0L,
  round(args$max_month * 365.25 / 12)
)))
grid_ticks <- sort(unique(c(metadata$time_tick, prediction_ticks)))
t_grid <- grid_ticks / 365.25
metadata[, time_index := match(time_tick, grid_ticks)]
prediction_years <- prediction_ticks / 365.25

fwrite(metadata, file.path(args$output_dir, "analysis_manifest.csv"))
metadata_hash <- unname(tools::md5sum(args$metadata))
latent_file_information <- file.info(metadata$latent_path)
if (anyNA(latent_file_information[, c("size", "mtime")])) {
  stop("Could not fingerprint one or more latent files.")
}
latent_file_hashes <- unname(tools::md5sum(metadata$latent_path))
if (anyNA(latent_file_hashes)) {
  stop("Could not hash one or more latent files.")
}
analysis_code_files <- c(
  file.path("application", "adni-latent", "adni_latent_matfpca.R"),
  file.path("application", "adni-latent", "adni_latent_utils.R"),
  file.path("application", "adni-latent", "latent_io.R"),
  file.path("R", "matfpca.R"),
  file.path("R", "matfdata.R")
)
analysis_code_hashes <- unname(tools::md5sum(analysis_code_files))
if (anyNA(analysis_code_hashes)) {
  stop("Could not fingerprint one or more analysis source files.")
}
configuration_signature <- paste(
  "adni_brainiac_patch_matfpca_v4_subject_weighted_feature_mean",
  metadata_hash,
  paste(analysis_code_hashes, collapse = ","),
  paste(latent_file_information$size, latent_file_hashes, collapse = ","),
  paste(metadata$visit_id, metadata$split, collapse = ","),
  args$seed, args$holdout_fraction, args$max_month,
  args$feature_rank, args$mean_df, args$weight_scheme,
  args$weight_folds, args$opw_iter, args$nbasis, args$mode_pve,
  args$row_rank_cap, args$column_rank_cap, args$time_rank_cap,
  args$prediction_step_months,
  sep = "|"
)

preprocessing_cache_path <- file.path(cache_dir, "latent_preprocessing.rds")
preprocessing <- NULL
if (args$cache && file.exists(preprocessing_cache_path)) {
  candidate <- readRDS(preprocessing_cache_path)
  if (identical(candidate$signature, configuration_signature) &&
      all(file.exists(candidate$subject_tensor_paths))) {
    preprocessing <- candidate
    log_message("Reusing latent preprocessing cache.")
  }
}

if (is.null(preprocessing)) {
  training_metadata <- metadata[split == "population_train"]
  log_message(
    "Fitting population-training-only feature PCA (", nrow(training_metadata),
    " visits / ", length(population_subject_ids), " subjects, target rank ",
    args$feature_rank, ")..."
  )
  feature_projection <- fit_feature_projection(
    training_metadata, args$feature_rank
  )
  log_message(
    "Feature PCA retained ", feature_projection$retained_rank,
    " / ", feature_projection$original_feature_count, " coordinates (",
    sprintf("%.2f%% PVE", 100 * feature_projection$retained_pve), ")."
  )
  mean_boundary <- range(t_grid)
  log_message("Fitting the population-training-only subject-weighted temporal mean...")
  mean_model <- fit_latent_mean_streaming(
    training_metadata, feature_projection, args$mean_df, mean_boundary
  )

  conditioning_by_subject <- split(
    which(metadata$split != "validation_test"),
    factor(
      metadata$subject_id[metadata$split != "validation_test"],
      levels = subject_ids
    )
  )
  subject_tensor_paths <- file.path(
    tensor_dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", subject_ids), ".rds")
  )
  log_message("Centering and caching one tensor per subject...")
  for (subject_position in seq_along(subject_ids)) {
    rows <- conditioning_by_subject[[subject_position]]
    subject_array <- array(
      NA_real_,
      dim = c(
        feature_projection$token_count,
        feature_projection$retained_rank,
        length(rows)
      )
    )
    for (visit_position in seq_along(rows)) {
      row <- rows[[visit_position]]
      projected <- project_latent_matrix(
        read_npy_matrix(metadata$latent_path[[row]]), feature_projection
      )
      subject_array[, , visit_position] <- projected -
        predict_latent_mean(mean_model, metadata$time_year[[row]])
    }
    saveRDS(rTensor::as.tensor(subject_array), subject_tensor_paths[[subject_position]],
            compress = FALSE)
  }
  t_idx_list <- lapply(
    conditioning_by_subject,
    function(rows) as.integer(metadata$time_index[rows])
  )
  preprocessing <- list(
    signature = configuration_signature,
    feature_projection = feature_projection,
    mean_model = mean_model,
    subject_tensor_paths = subject_tensor_paths,
    t_idx_list = t_idx_list
  )
  if (args$cache) saveRDS(preprocessing, preprocessing_cache_path, compress = FALSE)
}

feature_projection <- preprocessing$feature_projection
mean_model <- preprocessing$mean_model
subject_tensor_paths <- preprocessing$subject_tensor_paths
t_idx_list <- preprocessing$t_idx_list
mat_data <- create_matrix_funcdata_irreg(
  X_list = subject_tensor_paths,
  t_grid = t_grid,
  t_idx_list = t_idx_list
)
population_positions <- match(population_subject_ids, subject_ids)
population_mat_data <- create_matrix_funcdata_irreg(
  X_list = subject_tensor_paths[population_positions],
  t_grid = t_grid,
  t_idx_list = t_idx_list[population_positions]
)

fit_cache_path <- file.path(cache_dir, "matfpca_fit.rds")
analysis_fit <- NULL
if (args$cache && file.exists(fit_cache_path)) {
  candidate <- readRDS(fit_cache_path)
  if (identical(candidate$signature, configuration_signature)) {
    analysis_fit <- candidate
    log_message("Reusing MatFPCA fit cache.")
  }
}

if (is.null(analysis_fit)) {
  fit_nbasis <- if (grepl("^[0-9]+$", args$nbasis)) as.integer(args$nbasis) else args$nbasis
  fold_id <- NULL
  if (args$weight_scheme == "gaussian_refined") {
    args$weight_folds <- min(args$weight_folds, length(population_subject_ids))
    set.seed(args$seed + 2L)
    fold_id <- sample(rep(
      seq_len(args$weight_folds), length.out = length(population_subject_ids)
    ))
  }
  log_message(
    "Fitting MatFPCA (", args$weight_scheme, ", representation ",
    paste(population_mat_data$vdim, collapse = " x "), ", grid ",
    length(t_grid), ", population subjects ", length(population_subject_ids), ")..."
  )
  timing <- system.time({
    covariance_fit <- fit_matfcov(
      population_mat_data,
      wScheme = args$weight_scheme,
      noise = TRUE,
      nbasis = fit_nbasis,
      opwIter = args$opw_iter,
      temporal_method = "two_step",
      nfolds = args$weight_folds,
      fold_id = fold_id,
      crossfit_diagonal = "refined"
    )
  })
  eig_full <- eigen_matfcov_quadrature(covariance_fit$C)
  ranks <- c(
    row = rank_for_pve(eig_full$values$lam_phi, args$mode_pve, args$row_rank_cap),
    column = rank_for_pve(
      eig_full$values$lam_psi, args$mode_pve, args$column_rank_cap
    ),
    time = rank_for_pve(
      eig_full$values$lam_theta, args$mode_pve, args$time_rank_cap
    )
  )
  eig <- trunc.matfeigen(
    eig_full, I = ranks[["row"]], J = ranks[["column"]],
    K = ranks[["time"]], verbose = FALSE
  )
  log_message(
    "Estimating subject scores from population or validation-conditioning visits at ranks ",
    paste(ranks, collapse = " x "), "..."
  )
  scores <- get_matfpca_scores(
    mat_data, eig,
    sig2 = covariance_fit$C$sig2,
    K_grid = covariance_fit$C$K
  )
  analysis_fit <- list(
    signature = configuration_signature,
    covariance = covariance_fit$C,
    eig_full = eig_full,
    eig = eig,
    ranks = ranks,
    scores = scores,
    elapsed_seconds = unname(timing[["elapsed"]])
  )
  if (args$cache) saveRDS(analysis_fit, fit_cache_path, compress = FALSE)
}

covariance <- analysis_fit$covariance
eig_full <- analysis_fit$eig_full
eig <- analysis_fit$eig
ranks <- analysis_fit$ranks
scores <- analysis_fit$scores

model <- list(
  model_type = "ADNI frozen BrainIAC patch-token MatFPCA",
  subject_ids = subject_ids,
  population_subject_ids = population_subject_ids,
  validation_subject_ids = validation_subjects,
  t_grid = t_grid,
  prediction_years = prediction_years,
  covariance = covariance,
  eig = eig,
  ranks = ranks,
  scores = scores@data,
  mean_model = mean_model,
  feature_projection = feature_projection,
  representation_shape = c(
    feature_projection$token_count, feature_projection$retained_rank
  ),
  original_representation_shape = c(
    feature_projection$token_count, feature_projection$original_feature_count
  ),
  token_grid = c(6L, 6L, 6L),
  metadata_path = args$metadata,
  configuration = args,
  configuration_signature = configuration_signature
)

# -------------------------------------------------------------------------
# Representation diagnostics before covariance interpretation

load_projected_visit <- function(row) {
  project_latent_matrix(
    read_npy_matrix(metadata$latent_path[[row]]), feature_projection
  )
}

normalized_distance <- function(left, right) {
  sqrt(mean((left - right)^2))
}

distance_rows <- list()
distance_position <- 0L
for (subject in subject_ids) {
  rows <- which(metadata$subject_id == subject)
  visits <- lapply(rows, load_projected_visit)
  if (length(rows) < 2L) next
  pairs <- utils::combn(seq_along(rows), 2L)
  for (pair_position in seq_len(ncol(pairs))) {
    left_position <- pairs[1L, pair_position]
    right_position <- pairs[2L, pair_position]
    left <- visits[[left_position]]
    right <- visits[[right_position]]
    lag <- abs(
      metadata$time_year[rows[[left_position]]] -
        metadata$time_year[rows[[right_position]]]
    )
    representations <- list(
      patch_matrix = c(left),
      mean_pooled = colMeans(left),
      upstream_first_patch = left[1L, ]
    )
    right_representations <- list(
      patch_matrix = c(right),
      mean_pooled = colMeans(right),
      upstream_first_patch = right[1L, ]
    )
    for (representation in names(representations)) {
      distance_position <- distance_position + 1L
      distance_rows[[distance_position]] <- data.table(
        pair_type = "within_subject",
        representation = representation,
        subject_left = subject,
        subject_right = subject,
        lag_years = lag,
        distance = normalized_distance(
          representations[[representation]],
          right_representations[[representation]]
        )
      )
    }
  }
}

within_pair_count <- sum(vapply(
  split(metadata$analysis_row, metadata$subject_id),
  function(rows) choose(length(rows), 2L),
  numeric(1L)
))
set.seed(args$seed + 3L)
between_pairs <- matrix(integer(), ncol = 2L)
attempts <- 0L
while (nrow(between_pairs) < within_pair_count && attempts < 20L * within_pair_count) {
  candidate <- sample(seq_len(nrow(metadata)), 2L)
  attempts <- attempts + 1L
  if (metadata$subject_id[candidate[[1L]]] == metadata$subject_id[candidate[[2L]]]) next
  between_pairs <- rbind(between_pairs, candidate)
}
if (nrow(between_pairs)) {
  between_pairs <- unique(between_pairs)
  between_pairs <- between_pairs[seq_len(min(nrow(between_pairs), within_pair_count)), , drop = FALSE]
  between_cache <- new.env(parent = emptyenv())
  get_between_visit <- function(row) {
    key <- as.character(row)
    if (!exists(key, envir = between_cache, inherits = FALSE)) {
      assign(key, load_projected_visit(row), envir = between_cache)
    }
    get(key, envir = between_cache, inherits = FALSE)
  }
  for (pair_position in seq_len(nrow(between_pairs))) {
    left_row <- between_pairs[pair_position, 1L]
    right_row <- between_pairs[pair_position, 2L]
    left <- get_between_visit(left_row)
    right <- get_between_visit(right_row)
    representations <- list(
      patch_matrix = c(left),
      mean_pooled = colMeans(left),
      upstream_first_patch = left[1L, ]
    )
    right_representations <- list(
      patch_matrix = c(right),
      mean_pooled = colMeans(right),
      upstream_first_patch = right[1L, ]
    )
    for (representation in names(representations)) {
      distance_position <- distance_position + 1L
      distance_rows[[distance_position]] <- data.table(
        pair_type = "between_subject",
        representation = representation,
        subject_left = metadata$subject_id[[left_row]],
        subject_right = metadata$subject_id[[right_row]],
        lag_years = NA_real_,
        distance = normalized_distance(
          representations[[representation]],
          right_representations[[representation]]
        )
      )
    }
  }
}
distance_table <- rbindlist(distance_rows)
fwrite(distance_table, file.path(args$output_dir, "latent_distance_pairs.csv"))
distance_summary <- distance_table[, .(
  pairs = .N,
  median_distance = median(distance),
  mean_distance = mean(distance),
  lag_distance_spearman = if (pair_type[[1L]] == "within_subject") {
    suppressWarnings(cor(lag_years, distance, method = "spearman"))
  } else NA_real_
), by = .(pair_type, representation)]
fwrite(distance_summary, file.path(args$output_dir, "latent_distance_summary.csv"))

png(file.path(args$output_dir, "latent_distance_diagnostics.png"),
    width = 1500, height = 650, res = 150)
par(mfrow = c(1L, 2L), mar = c(4.2, 4.4, 2.4, 1.2))
matrix_within <- distance_table[
  pair_type == "within_subject" & representation == "patch_matrix"
]
plot(
  matrix_within$lag_years, matrix_within$distance,
  pch = 16, col = grDevices::adjustcolor("#2166AC", 0.45),
  xlab = "Within-subject visit lag (years)",
  ylab = "Normalized patch-matrix distance",
  main = "Longitudinal signal"
)
if (nrow(matrix_within) >= 4L) {
  lines(stats::lowess(matrix_within$lag_years, matrix_within$distance),
        col = "#B2182B", lwd = 2.5)
}
boxplot(
  distance ~ interaction(pair_type, representation, drop = TRUE),
  data = distance_table,
  las = 2, outline = FALSE, col = c("#92C5DE", "#F4A582"),
  ylab = "Normalized latent distance", xlab = "", main = "Within vs between subjects"
)
dev.off()

# Sampled same-visit covariance check. This avoids ever materializing the
# (p*q) x (p*q) empirical covariance while still testing the fitted Kronecker
# form. The identity MatFPCA estimator pools visits/pairs, so this diagnostic
# deliberately uses the same visit-pooled estimand (feature PCA and the mean
# remain subject-weighted).
coordinate_count <- min(args$diagnostic_coordinates, prod(model$representation_shape))
set.seed(args$seed + 4L)
sampled_coordinates <- sort(sample(seq_len(prod(model$representation_shape)), coordinate_count))
sampled_rows <- ((sampled_coordinates - 1L) %% model$representation_shape[[1L]]) + 1L
sampled_columns <- ((sampled_coordinates - 1L) %/% model$representation_shape[[1L]]) + 1L
training_rows <- which(metadata$split == "population_train")
sampled_values <- matrix(NA_real_, nrow = length(training_rows), ncol = coordinate_count)
for (position in seq_along(training_rows)) {
  row <- training_rows[[position]]
  centered <- load_projected_visit(row) -
    predict_latent_mean(mean_model, metadata$time_year[[row]])
  sampled_values[position, ] <- centered[sampled_coordinates]
}
training_subject_counts <- table(metadata$subject_id[training_rows])
sampled_center <- colMeans(sampled_values)
sampled_values <- sweep(sampled_values, 2L, sampled_center, "-")
empirical_sampled_covariance <- crossprod(sampled_values) /
  max(1L, nrow(sampled_values) - 1L)
mean_temporal_variance <- mean(diag(covariance$K)[metadata$time_index[training_rows]])
model_sampled_covariance <- mean_temporal_variance *
  covariance$U[sampled_rows, sampled_rows, drop = FALSE] *
  covariance$V[sampled_columns, sampled_columns, drop = FALSE]
diag(model_sampled_covariance) <- diag(model_sampled_covariance) + covariance$sig2
off_diagonal <- upper.tri(empirical_sampled_covariance)
separability_metrics <- data.table(
  sampled_coordinates = coordinate_count,
  training_visits = length(training_rows),
  training_subjects = length(training_subject_counts),
  empirical_weighting = "visit_pooled_to_match_identity_fit",
  empirical_model_offdiag_correlation = suppressWarnings(cor(
    empirical_sampled_covariance[off_diagonal],
    model_sampled_covariance[off_diagonal]
  )),
  relative_frobenius_error = sqrt(
    sum((empirical_sampled_covariance - model_sampled_covariance)^2) /
      sum(empirical_sampled_covariance^2)
  ),
  mean_temporal_diagonal = mean_temporal_variance
)
fwrite(
  separability_metrics,
  file.path(args$output_dir, "sampled_separability_diagnostic.csv")
)

# -------------------------------------------------------------------------
# Eigenspectra, scores, and held-out interpolation

spectra <- rbindlist(list(
  data.table(mode = "token", component = seq_along(eig_full$values$lam_phi),
             eigenvalue = eig_full$values$lam_phi),
  data.table(mode = "feature", component = seq_along(eig_full$values$lam_psi),
             eigenvalue = eig_full$values$lam_psi),
  data.table(mode = "time", component = seq_along(eig_full$values$lam_theta),
             eigenvalue = eig_full$values$lam_theta)
))
spectra[, positive_eigenvalue := pmax(eigenvalue, 0)]
spectra[, marginal_pve := positive_eigenvalue / sum(positive_eigenvalue), by = mode]
spectra[, cumulative_pve := cumsum(marginal_pve), by = mode]
fwrite(spectra, file.path(args$output_dir, "marginal_eigenspectra.csv"))
feature_pca_spectrum <- data.table(
  component = seq_along(feature_projection$eigenvalues),
  eigenvalue = feature_projection$eigenvalues
)
feature_pca_spectrum[, positive_eigenvalue := pmax(eigenvalue, 0)]
feature_pca_spectrum[, cumulative_pve :=
  cumsum(positive_eigenvalue) / sum(positive_eigenvalue)]
fwrite(feature_pca_spectrum, file.path(args$output_dir, "feature_pca_spectrum.csv"))

effective_ranks <- data.table(
  mode = c("raw_encoder_feature_pca", "token", "compressed_feature", "time"),
  nominal_rank = c(
    length(feature_projection$eigenvalues), nrow(covariance$U),
    nrow(covariance$V), nrow(covariance$K)
  ),
  effective_rank = c(
    effective_rank(feature_projection$eigenvalues),
    effective_rank(eig_full$values$lam_phi),
    effective_rank(eig_full$values$lam_psi),
    effective_rank(eig_full$values$lam_theta)
  ),
  retained_rank = c(
    feature_projection$retained_rank, ranks[["row"]], ranks[["column"]],
    ranks[["time"]]
  )
)
fwrite(effective_ranks, file.path(args$output_dir, "effective_ranks.csv"))

png(file.path(args$output_dir, "marginal_eigenspectra.png"),
    width = 1500, height = 500, res = 150)
par(mfrow = c(1L, 3L), mar = c(4, 4, 2.2, 1))
for (mode_name in c("token", "feature", "time")) {
  values <- spectra[mode == mode_name]$positive_eigenvalue
  display <- seq_len(min(length(values), 50L))
  plot(display, values[display] / sum(values), type = "b", pch = 16,
       xlab = "Component", ylab = "Marginal PVE", main = paste(mode_name, "mode"))
  abline(v = ranks[[switch(mode_name, token = "row", feature = "column", time = "time")]],
         col = "#B2182B", lty = 2, lwd = 2)
}
dev.off()

temporal_display_count <- min(4L, ncol(eig_full$funcs$Theta))
temporal_table <- data.table(time_year = t_grid)
for (component in seq_len(temporal_display_count)) {
  temporal_table[[paste0("theta", component)]] <- eig_full$funcs$Theta[, component]
}
fwrite(temporal_table, file.path(args$output_dir, "temporal_eigenfunctions.csv"))
png(file.path(args$output_dir, "temporal_covariance_and_eigenfunctions.png"),
    width = 1400, height = 650, res = 150)
par(mfrow = c(1L, 2L), mar = c(4.2, 4.4, 2.4, 1.2))
image(t_grid, t_grid, covariance$K, col = hcl.colors(80, "Blue-Red 3"),
      xlab = "Years", ylab = "Years", main = "Estimated temporal covariance")
matplot(t_grid, eig_full$funcs$Theta[, seq_len(temporal_display_count), drop = FALSE],
        type = "l", lty = 1, lwd = 2, xlab = "Years", ylab = "Eigenfunction",
        main = "Quadrature-weighted temporal modes")
legend("topright", legend = paste("Mode", seq_len(temporal_display_count)),
       col = seq_len(temporal_display_count), lty = 1, bty = "n")
dev.off()

product <- ordered_product_scores(scores, eig)
product_index <- product$index
score_matrix <- product$scores
colnames(score_matrix) <- sprintf("FPC%03d", seq_len(ncol(score_matrix)))
fwrite(product_index, file.path(args$output_dir, "product_fpc_eigenvalues.csv"))

subject_metadata_columns <- intersect(
  c(
    "baseline_age", "baseline_dx", "dx_at_or_near_mri", "outcome_type",
    "incident_ad_conversion", "prevalent_ad_at_baseline", "sex",
    "baseline_mmse", "baseline_cdrsb", "apoe_e4_count"
  ),
  names(metadata)
)
first_nonmissing <- function(values) {
  available <- which(!is.na(values) & as.character(values) != "")
  if (length(available)) values[[available[[1L]]]] else values[[1L]]
}
subject_table <- metadata[, lapply(.SD, first_nonmissing), by = subject_id,
                          .SDcols = subject_metadata_columns]
score_table <- cbind(
  subject_table[match(subject_ids, subject_id)],
  as.data.table(score_matrix)
)

score_input_count <- min(args$top_scores, ncol(score_matrix), length(subject_ids) - 2L)
valid_score_columns <- which(vapply(
  seq_len(score_input_count),
  function(column) all(is.finite(score_matrix[, column])) &&
    stats::sd(score_matrix[, column]) > 1e-10,
  logical(1L)
))
if (length(valid_score_columns) >= 2L) {
  score_pca <- prcomp(
    score_matrix[, valid_score_columns, drop = FALSE],
    center = TRUE, scale. = FALSE, rank. = 2L
  )
  score_table[, `:=`(score_PCA1 = score_pca$x[, 1L], score_PCA2 = score_pca$x[, 2L])]
} else {
  score_table[, `:=`(score_PCA1 = score_matrix[, 1L], score_PCA2 = 0)]
}
fwrite(score_table, file.path(args$output_dir, "subject_fpc_scores.csv"))

group_column <- if ("baseline_dx" %in% names(score_table)) "baseline_dx" else NULL
group <- if (is.null(group_column)) rep("all", nrow(score_table)) else
  as.character(score_table[[group_column]])
group[is.na(group) | !nzchar(group)] <- "unknown"
group_colors <- setNames(hcl.colors(length(unique(group)), "Dark 3"), sort(unique(group)))
png(file.path(args$output_dir, "subject_score_pca.png"),
    width = 850, height = 700, res = 150)
plot(score_table$score_PCA1, score_table$score_PCA2, pch = 19,
     col = group_colors[group], xlab = "PCA 1 of leading MatFPCA scores",
     ylab = "PCA 2 of leading MatFPCA scores", main = "Subject score geometry")
legend("topright", legend = names(group_colors), col = group_colors,
       pch = 19, bty = "n", cex = 0.85)
dev.off()

association_rows <- list()
association_position <- 0L
association_components <- seq_len(min(10L, ncol(score_matrix)))
if ("baseline_age" %in% names(score_table)) {
  age <- as.numeric(score_table$baseline_age)
  for (component in association_components) {
    complete <- is.finite(score_matrix[, component]) & is.finite(age)
    if (sum(complete) < 4L || length(unique(age[complete])) < 2L) next
    test <- suppressWarnings(cor.test(
      score_matrix[complete, component], age[complete],
      method = "spearman", exact = FALSE
    ))
    association_position <- association_position + 1L
    association_rows[[association_position]] <- data.table(
      product_pc = component, phenotype = "baseline_age", method = "Spearman",
      estimate = unname(test$estimate), p_value = test$p.value
    )
  }
}
if ("incident_ad_conversion" %in% names(score_table)) {
  conversion_text <- tolower(as.character(score_table$incident_ad_conversion))
  conversion <- rep(NA, length(conversion_text))
  conversion[conversion_text %in% c("true", "1", "yes")] <- TRUE
  conversion[conversion_text %in% c("false", "0", "no")] <- FALSE
  complete <- !is.na(conversion)
  if (sum(complete) >= 4L && length(unique(conversion[complete])) == 2L) {
    for (component in association_components) {
      fit <- summary(lm(score_matrix[complete, component] ~ conversion[complete]))
      association_position <- association_position + 1L
      association_rows[[association_position]] <- data.table(
        product_pc = component, phenotype = "incident_ad_conversion",
        method = "unadjusted mean difference", estimate = coef(fit)[2L, 1L],
        p_value = coef(fit)[2L, 4L]
      )
    }
  }
}
association_table <- if (length(association_rows)) rbindlist(association_rows) else
  data.table(product_pc = integer(), phenotype = character(), method = character(),
             estimate = numeric(), p_value = numeric())
association_table[, q_value := p.adjust(p_value, method = "BH"), by = phenotype]
fwrite(association_table, file.path(args$output_dir, "exploratory_fpc_associations.csv"))

reconstruction_summary <- data.table(split = character())
if (args$reconstruction_analysis) {
log_message("Evaluating training reconstruction and held-out interior-visit interpolation...")
metric_rows <- vector("list", nrow(metadata))
for (row in seq_len(nrow(metadata))) {
  observed <- load_projected_visit(row)
  reconstructed <- reconstruct_subject_latent(
    model, metadata$subject_id[[row]], metadata$time_year[[row]]
  )
  mean_at_visit <- predict_latent_mean(mean_model, metadata$time_year[[row]])
  error <- reconstructed - observed
  centered_observed <- observed - mean_at_visit
  denominator <- sum(centered_observed^2)
  vector_observed <- c(observed)
  vector_reconstructed <- c(reconstructed)
  metric_rows[[row]] <- data.table(
    analysis_row = row,
    subject_id = metadata$subject_id[[row]],
    visit_id = metadata$visit_id[[row]],
    time_year = metadata$time_year[[row]],
    split = metadata$split[[row]],
    rmse = sqrt(mean(error^2)),
    relative_rmse = if (denominator > 0) sqrt(sum(error^2) / denominator) else NA_real_,
    correlation = suppressWarnings(cor(vector_observed, vector_reconstructed)),
    cosine_similarity = sum(vector_observed * vector_reconstructed) /
      sqrt(sum(vector_observed^2) * sum(vector_reconstructed^2))
  )
}
reconstruction_metrics <- rbindlist(metric_rows)
reconstruction_summary <- reconstruction_metrics[, .(
  visits = .N,
  subjects = uniqueN(subject_id),
  mean_rmse = mean(rmse),
  median_relative_rmse = median(relative_rmse, na.rm = TRUE),
  mean_correlation = mean(correlation, na.rm = TRUE),
  mean_cosine_similarity = mean(cosine_similarity, na.rm = TRUE)
), by = split]
fwrite(reconstruction_metrics, file.path(args$output_dir, "reconstruction_metrics.csv"))
fwrite(reconstruction_summary, file.path(args$output_dir, "reconstruction_summary.csv"))

# -------------------------------------------------------------------------
# Complete-grid latent reconstructions and encoder-only image proxies

heldout_subjects <- unique(metadata[split == "validation_test"]$subject_id)
remaining_subjects <- setdiff(subject_ids, heldout_subjects)
selected_reconstruction_subjects <- head(
  c(heldout_subjects, remaining_subjects), args$reconstructed_subjects
)
reconstructed_dir <- file.path(args$output_dir, "reconstructed_latents")
dir.create(reconstructed_dir, recursive = TRUE, showWarnings = FALSE)

regular_grid_rows <- list()
change_map_rows <- list()
change_maps <- list()
change_position <- 0L
for (subject_position in seq_along(subject_ids)) {
  subject <- subject_ids[[subject_position]]
  subject_metadata <- metadata[subject_id == subject]
  observed_range <- range(subject_metadata$time_year)
  # The regular grid is common from baseline onward. Preserve screening-stage
  # information only at this subject's actually observed negative times,
  # avoiding extrapolation of every subject to another subject's SC date.
  subject_prediction_years <- sort(unique(c(
    prediction_years,
    subject_metadata[time_tick < 0L]$time_year
  )))
  reconstructed_list <- reconstruct_subject_latent(
    model, subject, subject_prediction_years
  )
  reconstructed_array <- array(
    unlist(reconstructed_list, use.names = FALSE),
    dim = c(model$representation_shape, length(subject_prediction_years))
  )
  within_range <- subject_prediction_years >= observed_range[[1L]] &
    subject_prediction_years <= observed_range[[2L]]
  saveRDS(
    list(
      subject_id = subject,
      time_year = subject_prediction_years,
      grid_point_type = ifelse(
        subject_prediction_years < 0,
        "observed_screening_time",
        "common_regular_grid"
      ),
      within_observed_range = within_range,
      compressed_latent = reconstructed_array,
      feature_projection_reference = "matfpca_model.rds$feature_projection"
    ),
    file.path(reconstructed_dir, paste0(subject, ".rds")),
    compress = "xz"
  )
  regular_grid_rows[[subject_position]] <- data.table(
    subject_id = subject,
    time_year = subject_prediction_years,
    grid_point_type = ifelse(
      subject_prediction_years < 0,
      "observed_screening_time",
      "common_regular_grid"
    ),
    within_observed_range = within_range,
    latent_rms = vapply(
      reconstructed_list, function(value) sqrt(mean(value^2)), numeric(1L)
    )
  )

  if (!subject %in% selected_reconstruction_subjects) next

  # Summarize spatially localized change on the native 6 x 6 x 6 token grid.
  change_position <- change_position + 1L
  endpoints <- observed_range
  endpoint_latents <- reconstruct_subject_latent(model, subject, endpoints)
  token_change <- sqrt(rowMeans((endpoint_latents[[2L]] - endpoint_latents[[1L]])^2))
  token_zero <- seq_len(length(token_change)) - 1L
  change_map_rows[[change_position]] <- data.table(
    subject_id = subject,
    from_year = endpoints[[1L]],
    to_year = endpoints[[2L]],
    token = token_zero + 1L,
    depth = token_zero %/% 36L + 1L,
    height = (token_zero %% 36L) %/% 6L + 1L,
    width = token_zero %% 6L + 1L,
    compressed_latent_change_rms = token_change
  )
  change_maps[[change_position]] <- token_vector_to_cube(token_change)
}
regular_grid_table <- rbindlist(regular_grid_rows)
token_change_table <- rbindlist(change_map_rows)
fwrite(regular_grid_table, file.path(args$output_dir, "regular_grid_latent_summary.csv"))
fwrite(token_change_table, file.path(args$output_dir, "token_change_maps.csv"))

plot_cube_midplanes <- function(cube, title, palette = "Blue-Red 3") {
  stopifnot(length(dim(cube)) == 3L)
  planes <- list(
    sagittal = cube[ceiling(dim(cube)[[1L]] / 2), , ],
    coronal = cube[, ceiling(dim(cube)[[2L]] / 2), ],
    axial = cube[, , ceiling(dim(cube)[[3L]] / 2)]
  )
  limit <- max(abs(unlist(planes)), na.rm = TRUE)
  if (!is.finite(limit) || limit == 0) limit <- 1
  for (plane_name in names(planes)) {
    image_arguments <- list(
      x = t(planes[[plane_name]])[
        , nrow(planes[[plane_name]]):1L, drop = FALSE
      ],
      axes = FALSE, useRaster = TRUE, col = hcl.colors(80, palette),
      main = paste(title, plane_name)
    )
    if (palette == "Blue-Red 3") image_arguments$zlim <- c(-limit, limit)
    do.call(image, image_arguments)
    box()
  }
}

token_display_count <- min(3L, ncol(eig$funcs$Phi))
png(file.path(args$output_dir, "token_eigenmaps.png"),
    width = 1350, height = 390 * token_display_count, res = 150)
par(mfrow = c(token_display_count, 3L), mar = c(1, 1, 2.2, 1))
for (component in seq_len(token_display_count)) {
  plot_cube_midplanes(
    token_vector_to_cube(eig$funcs$Phi[, component]),
    paste("Token mode", component)
  )
}
dev.off()

if (length(change_maps)) {
  png(file.path(args$output_dir, "subject_token_change_maps.png"),
      width = 1350, height = 390 * length(change_maps), res = 150)
  par(mfrow = c(length(change_maps), 3L), mar = c(1, 1, 2.2, 1))
  for (position in seq_along(change_maps)) {
    plot_cube_midplanes(
      change_maps[[position]], selected_reconstruction_subjects[[position]],
      palette = "YlOrRd"
    )
  }
  dev.off()
}

# BrainIAC is an encoder, not a decoder.  The honest image-space proxy is the
# nearest observed scan to each reconstructed latent matrix.  Excluding the
# query subject prevents an observed scan from retrieving itself.
candidate_rows <- which(metadata$split == "population_train")
if (length(candidate_rows) > args$retrieval_candidates) {
  set.seed(args$seed + 5L)
  candidate_rows <- sort(sample(candidate_rows, args$retrieval_candidates))
}
candidate_matrix <- t(vapply(
  candidate_rows, function(row) c(load_projected_visit(row)),
  numeric(prod(model$representation_shape))
))
retrieval_rows <- list()
retrieval_position <- 0L
for (subject in selected_reconstruction_subjects) {
  subject_prediction_years <- sort(unique(c(
    prediction_years,
    metadata[subject_id == subject & time_tick < 0L]$time_year
  )))
  for (time in subject_prediction_years) {
    eligible <- metadata$subject_id[candidate_rows] != subject
    if (!any(eligible)) next
    query <- c(reconstruct_subject_latent(model, subject, time))
    differences <- sweep(candidate_matrix[eligible, , drop = FALSE], 2L, query, "-")
    distances <- sqrt(rowMeans(differences^2))
    eligible_rows <- candidate_rows[eligible]
    best_row <- eligible_rows[[which.min(distances)]]
    retrieval_position <- retrieval_position + 1L
    retrieval_rows[[retrieval_position]] <- data.table(
      query_subject_id = subject,
      query_time_year = time,
      query_within_observed_range = time >= min(metadata[subject_id == subject]$time_year) &
        time <= max(metadata[subject_id == subject]$time_year),
      retrieved_subject_id = metadata$subject_id[[best_row]],
      retrieved_visit_id = metadata$visit_id[[best_row]],
      retrieved_image_id = metadata$image_id[[best_row]],
      retrieved_time_year = metadata$time_year[[best_row]],
      latent_rmse = min(distances),
      retrieved_nifti_path = metadata$nifti_path[[best_row]]
    )
  }
}
retrieval_table <- if (length(retrieval_rows)) rbindlist(retrieval_rows) else
  data.table(
    query_subject_id = character(), query_time_year = numeric(),
    query_within_observed_range = logical(), retrieved_subject_id = character(),
    retrieved_visit_id = character(), retrieved_image_id = integer(),
    retrieved_time_year = numeric(), latent_rmse = numeric(),
    retrieved_nifti_path = character()
  )
fwrite(retrieval_table, file.path(args$output_dir, "nearest_observed_mri_proxies.csv"))

robust_scale_image <- function(image) {
  image_dimensions <- dim(image)
  finite <- is.finite(image)
  foreground <- finite & image != 0
  reference <- image[foreground]
  if (length(reference) < 2L) reference <- image[finite]
  limits <- as.numeric(quantile(reference, c(0.01, 0.99), na.rm = TRUE))
  if (!all(is.finite(limits)) || limits[[2L]] <= limits[[1L]]) {
    limits <- range(reference, finite = TRUE)
  }
  scaled <- pmin(1, pmax(0, (image - limits[[1L]]) / diff(limits)))
  matrix(scaled, nrow = image_dimensions[[1L]], ncol = image_dimensions[[2L]])
}

if (nrow(retrieval_table) && requireNamespace("RNifti", quietly = TRUE)) {
  montage_rows <- retrieval_table[
    query_within_observed_range == TRUE,
    head(.SD, 3L), by = query_subject_id
  ]
  if (nrow(montage_rows)) {
    png(file.path(args$output_dir, "nearest_observed_mri_proxies.png"),
        width = 1050, height = 360 * ceiling(nrow(montage_rows) / 3), res = 150)
    par(mfrow = c(ceiling(nrow(montage_rows) / 3), 3L), mar = c(0, 0, 3.2, 0))
    for (row in seq_len(nrow(montage_rows))) {
      volume <- as.array(RNifti::readNifti(montage_rows$retrieved_nifti_path[[row]]))
      coronal_index <- min(110L, dim(volume)[[2L]])
      coronal <- volume[, coronal_index, , drop = TRUE]
      # Match the corrected coronal PNG convention: only a vertical display
      # flip relative to the existing exported files.
      display <- t(coronal[nrow(coronal):1L, ncol(coronal):1L, drop = FALSE])
      display <- robust_scale_image(display)
      plot.new()
      rasterImage(as.raster(display), 0, 0, 1, 1, interpolate = TRUE)
      title(main = sprintf(
        "%s at %.1fy -> %s at %.1fy\nnearest observed proxy (I%s)",
        montage_rows$query_subject_id[[row]], montage_rows$query_time_year[[row]],
        montage_rows$retrieved_subject_id[[row]], montage_rows$retrieved_time_year[[row]],
        montage_rows$retrieved_image_id[[row]]
      ), cex.main = 0.72)
    }
    dev.off()
  }
}

png(file.path(args$output_dir, "reconstruction_error.png"),
    width = 850, height = 650, res = 150)
reconstruction_metrics[, split_display := factor(
  split,
  levels = c("population_train", "validation_conditioning", "validation_test")
)]
boxplot(
  relative_rmse ~ split_display, data = reconstruction_metrics,
  outline = FALSE, col = c("#92C5DE", "#D9D9D9", "#F4A582"),
  xlab = "Visit split", ylab = "Relative latent RMSE",
  main = "Observed-visit latent reconstruction"
)
stripchart(
  relative_rmse ~ split_display, data = reconstruction_metrics,
  vertical = TRUE, method = "jitter", pch = 16,
  col = grDevices::adjustcolor("#333333", 0.35), add = TRUE
)
dev.off()
reconstruction_metrics[, split_display := NULL]
} else {
  log_message(
    "Skipping reconstruction analysis; interpretation is handled by ",
    "interpret_brainiac_matfpca.R."
  )
}

saveRDS(model, file.path(args$output_dir, "matfpca_model.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(args$output_dir, "session_info.txt"))
run_summary <- list(
  status = "complete",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  subjects = length(subject_ids),
  visits = nrow(metadata),
  population_training_subjects = length(population_subject_ids),
  population_training_visits = sum(metadata$split == "population_train"),
  validation_subjects = length(validation_subjects),
  validation_conditioning_visits = sum(metadata$split == "validation_conditioning"),
  heldout_interior_visits = sum(metadata$split == "validation_test"),
  original_latent_shape = unname(model$original_representation_shape),
  compressed_latent_shape = unname(model$representation_shape),
  selected_matfpca_ranks = unname(ranks),
  feature_pca_retained_pve = feature_projection$retained_pve,
  covariance_estimand = "visit/pair pooled identity MatFPCA estimator",
  common_prediction_grid_years = prediction_years,
  screening_predictions = paste(
    "Each subject includes only that subject's observed negative-time",
    "screening visits in addition to the common baseline-to-max grid."
  ),
  elapsed_matfpca_seconds = analysis_fit$elapsed_seconds,
  holdout_metrics = if ("validation_test" %in% reconstruction_summary$split) {
    as.list(reconstruction_summary[split == "validation_test", -"split"])
  } else NULL,
  reconstruction_analysis = args$reconstruction_analysis,
  scientific_focus = paste(
    "Product-FPC morphology and phenotype interpretation; BrainIAC image",
    "recovery is not an estimand."
  ),
  configuration_signature = configuration_signature
)
jsonlite::write_json(
  run_summary, run_summary_path,
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null"
)
log_message(
  "Complete. Results: ", args$output_dir,
  "; held-out visits: ", sum(metadata$split == "validation_test"),
  "; retained ranks: ", paste(ranks, collapse = " x "), "."
)
