#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(png)
  library(rTensor)
})
if (!nzchar(Sys.getenv("TZ"))) Sys.setenv(TZ = "America/Vancouver")

required_packages <- c(
  "mgcv", "pracma", "foreach", "splines", "uwot", "Rtsne"
)
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
  normalizePath(file.path("application", "adni", "adni_matfpca.R"))
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
setwd(project_root)
source(file.path("R", "matfpca.R"))
source(file.path("application", "adni", "adni_continuum_utils.R"))

defaults <- list(
  manifest = "",
  image_dir = file.path("data", "ADNI", "affine_axial_slice094"),
  output_dir = file.path("application", "adni", "results"),
  seed = 20260817L,
  max_month = 48,
  time_step_months = 12 / 365.25,
  holdout_fraction = 0.20,
  pool = 2L,
  mask_intensity = 1 / 255,
  mask_min_fraction = 0.50,
  mask_dilation = 0L,
  mean_df = 5L,
  weight_scheme = "gaussian_refined",
  weight_folds = 3L,
  opw_iter = 1L,
  nbasis = "theory",
  mode_pve = 0.95,
  score_pcs = 30L,
  continuum_df = 5,
  continuum_max_iterations = 20L,
  continuum_landmarks = 6L,
  umap_neighbors = 30L,
  umap_min_dist = 0.20,
  tsne_perplexity = 30,
  top_fpcs = 4L,
  max_subjects = Inf,
  cache = TRUE
)

usage <- function() {
  cat(
    "Usage: Rscript application/adni/adni_matfpca.R [--option=value ...]\n\n",
    "Important options:\n",
    "  --output-dir=PATH           Result directory\n",
    "  --pool=2                    Block-average factor; gives 114 x 96 images\n",
    "  --mask-min-fraction=0.50    Training-image support needed for brain mask\n",
    "  --weight-scheme=SCHEME      gaussian_refined (default) or identity\n",
    "  --weight-folds=3            Cross-fitting folds for refined weights\n",
    "  --holdout-fraction=0.20     Random fraction of subjects tested\n",
    "  --time-step-months=VALUE    Visit-time resolution (default: one day)\n",
    "  --score-pcs=30              Leading raw FPC scores used for the continuum\n",
    "  --continuum-df=5            Principal-curve spline degrees of freedom\n",
    "  --continuum-landmarks=6     MRI positions shown along the continuum\n",
    "  --umap-neighbors=30         UMAP neighborhood size\n",
    "  --tsne-perplexity=30        t-SNE perplexity (reduced if needed)\n",
    "  --max-subjects=Inf          Random subject cap for smoke tests\n",
    "  --cache=true                Reuse image and fitted-model caches\n",
    sep = ""
  )
}

coerce_option <- function(value, template) {
  if (is.logical(template)) {
    value_lower <- tolower(value)
    if (!value_lower %in% c("true", "false", "1", "0", "yes", "no")) {
      stop("Logical option must be true or false, not: ", value)
    }
    return(value_lower %in% c("true", "1", "yes"))
  }
  if (is.integer(template)) return(as.integer(value))
  if (is.numeric(template)) return(as.numeric(value))
  value
}

parse_options <- function(default_values) {
  args <- commandArgs(trailingOnly = TRUE)
  if (any(args %in% c("--help", "-h"))) {
    usage()
    quit(save = "no", status = 0)
  }
  result <- default_values
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) {
      stop("Options must use --name=value syntax: ", arg)
    }
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", pair[[1L]], fixed = TRUE)
    value <- paste(pair[-1L], collapse = "=")
    if (!key %in% names(result)) stop("Unknown option: --", pair[[1L]])
    result[[key]] <- coerce_option(value, default_values[[key]])
  }
  result
}

args <- parse_options(defaults)
if (!args$weight_scheme %in% c("identity", "gaussian_refined")) {
  stop("--weight-scheme must be identity or gaussian_refined.")
}
if (!is.finite(args$max_month) || args$max_month <= 0) stop("Invalid max month.")
if (!is.finite(args$time_step_months) || args$time_step_months <= 0) {
  stop("--time-step-months must be positive.")
}
if (abs(args$max_month / args$time_step_months -
        round(args$max_month / args$time_step_months)) > 1e-8) {
  stop("--time-step-months must divide --max-month exactly.")
}
if (args$holdout_fraction <= 0 || args$holdout_fraction >= 1) {
  stop("--holdout-fraction must lie in (0, 1).")
}
if (args$pool < 1L || args$mean_df < 3L || args$score_pcs < 2L) {
  stop("pool >= 1, mean-df >= 3, and score-pcs >= 2 are required.")
}
if (args$top_fpcs < 1L) stop("--top-fpcs must be positive.")
if (args$mask_intensity < 0 || args$mask_min_fraction <= 0 ||
    args$mask_min_fraction >= 1 || args$mask_dilation < 0L) {
  stop("Invalid brain-mask settings.")
}
if (args$umap_neighbors < 2L || args$continuum_df < 2 ||
    args$continuum_max_iterations < 1L ||
    args$continuum_landmarks < 3L || args$umap_min_dist < 0 ||
    args$tsne_perplexity < 1) {
  stop("Invalid continuum or embedding settings.")
}

if (!nzchar(args$manifest)) {
  manifests <- sort(Sys.glob(file.path("data", "ADNI", "filtered_manifest_*.csv")))
  if (!length(manifests)) stop("No data/ADNI/filtered_manifest_*.csv was found.")
  args$manifest <- tail(manifests, 1L)
}
args$manifest <- normalizePath(args$manifest)
args$image_dir <- normalizePath(args$image_dir)
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
args$output_dir <- normalizePath(args$output_dir)
cache_dir <- file.path(args$output_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Remove only exact generated artifacts from superseded label- and
# cluster-based analyses so stale results cannot be mistaken for current output.
legacy_outputs <- file.path(args$output_dir, c(
  "fpc_score_clustering.png", "fpc_scores_and_clusters.csv",
  "clustering_metrics.csv", "clustering_contingency.csv",
  "fpc_score_umap_by_patient_group.png", "fpc_scores_and_umap.csv",
  "reconstruction_examples_by_group.png", "umap_representation_diagnostics.csv",
  "umap_neighbor_diagnostics.csv", "knn_probe_group_recall.csv",
  "knn_probe_confusion.csv", "top_fpc_spatial_panels.png",
  "top_fpc_temporal_panels.png", "fpc_scores_clusters_embeddings.csv",
  "clustering_candidate_metrics.csv", "cluster_summary.csv",
  "cluster_representative_subjects.csv", "fpc_score_umap_by_cluster.png",
  "fpc_score_embedding_comparison_by_cluster.png",
  "representative_mris_by_cluster.png", "age_summary_by_cluster.csv",
  "age_distribution_by_cluster.png", "age_cluster_association.csv",
  "age_pairwise_cluster_tests.csv", "fpc_scores_embeddings.csv",
  "embedding_quality_metrics.csv", "latent_coordinate_diagnostics.csv",
  "latent_coordinate_landmarks.csv", "latent_principal_curve_scores.csv",
  "fpc_score_embedding_comparison_by_continuum.png",
  "fpc_score_umap_by_continuum.png", "reconstructed_mri_continuum.png",
  "smoothed_observed_mri_continuum.png", "observed_mri_continuum.png",
  "latent_coordinate_vs_age.png", "adni_matfpca_model.rds",
  "run_summary.txt", "session_info.txt"
))
invisible(file.remove(legacy_outputs[file.exists(legacy_outputs)]))

log_message <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
  flush.console()
}

pool_image <- function(path, factor) {
  image <- png::readPNG(path)
  if (length(dim(image)) == 3L) image <- rowMeans(image, dims = 2L)
  if (!is.matrix(image)) stop("Not a grayscale-compatible PNG: ", path)
  if (factor == 1L) return(image)

  retained_rows <- floor(nrow(image) / factor) * factor
  retained_cols <- floor(ncol(image) / factor) * factor
  image <- image[seq_len(retained_rows), seq_len(retained_cols), drop = FALSE]
  row_groups <- rep(seq_len(retained_rows / factor), each = factor)
  image <- rowsum(image, row_groups, reorder = FALSE) / factor
  col_groups <- rep(seq_len(retained_cols / factor), each = factor)
  t(rowsum(t(image), col_groups, reorder = FALSE) / factor)
}

flood_component <- function(target, seeds) {
  nr <- nrow(target)
  nc <- ncol(target)
  visited <- matrix(FALSE, nr, nc)
  seeds <- unique(as.integer(seeds[target[seeds]]))
  if (!length(seeds)) return(visited)
  queue <- integer(length(target))
  queue[seq_along(seeds)] <- seeds
  visited[seeds] <- TRUE
  head <- 1L
  tail <- length(seeds)
  while (head <= tail) {
    index <- queue[[head]]
    head <- head + 1L
    row <- ((index - 1L) %% nr) + 1L
    column <- ((index - 1L) %/% nr) + 1L
    neighbors <- c(
      if (row > 1L) index - 1L,
      if (row < nr) index + 1L,
      if (column > 1L) index - nr,
      if (column < nc) index + nr
    )
    new_neighbors <- neighbors[target[neighbors] & !visited[neighbors]]
    if (length(new_neighbors)) {
      visited[new_neighbors] <- TRUE
      positions <- tail + seq_along(new_neighbors)
      queue[positions] <- new_neighbors
      tail <- tail + length(new_neighbors)
    }
  }
  visited
}

largest_connected_component <- function(mask) {
  remaining <- mask
  largest <- matrix(FALSE, nrow(mask), ncol(mask))
  while (any(remaining)) {
    seed <- which(remaining)[[1L]]
    component <- flood_component(remaining, seed)
    if (sum(component) > sum(largest)) largest <- component
    remaining[component] <- FALSE
  }
  largest
}

fill_binary_holes <- function(mask) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  boundary <- unique(c(
    seq_len(nr),
    seq.int(1L, by = nr, length.out = nc),
    seq.int(nr, by = nr, length.out = nc),
    (nc - 1L) * nr + seq_len(nr)
  ))
  outside <- flood_component(!mask, boundary)
  mask | (!mask & !outside)
}

dilate_binary_mask <- function(mask, iterations = 1L) {
  if (iterations <= 0L) return(mask)
  for (iteration in seq_len(iterations)) {
    padded <- matrix(FALSE, nrow(mask) + 2L, ncol(mask) + 2L)
    padded[2:(nrow(mask) + 1L), 2:(ncol(mask) + 1L)] <- mask
    expanded <- matrix(FALSE, nrow(mask), ncol(mask))
    for (row_shift in 0:2) {
      for (column_shift in 0:2) {
        expanded <- expanded | padded[
          1:nrow(mask) + row_shift,
          1:ncol(mask) + column_shift,
          drop = FALSE
        ]
      }
    }
    mask <- expanded
  }
  mask
}

infer_brain_mask <- function(
  images, image_dim, intensity_threshold, minimum_fraction, dilation
) {
  foreground_frequency <- matrix(
    colMeans(images > intensity_threshold), nrow = image_dim[[1L]]
  )
  support <- foreground_frequency >= minimum_fraction
  support <- largest_connected_component(support)
  support <- fill_binary_holes(support)
  support <- dilate_binary_mask(support, dilation)
  coverage <- mean(support)
  if (!is.finite(coverage) || coverage < 0.10 || coverage > 0.90) {
    stop(sprintf(
      "Inferred brain-mask coverage %.3f is implausible; adjust mask settings.",
      coverage
    ))
  }
  list(mask = support, foreground_frequency = foreground_frequency)
}

random_subject_sample <- function(subject_ids, size, seed) {
  if (length(subject_ids) < 8L) stop("At least eight eligible subjects are required.")
  if (!is.finite(size) || size >= length(subject_ids)) return(subject_ids)
  size <- min(length(subject_ids), max(8L, as.integer(size)))
  set.seed(seed)
  sample(subject_ids, size)
}

prepare_manifest <- function() {
  manifest <- data.table::fread(args$manifest, na.strings = c("", "NA"))
  required <- c(
    "PTID", "image_id", "image_time_month", "mri_VISCODE2",
    "image_retained_for_analysis", "baseline_age"
  )
  absent <- setdiff(required, names(manifest))
  if (length(absent)) stop("Manifest is missing: ", paste(absent, collapse = ", "))

  files <- Sys.glob(file.path(args$image_dir, "*.png"))
  image_ids <- suppressWarnings(as.integer(sub(
    ".*_I([0-9]+)_Warped\\.png$", "\\1", basename(files)
  )))
  if (anyNA(image_ids) || anyDuplicated(image_ids)) {
    stop("PNG filenames do not provide unique _I<image_id>_ identifiers.")
  }
  manifest[, image_file := files[match(as.integer(image_id), image_ids)]]
  manifest[, exclusion := fifelse(
    is.na(image_file), "missing_png",
    fifelse(
      image_retained_for_analysis != TRUE, "manifest_not_retained",
      fifelse(
        !is.finite(image_time_month), "missing_time",
        fifelse(image_time_month > args$max_month, "after_m48", "candidate")
      )
    )
  )]

  candidates <- manifest[exclusion == "candidate"]
  candidates[, available_visits := .N, by = PTID]
  candidates[available_visits < 2L, exclusion := "subject_fewer_than_2_images"]
  ineligible_subjects <- unique(candidates[
    exclusion == "subject_fewer_than_2_images", PTID
  ])
  manifest[
    exclusion == "candidate" & PTID %in% ineligible_subjects,
    exclusion := "subject_fewer_than_2_images"
  ]
  retained <- candidates[exclusion == "candidate"]
  if (!nrow(retained)) stop("No eligible images remain after filtering.")

  subject_ages <- unique(retained[, .(PTID, baseline_age)])
  if (anyDuplicated(subject_ages$PTID)) {
    stop("A subject has inconsistent baseline ages in the manifest.")
  }
  selected_subjects <- random_subject_sample(
    subject_ages$PTID, args$max_subjects, args$seed + 10L
  )
  retained <- retained[PTID %in% selected_subjects]

  retained[, time_tick := as.integer(round(image_time_month / args$time_step_months))]
  # The clinical visit code, rather than a small acquisition-date discrepancy,
  # defines baseline for this analysis.
  retained[mri_VISCODE2 == "bl", time_tick := 0L]
  duplicate_ticks <- retained[, .N, by = .(PTID, time_tick)][N > 1L]
  if (nrow(duplicate_ticks)) {
    stop(
      nrow(duplicate_ticks),
      " subject/time collisions after snapping. Use a smaller --time-step-months."
    )
  }
  setorder(retained, PTID, time_tick, image_id)
  retained[, analysis_row := .I]

  exclusion_counts <- manifest[, .N, by = exclusion]
  exclusion_counts[
    exclusion == "candidate",
    exclusion := "eligible_before_optional_subject_cap"
  ]
  exclusion_counts <- rbind(
    exclusion_counts,
    data.table(exclusion = "analysis_retained", N = nrow(retained))
  )
  list(data = retained, exclusions = exclusion_counts, original = manifest)
}

prepared <- prepare_manifest()
meta <- prepared$data

manifest_hash <- unname(tools::md5sum(args$manifest))
image_cache <- file.path(
  cache_dir,
  sprintf("images_%s_pool%02d.rds", substr(manifest_hash, 1L, 12L), args$pool)
)
cached_images <- NULL
if (args$cache && file.exists(image_cache)) {
  candidate <- readRDS(image_cache)
  if (identical(candidate$image_ids, meta$image_id) &&
      identical(candidate$image_files, meta$image_file)) {
    cached_images <- candidate
    log_message("Reusing pooled image cache: ", image_cache)
  }
}

if (is.null(cached_images)) {
  log_message("Reading and block-averaging ", nrow(meta), " PNG slices...")
  first_image <- pool_image(meta$image_file[[1L]], args$pool)
  image_dim <- dim(first_image)
  image_matrix <- matrix(NA_real_, nrow(meta), prod(image_dim))
  image_matrix[1L, ] <- as.vector(first_image)
  if (nrow(meta) > 1L) {
    for (i in 2:nrow(meta)) {
      current <- pool_image(meta$image_file[[i]], args$pool)
      if (!identical(dim(current), image_dim)) stop("Inconsistent PNG dimensions.")
      image_matrix[i, ] <- as.vector(current)
      if (i %% 250L == 0L) log_message("Read ", i, " / ", nrow(meta), " images")
    }
  }
  cached_images <- list(
    image_ids = meta$image_id,
    image_files = meta$image_file,
    image_dim = image_dim,
    image_matrix = image_matrix
  )
  if (args$cache) saveRDS(cached_images, image_cache, compress = FALSE)
}
image_dim <- cached_images$image_dim
image_matrix <- cached_images$image_matrix

# Select subjects for one-visit holdout without using any clinical outcome label.
subject_table <- unique(meta[, .(PTID, baseline_age)])
if (anyDuplicated(subject_table$PTID)) stop("A subject has inconsistent baseline ages.")
set.seed(args$seed)
holdout_count <- max(1L, as.integer(round(args$holdout_fraction * nrow(subject_table))))
holdout_subjects <- sample(subject_table$PTID, min(holdout_count, nrow(subject_table)))
meta[, split := "train"]
set.seed(args$seed + 1L)
test_rows <- meta[PTID %in% holdout_subjects, sample(.I, 1L), by = PTID]$V1
meta[test_rows, split := "test"]
if (any(meta[split == "train", .N, by = PTID]$N < 1L)) {
  stop("A held-out subject has no remaining training image.")
}

# Infer a common support from training images only. The largest connected
# foreground component is hole-filled to recover a robust common brain
# silhouette. Optional dilation is available, but is off by default.
train_rows <- which(meta$split == "train")
mask_fit <- infer_brain_mask(
  image_matrix[train_rows, , drop = FALSE],
  image_dim = image_dim,
  intensity_threshold = args$mask_intensity,
  minimum_fraction = args$mask_min_fraction,
  dilation = args$mask_dilation
)
brain_mask <- mask_fit$mask
mask_vector <- as.vector(brain_mask)
image_matrix[, !mask_vector] <- 0
log_message(
  "Inferred training-only brain mask covering ",
  sprintf("%.1f%%", 100 * mean(brain_mask)), " of image pixels"
)

# Use integer day-resolution ticks internally, but pass years to MatFPCA. The
# covariance grid needs only observed locations plus baseline and the m48 end;
# omitting unobserved daily locations keeps the temporal eigendecomposition
# small without changing any subject's visit time.
grid_min_tick <- min(meta$time_tick)
grid_max_tick <- as.integer(round(args$max_month / args$time_step_months))
grid_ticks <- sort(unique(c(meta$time_tick, 0L, grid_max_tick)))
t_grid <- grid_ticks * args$time_step_months / 12
meta[, time_year := time_tick * args$time_step_months / 12]
meta[, time_index := match(time_tick, grid_ticks)]
stopifnot(all(is.finite(meta$time_index)), any(t_grid < 0), any(t_grid == 0))

fit_mean_function <- function(times, values, subject_ids, df, boundary) {
  spline_basis <- splines::ns(
    times,
    df = df - 1L,
    intercept = FALSE,
    Boundary.knots = boundary
  )
  design <- cbind(intercept = 1, spline_basis)
  subject_count <- table(subject_ids)
  weights <- 1 / as.numeric(subject_count[subject_ids])
  normal_matrix <- crossprod(design, design * weights)
  ridge <- 1e-8 * max(mean(diag(normal_matrix)), 1)
  coefficients <- solve(
    normal_matrix + diag(ridge, ncol(design)),
    crossprod(design, values * weights)
  )
  list(
    coefficients = coefficients,
    knots = attr(spline_basis, "knots"),
    Boundary.knots = attr(spline_basis, "Boundary.knots")
  )
}

mean_design <- function(model, times) {
  cbind(
    intercept = 1,
    splines::ns(
      times,
      knots = model$knots,
      intercept = FALSE,
      Boundary.knots = model$Boundary.knots
    )
  )
}

log_message(
  "Fitting training-only temporal mean to ", length(train_rows),
  " images from ", uniqueN(meta$PTID), " subjects..."
)
mean_model <- fit_mean_function(
  meta$time_year[train_rows],
  image_matrix[train_rows, , drop = FALSE],
  meta$PTID[train_rows],
  args$mean_df,
  range(t_grid)
)
mean_all <- mean_design(mean_model, meta$time_year) %*% mean_model$coefficients
centered_images <- image_matrix - mean_all
rm(mean_all)

subject_ids <- sort(unique(meta$PTID))
train_by_subject <- split(train_rows, factor(meta$PTID[train_rows], levels = subject_ids))
X_list <- lapply(train_by_subject, function(rows) {
  subject_array <- array(NA_real_, dim = c(image_dim, length(rows)))
  for (j in seq_along(rows)) {
    subject_array[, , j] <- matrix(centered_images[rows[[j]], ], nrow = image_dim[[1L]])
  }
  rTensor::as.tensor(subject_array)
})
t_idx_list <- lapply(train_by_subject, function(rows) as.integer(meta$time_index[rows]))
mat_data <- create_matrix_funcdata_irreg(
  X_list = X_list,
  t_grid = t_grid,
  t_idx_list = t_idx_list
)

configuration_signature <- paste(
  "adni_pipeline_v3_masked_highres",
  manifest_hash,
  paste(meta$image_id, meta$time_tick, sep = ":", collapse = ","),
  paste(meta$image_id[meta$split == "test"], collapse = ","),
  args$pool, args$mask_intensity, args$mask_min_fraction, args$mask_dilation,
  args$time_step_months, args$mean_df, args$weight_scheme,
  args$weight_folds, args$opw_iter, args$nbasis, args$mode_pve,
  sep = "|"
)
fit_cache <- file.path(cache_dir, "matfpca_fit.rds")
analysis_fit <- NULL
if (args$cache && file.exists(fit_cache)) {
  candidate <- readRDS(fit_cache)
  if (identical(candidate$signature, configuration_signature)) {
    analysis_fit <- candidate
    log_message("Reusing MatFPCA fit cache: ", fit_cache)
  }
}

rank_for_pve <- function(values, pve) {
  values <- pmax(as.numeric(values), 0)
  if (!length(values) || sum(values) <= 0) stop("Non-positive marginal eigensystem.")
  which(cumsum(values) / sum(values) >= pve)[[1L]]
}

if (is.null(analysis_fit)) {
  fit_nbasis <- if (grepl("^[0-9]+$", args$nbasis)) {
    as.integer(args$nbasis)
  } else {
    args$nbasis
  }
  fold_id <- NULL
  if (args$weight_scheme == "gaussian_refined") {
    args$weight_folds <- min(args$weight_folds, length(subject_ids))
    set.seed(args$seed + 2L)
    fold_id <- sample(rep(seq_len(args$weight_folds), length.out = length(subject_ids)))
  }
  log_message(
    "Fitting MatFPCA (", args$weight_scheme, ", ",
    if (args$weight_scheme == "identity") 1L else args$weight_folds,
    " fold(s), image ", paste(image_dim, collapse = " x "),
    ", grid ", length(t_grid), ")..."
  )
  timing <- system.time({
    covariance_fit <- fit_matfcov(
      mat_data,
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
  eig_full <- eigen.matfcov(covariance_fit$C)
  ranks <- c(
    row = rank_for_pve(eig_full$values$lam_phi, args$mode_pve),
    column = rank_for_pve(eig_full$values$lam_psi, args$mode_pve),
    time = rank_for_pve(eig_full$values$lam_theta, args$mode_pve)
  )
  eig <- trunc.matfeigen(
    eig_full,
    I = ranks[["row"]], J = ranks[["column"]], K = ranks[["time"]],
    verbose = FALSE
  )
  log_message(
    "Estimating subject scores at retained ranks ",
    paste(ranks, collapse = " x "), "..."
  )
  scores <- get_matfpca_scores(
    mat_data,
    eig,
    sig2 = covariance_fit$C$sig2,
    K_grid = covariance_fit$C$K
  )
  analysis_fit <- list(
    signature = configuration_signature,
    covariance_fit = covariance_fit,
    eig = eig,
    ranks = ranks,
    scores = scores,
    elapsed_seconds = unname(timing[["elapsed"]])
  )
  if (args$cache) saveRDS(analysis_fit, fit_cache, compress = FALSE)
}

covariance_fit <- analysis_fit$covariance_fit
eig <- analysis_fit$eig
scores <- analysis_fit$scores
ranks <- analysis_fit$ranks

# The retained eigensystem can contain only one temporal mode when that mode
# clears the marginal PVE threshold. Recover the full spatial eigensystem for
# product-rank diagnostics. The temporal diagnostic below uses quadrature
# weights because the compact ADNI grid is irregular rather than equally spaced.
eig_full <- eigen.matfcov(covariance_fit$C)

trapezoid_weights <- function(grid) {
  if (length(grid) < 3L || is.unsorted(grid, strictly = TRUE)) {
    stop("Temporal diagnostics require a strictly increasing grid.")
  }
  c(
    (grid[[2L]] - grid[[1L]]) / 2,
    (grid[3:length(grid)] - grid[1:(length(grid) - 2L)]) / 2,
    (grid[[length(grid)]] - grid[[length(grid) - 1L]]) / 2
  )
}

temporal_weights <- trapezoid_weights(t_grid)
temporal_sqrt_weights <- sqrt(temporal_weights)
temporal_weighted_covariance <- sweep(
  sweep(covariance_fit$C$K, 1L, temporal_sqrt_weights, "*"),
  2L, temporal_sqrt_weights, "*"
)
temporal_integral_eigen <- eigen(
  0.5 * (temporal_weighted_covariance + t(temporal_weighted_covariance)),
  symmetric = TRUE
)
temporal_values <- pmax(as.numeric(temporal_integral_eigen$values), 0)
temporal_functions <- sweep(
  temporal_integral_eigen$vectors, 1L, temporal_sqrt_weights, "/"
)
temporal_positive_total <- sum(temporal_values)
temporal_relative_tolerance <- 1e-10
temporal_mode_count <- max(
  2L,
  sum(temporal_values >= max(temporal_values) * temporal_relative_tolerance)
)
temporal_mode_count <- min(temporal_mode_count, ncol(temporal_functions))
temporal_baseline_index <- which.min(abs(t_grid))

temporal_diagnostics <- rbindlist(lapply(seq_len(temporal_mode_count), function(k) {
  profile <- temporal_functions[, k]
  weighted_mean <- sum(temporal_weights * profile) / sum(temporal_weights)
  total_energy <- sum(temporal_weights * profile^2)
  nonconstant_energy <- sum(temporal_weights * (profile - weighted_mean)^2)
  rms <- sqrt(total_energy / sum(temporal_weights))
  data.table(
    time_component = k,
    eigenvalue = temporal_values[[k]],
    marginal_variance_fraction = temporal_values[[k]] / temporal_positive_total,
    cumulative_marginal_variance_fraction =
      sum(temporal_values[seq_len(k)]) / temporal_positive_total,
    nonconstant_fraction = nonconstant_energy / total_energy,
    range_over_rms = diff(range(profile)) / rms,
    baseline_to_end_change_over_rms =
      abs(profile[[length(profile)]] - profile[[temporal_baseline_index]]) / rms
  )
}))

# This screen separates a visibly dynamic, variance-bearing component from
# high-frequency numerical modes: at least 0.1% marginal PVE and at least half
# of the eigenfunction energy outside its best constant approximation.
temporal_diagnostics[, `:=`(
  variance_bearing = marginal_variance_fraction >= 0.001,
  highly_nonconstant = nonconstant_fraction >= 0.50
)]
temporal_diagnostics[, dynamic_screen := variance_bearing & highly_nonconstant]
temporal_diagnostics[, interpretation := fifelse(
  dynamic_screen,
  "material non-constant temporal mode",
  fifelse(
    variance_bearing,
    "material but nearly constant temporal mode",
    "non-constant shape with negligible variance"
  )
)]

predict_mean_matrix <- function(time_year) {
  vector <- as.numeric(mean_design(mean_model, time_year) %*% mean_model$coefficients)
  matrix(vector, nrow = image_dim[[1L]])
}

reconstruct_centered <- function(subject_position, time_indices) {
  score_array <- scores@data[subject_position, , , , drop = TRUE]
  dim(score_array) <- c(ranks[["row"]], ranks[["column"]], ranks[["time"]])
  theta <- eig$funcs$Theta[time_indices, , drop = FALSE]
  rTensor::ttl(
    rTensor::as.tensor(score_array),
    list(eig$funcs$Phi, eig$funcs$Psi, theta),
    ms = 1:3
  )@data
}

image_metrics <- function(observed, reconstructed, centered_observed, mask) {
  observed <- observed[mask]
  reconstructed <- reconstructed[mask]
  centered_observed <- centered_observed[mask]
  error <- reconstructed - observed
  rmse <- sqrt(mean(error^2))
  centered_norm <- sum(centered_observed^2)
  correlation <- suppressWarnings(cor(observed, reconstructed))
  c(
    rmse = rmse,
    relative_rmse = if (centered_norm > 0) sqrt(sum(error^2) / centered_norm) else NA_real_,
    correlation = correlation,
    psnr = if (rmse > 0) 20 * log10(1 / rmse) else Inf
  )
}

log_message("Evaluating training reconstructions and held-out predictions...")
metric_rows <- vector("list", nrow(meta))
reconstruction_store <- vector("list", length(subject_ids))
names(reconstruction_store) <- subject_ids
for (subject_position in seq_along(subject_ids)) {
  subject <- subject_ids[[subject_position]]
  rows <- which(meta$PTID == subject)
  reconstructed_centered <- reconstruct_centered(subject_position, meta$time_index[rows])
  if (length(rows) == 1L) dim(reconstructed_centered) <- c(image_dim, 1L)
  subject_reconstructions <- vector("list", length(rows))
  names(subject_reconstructions) <- as.character(rows)
  for (j in seq_along(rows)) {
    row <- rows[[j]]
    observed <- matrix(image_matrix[row, ], nrow = image_dim[[1L]])
    centered_observed <- matrix(centered_images[row, ], nrow = image_dim[[1L]])
    reconstructed <- (
      reconstructed_centered[, , j] + predict_mean_matrix(meta$time_year[[row]])
    ) * brain_mask
    values <- image_metrics(observed, reconstructed, centered_observed, brain_mask)
    metric_rows[[row]] <- data.table(
      analysis_row = row,
      PTID = subject,
      image_id = meta$image_id[[row]],
      split = meta$split[[row]],
      time_year = meta$time_year[[row]],
      mask_pixels = sum(brain_mask),
      rmse = values[["rmse"]],
      relative_rmse = values[["relative_rmse"]],
      correlation = values[["correlation"]],
      psnr = values[["psnr"]]
    )
    subject_reconstructions[[as.character(row)]] <- reconstructed
  }
  reconstruction_store[[subject]] <- subject_reconstructions
}
reconstruction_metrics <- rbindlist(metric_rows)
reconstruction_summary <- reconstruction_metrics[, .(
  images = .N,
  subjects = uniqueN(PTID),
  rmse_mean = mean(rmse),
  rmse_median = median(rmse),
  relative_rmse_mean = mean(relative_rmse, na.rm = TRUE),
  relative_rmse_median = median(relative_rmse, na.rm = TRUE),
  correlation_mean = mean(correlation, na.rm = TRUE),
  psnr_mean = mean(psnr)
), by = split]

# Order product FPCs by their product eigenvalue and expose their score columns.
product_values_array <- eig$values$lam_phi %o%
  eig$values$lam_psi %o% eig$values$lam_theta
product_index <- as.data.table(arrayInd(
  seq_along(product_values_array), dim(product_values_array)
))
setnames(product_index, c("row_component", "column_component", "time_component"))
product_index[, eigenvalue := as.vector(product_values_array)]
product_index <- product_index[is.finite(eigenvalue) & eigenvalue > 0][order(-eigenvalue)]
product_index[, product_pc := .I]
product_index[, pve := cumsum(eigenvalue) / sum(eigenvalue)]

# Locate the leading product FPC associated with each diagnosed time mode. This
# makes clear why a dynamic marginal mode may not appear among the leading
# product components even when its temporal shape is pronounced.
diagnostic_product_values <-
  pmax(eig_full$values$lam_phi, 0) %o%
  pmax(eig_full$values$lam_psi, 0) %o%
  temporal_values[seq_len(temporal_mode_count)]
diagnostic_product_order <- order(
  as.vector(diagnostic_product_values), decreasing = TRUE
)
diagnostic_product_indices <- arrayInd(
  diagnostic_product_order, dim(diagnostic_product_values)
)
temporal_diagnostics[, leading_full_product_rank := vapply(
  time_component,
  function(k) which(diagnostic_product_indices[, 3L] == k)[[1L]],
  integer(1L)
)]
temporal_diagnostics[, leading_product_eigenvalue := vapply(
  time_component,
  function(k) max(diagnostic_product_values[, , k]),
  numeric(1L)
)]

score_matrix_all <- matrix(scores@data, nrow = length(subject_ids))
score_columns <- with(product_index, {
  row_component +
    (column_component - 1L) * ranks[["row"]] +
    (time_component - 1L) * ranks[["row"]] * ranks[["column"]]
})
score_matrix_ordered <- score_matrix_all[, score_columns, drop = FALSE]
colnames(score_matrix_ordered) <- sprintf("FPC%03d", seq_len(ncol(score_matrix_ordered)))
score_table <- data.table(
  PTID = subject_ids,
  baseline_age = subject_table$baseline_age[match(subject_ids, subject_table$PTID)]
)
score_table <- cbind(score_table, as.data.table(score_matrix_ordered))

n_score_pcs <- min(
  args$score_pcs, ncol(score_matrix_ordered), nrow(score_matrix_ordered) - 2L
)
score_input <- score_matrix_ordered[, seq_len(n_score_pcs), drop = FALSE]
valid_dimensions <- apply(
  score_input, 2L,
  function(x) all(is.finite(x)) && stats::sd(x) > 1e-10
)
score_input_positions <- which(valid_dimensions)
score_input <- score_input[, valid_dimensions, drop = FALSE]
if (ncol(score_input) < 2L) stop("Too few non-degenerate FPC score dimensions.")
score_dimension <- ncol(score_input)
score_dimension_suffix <- sprintf("_pcs%d", score_dimension)
continuum_output_file <- function(stem, extension) {
  file.path(
    args$output_dir,
    paste0(stem, score_dimension_suffix, ".", extension)
  )
}

log_message(
  "Fitting a one-dimensional principal curve through the raw FPC scores..."
)
continuum_fit <- fit_score_principal_curve(
  score_input,
  spline_df = args$continuum_df,
  max_iterations = args$continuum_max_iterations
)
# A one-dimensional curve has no intrinsic direction. Age is used only to name
# the increasing end; reversing an axis changes neither its fit nor ordering.
continuum_fit <- orient_principal_curve(continuum_fit, score_table$baseline_age)
linear_latent_coordinate <- continuum_fit$linear_coordinate
linear_latent_coordinate <- (
  linear_latent_coordinate - min(linear_latent_coordinate)
) / diff(range(linear_latent_coordinate))
score_table[, `:=`(
  latent_coordinate = continuum_fit$coordinate,
  linear_latent_coordinate = linear_latent_coordinate,
  latent_curve_residual = continuum_fit$residual_distance
)]
continuum_diagnostics <- data.table(
  method = "fixed-complexity principal curve in raw FPC-score space",
  requested_input_fpcs = args$score_pcs,
  input_fpcs = score_dimension,
  input_retained_product_variance_fraction = product_index$pve[[n_score_pcs]],
  spline_df = continuum_fit$spline_df,
  iterations = continuum_fit$iterations,
  principal_curve_variance_explained = continuum_fit$variance_explained,
  linear_axis_variance_explained = continuum_fit$linear_variance_explained,
  nonlinear_improvement = continuum_fit$variance_explained -
    continuum_fit$linear_variance_explained,
  curve_vs_linear_spearman = cor(
    continuum_fit$coordinate, linear_latent_coordinate,
    method = "spearman"
  ),
  coordinate_vs_age_spearman = cor(
    continuum_fit$coordinate, score_table$baseline_age,
    method = "spearman", use = "complete.obs"
  ),
  coordinate_vs_age_pearson = cor(
    continuum_fit$coordinate, score_table$baseline_age,
    method = "pearson", use = "complete.obs"
  ),
  orientation_rule = "direction only chosen for nonnegative age rank correlation"
)

log_message("Comparing PCA, UMAP, and t-SNE embeddings...")
embedding_fit <- compare_score_embeddings(
  score_input,
  neighbors = args$umap_neighbors,
  min_dist = args$umap_min_dist,
  tsne_perplexity = args$tsne_perplexity,
  seed = args$seed + 3L
)
score_table[, `:=`(
  PCA1 = embedding_fit$embeddings$PCA[, 1L],
  PCA2 = embedding_fit$embeddings$PCA[, 2L],
  UMAP1 = embedding_fit$embeddings$UMAP[, 1L],
  UMAP2 = embedding_fit$embeddings$UMAP[, 2L],
  TSNE1 = embedding_fit$embeddings$`t-SNE`[, 1L],
  TSNE2 = embedding_fit$embeddings$`t-SNE`[, 2L]
)]

latent_landmark_percentiles <- seq(
  0.05, 0.95, length.out = args$continuum_landmarks
)
latent_landmark_positions <- as.numeric(stats::quantile(
  continuum_fit$coordinate,
  probs = latent_landmark_percentiles,
  names = FALSE
))
latent_landmark_scores <- interpolate_curve_scores(
  continuum_fit, latent_landmark_positions
)
baseline_time_index <- which.min(abs(t_grid))
reconstruct_latent_scores <- function(ordered_scores) {
  complete_scores <- numeric(ncol(score_matrix_all))
  complete_scores[score_columns[score_input_positions]] <- ordered_scores
  score_array <- array(
    complete_scores,
    dim = c(ranks[["row"]], ranks[["column"]], ranks[["time"]])
  )
  centered <- rTensor::ttl(
    rTensor::as.tensor(score_array),
    list(
      eig$funcs$Phi,
      eig$funcs$Psi,
      eig$funcs$Theta[baseline_time_index, , drop = FALSE]
    ),
    ms = 1:3
  )@data
  dim(centered) <- image_dim
  (centered + predict_mean_matrix(t_grid[[baseline_time_index]])) * brain_mask
}
latent_reconstructed_images <- lapply(
  seq_along(latent_landmark_positions),
  function(index) reconstruct_latent_scores(latent_landmark_scores[index, ])
)

chosen_subject_rows <- integer()
observed_landmark_rows <- vapply(latent_landmark_positions, function(target) {
  candidates <- setdiff(seq_len(nrow(score_table)), chosen_subject_rows)
  ordering <- order(
    abs(score_table$latent_coordinate[candidates] - target),
    score_table$latent_curve_residual[candidates],
    score_table$PTID[candidates]
  )
  selected <- candidates[[ordering[[1L]]]]
  chosen_subject_rows <<- c(chosen_subject_rows, selected)
  selected
}, integer(1L))
latent_landmarks <- data.table(
  target_percentile = latent_landmark_percentiles,
  target_coordinate = latent_landmark_positions,
  PTID = score_table$PTID[observed_landmark_rows],
  subject_coordinate = score_table$latent_coordinate[observed_landmark_rows],
  baseline_age = score_table$baseline_age[observed_landmark_rows],
  curve_residual = score_table$latent_curve_residual[observed_landmark_rows]
)
latent_landmarks[, analysis_row := vapply(PTID, function(subject) {
  candidates <- meta[PTID == subject]
  candidates$analysis_row[[which.min(abs(candidates$time_year))]]
}, integer(1L))]
latent_landmarks[, `:=`(
  image_id = meta$image_id[analysis_row],
  image_time_year = meta$time_year[analysis_row],
  image_file = basename(meta$image_file[analysis_row]),
  source_image_path = meta$image_file[analysis_row]
)]

baseline_subject_image_rows <- vapply(subject_ids, function(subject) {
  candidates <- meta[PTID == subject]
  candidates$analysis_row[[which.min(abs(candidates$time_year))]]
}, integer(1L))
continuum_smoother_neighbors <- min(
  nrow(score_table), max(20L, as.integer(round(0.15 * nrow(score_table))))
)
continuum_kernel_bandwidth <- vapply(
  latent_landmark_positions,
  function(position) {
    sorted_distance <- sort(abs(score_table$latent_coordinate - position))
    max(sorted_distance[[continuum_smoother_neighbors]], 0.03)
  },
  numeric(1L)
)
continuum_effective_subjects <- numeric(length(latent_landmark_positions))
smoothed_observed_images <- lapply(seq_along(latent_landmark_positions), function(index) {
  distance <- score_table$latent_coordinate - latent_landmark_positions[[index]]
  weights <- exp(-0.5 * (distance / continuum_kernel_bandwidth[[index]])^2)
  weights <- weights / sum(weights)
  continuum_effective_subjects[[index]] <<- 1 / sum(weights^2)
  averaged <- as.numeric(crossprod(
    weights,
    image_matrix[baseline_subject_image_rows, , drop = FALSE]
  ))
  matrix(averaged, nrow = image_dim[[1L]]) * brain_mask
})
latent_landmarks[, `:=`(
  smoothing_neighbors = continuum_smoother_neighbors,
  kernel_bandwidth = continuum_kernel_bandwidth,
  kernel_effective_subjects = continuum_effective_subjects
)]

latent_curve_table <- data.table(
  latent_coordinate = continuum_fit$curve_coordinate
)
latent_curve_table <- cbind(
  latent_curve_table,
  as.data.table(continuum_fit$curve_scores)
)
setnames(
  latent_curve_table,
  old = names(latent_curve_table)[-1L],
  new = colnames(score_input)
)

write_image_panel <- function(z, palette, zlim = range(z, finite = TRUE), main = "") {
  if (!all(is.finite(zlim)) || diff(zlim) == 0) zlim <- zlim + c(-1, 1) * 1e-8
  image(
    seq_len(ncol(z)), seq_len(nrow(z)), t(z[nrow(z):1L, , drop = FALSE]),
    col = palette, zlim = zlim, axes = FALSE, useRaster = TRUE,
    xlab = "", ylab = "", main = main
  )
  box()
}

gray_palette <- gray.colors(256L, start = 0, end = 1)
diverging_palette <- colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(255L)

png(
  file.path(args$output_dir, "inferred_brain_mask.png"),
  width = 700L, height = 700L, res = 140L
)
par(mar = c(1, 1, 3, 1))
write_image_panel(
  brain_mask * 1, gray_palette, c(0, 1),
  sprintf("Consensus brain mask: training only (%.1f%%)", 100 * mean(brain_mask))
)
dev.off()

# Publication panel: the two variance-bearing temporal modes form the rows,
# while the four leading marginal spatial product eigenfunctions form the
# columns. No eigenvalue or temporal loading is folded into either function.
spatial_values <- pmax(eig_full$values$lam_phi, 0) %o%
  pmax(eig_full$values$lam_psi, 0)
spatial_order <- order(as.vector(spatial_values), decreasing = TRUE)
spatial_indices <- arrayInd(spatial_order, dim(spatial_values))
spatial_panel_count <- min(4L, args$top_fpcs, nrow(spatial_indices))
spatial_panel_indices <- spatial_indices[seq_len(spatial_panel_count), , drop = FALSE]
spatial_panel_fve <- as.vector(spatial_values)[
  spatial_order[seq_len(spatial_panel_count)]
] / sum(spatial_values)
spatial_panel_functions <- lapply(seq_len(spatial_panel_count), function(k) {
  outer(
    eig_full$funcs$Phi[, spatial_panel_indices[k, 1L]],
    eig_full$funcs$Psi[, spatial_panel_indices[k, 2L]]
  )
})
spatial_panel_limit <- max(abs(unlist(lapply(
  spatial_panel_functions, function(z) z[brain_mask]
))))

orient_temporal_for_display <- function(profile, mode) {
  weighted_mean <- sum(temporal_weights * profile) / sum(temporal_weights)
  if (mode == 1L && weighted_mean < 0) return(-profile)
  if (mode > 1L &&
      profile[[length(profile)]] < profile[[temporal_baseline_index]]) {
    return(-profile)
  }
  profile
}

publication_temporal_modes <- seq_len(min(2L, ncol(temporal_functions)))
publication_temporal_functions <- lapply(
  publication_temporal_modes,
  function(k) orient_temporal_for_display(temporal_functions[, k], k)
)
temporal_panel_ylim <- range(unlist(publication_temporal_functions))
temporal_panel_padding <- 0.04 * diff(temporal_panel_ylim)
temporal_panel_ylim <- temporal_panel_ylim + c(-1, 1) * temporal_panel_padding

write_masked_spatial_panel <- function(z, main) {
  displayed <- z
  displayed[!brain_mask] <- NA_real_
  image(
    seq_len(ncol(displayed)), seq_len(nrow(displayed)),
    t(displayed[nrow(displayed):1L, , drop = FALSE]),
    col = diverging_palette,
    zlim = c(-spatial_panel_limit, spatial_panel_limit),
    axes = FALSE, useRaster = TRUE, xlab = "", ylab = "", main = main
  )
  contour(
    seq_len(ncol(brain_mask)), seq_len(nrow(brain_mask)),
    t(brain_mask[nrow(brain_mask):1L, , drop = FALSE] * 1),
    levels = 0.5, add = TRUE, drawlabels = FALSE,
    col = "black", lwd = 1.1, lty = 1
  )
  box(col = "black", lwd = 0.8)
}

png(
  file.path(args$output_dir, "top_fpc_spatiotemporal_panels.png"),
  width = 2100L, height = 820L, res = 160L, bg = "white"
)
layout(matrix(
  seq_len(length(publication_temporal_modes) * (spatial_panel_count + 1L)),
  nrow = length(publication_temporal_modes), byrow = TRUE
))
par(
  mar = c(3.1, 3.2, 2.7, 0.8), mgp = c(1.8, 0.55, 0),
  tcl = -0.25, cex.axis = 0.72, cex.lab = 0.78, cex.main = 0.84,
  font.main = 1
)
for (row in seq_along(publication_temporal_modes)) {
  time_mode <- publication_temporal_modes[[row]]
  plot(
    t_grid, publication_temporal_functions[[row]],
    type = "l", lwd = 2.2, col = "black",
    xlim = range(t_grid), ylim = temporal_panel_ylim,
    xlab = "Years from baseline", ylab = "Eigenfunction",
    main = sprintf(
      "Temporal FPC %d (%.2f%%)",
      time_mode,
      100 * temporal_diagnostics[
        time_component == time_mode, marginal_variance_fraction
      ]
    )
  )
  abline(h = 0, col = "grey75", lty = 3, lwd = 0.8)
  box(col = "black", lwd = 0.8)
  for (spatial_mode in seq_len(spatial_panel_count)) {
    write_masked_spatial_panel(
      spatial_panel_functions[[spatial_mode]],
      sprintf(
        "Spatial FPC %d (%.2f%%)",
        spatial_mode, 100 * spatial_panel_fve[[spatial_mode]]
      )
    )
  }
}
dev.off()

temporal_plot_modes <- temporal_diagnostics[
  variance_bearing == TRUE, time_component
]
if (!length(temporal_plot_modes)) temporal_plot_modes <- 1L
png(
  file.path(args$output_dir, "temporal_eigencomponent_diagnostics.png"),
  width = 620L * length(temporal_plot_modes), height = 470L, res = 120L
)
par(mfrow = c(1L, length(temporal_plot_modes)), mar = c(4, 4.2, 4, 1))
for (k in temporal_plot_modes) {
  profile <- temporal_functions[, k]
  if (profile[[length(profile)]] < profile[[temporal_baseline_index]]) {
    profile <- -profile
  }
  rms <- sqrt(sum(temporal_weights * profile^2) / sum(temporal_weights))
  diagnostic_row <- temporal_diagnostics[time_component == k]
  plot(
    t_grid, profile / rms, type = "l", lwd = 2.4,
    col = if (diagnostic_row$dynamic_screen) "#D95F0E" else "#2C7FB8",
    xlab = "Years from baseline", ylab = "Eigenfunction (unit RMS)",
    main = sprintf(
      "Temporal mode %d: %.3f%% PVE\nnon-constant fraction = %.3f%%",
      k, 100 * diagnostic_row$marginal_variance_fraction,
      100 * diagnostic_row$nonconstant_fraction
    )
  )
  abline(v = 0, h = 0, lty = 3, col = "grey60")
}
dev.off()

continuum_palette <- grDevices::hcl.colors(256L, palette = "Viridis")
continuum_color <- function(coordinate) {
  index <- 1L + floor(pmin(1, pmax(0, coordinate)) * 255L)
  continuum_palette[index]
}
subject_continuum_colors <- continuum_color(score_table$latent_coordinate)

draw_continuum_embedding <- function(embedding, method_name) {
  plot(
    embedding,
    type = "n",
    xlab = paste(method_name, "1"), ylab = paste(method_name, "2"),
    main = method_name
  )
  points(
    embedding,
    pch = 19,
    col = grDevices::adjustcolor(subject_continuum_colors, alpha.f = 0.78),
    cex = 0.90
  )
}

png(
  continuum_output_file("fpc_score_embedding_comparison_by_continuum", "png"),
  width = 1800L, height = 620L, res = 140L
)
par(mfrow = c(1L, 3L), mar = c(4.2, 4.2, 4.2, 1), oma = c(0, 0, 2.5, 0))
for (method in names(embedding_fit$embeddings)) {
  draw_continuum_embedding(embedding_fit$embeddings[[method]], method)
  legend_values <- seq(0, 1, length.out = 5L)
  legend(
    "topright",
    legend = sprintf("%.2f", legend_values),
    col = continuum_color(legend_values),
    pch = 19, pt.cex = 0.9, bty = "n", cex = 0.72,
    title = "coordinate"
  )
}
mtext(
  sprintf(
    "Continuous morphology coordinate: principal curve through %d raw FPC scores",
    ncol(score_input)
  ),
  side = 3L, outer = TRUE, line = 0.6, cex = 0.95
)
dev.off()

png(
  continuum_output_file("fpc_score_umap_by_continuum", "png"),
  width = 1050L, height = 850L, res = 140L
)
par(mar = c(4.5, 4.5, 5, 1))
draw_continuum_embedding(embedding_fit$embeddings$UMAP, "UMAP")
mtext(
  sprintf(
    "Curve explains %.1f%% of score variation; coordinate-age Spearman %.3f",
    100 * continuum_fit$variance_explained,
    continuum_diagnostics$coordinate_vs_age_spearman
  ),
  side = 3L, line = 0.45, cex = 0.82
)
legend_values <- seq(0, 1, length.out = 6L)
legend(
  "topright", legend = sprintf("%.1f", legend_values),
  col = continuum_color(legend_values), pch = 19, bty = "n",
  title = "morphology coordinate"
)
dev.off()

png(
  continuum_output_file("reconstructed_mri_continuum", "png"),
  width = 390L * args$continuum_landmarks,
  height = 410L,
  res = 120L, bg = "white"
)
par(
  mfrow = c(1L, args$continuum_landmarks),
  mar = c(1, 1, 2.8, 1), cex.main = 1.05
)
for (index in seq_along(latent_landmark_positions)) {
  write_image_panel(
    latent_reconstructed_images[[index]], gray_palette, c(0, 1),
    sprintf("Coordinate %.2f", latent_landmark_positions[[index]])
  )
}
dev.off()

png(
  continuum_output_file("smoothed_observed_mri_continuum", "png"),
  width = 390L * args$continuum_landmarks,
  height = 410L,
  res = 120L, bg = "white"
)
par(
  mfrow = c(1L, args$continuum_landmarks),
  mar = c(1, 1, 2.8, 1), cex.main = 1.05
)
for (index in seq_along(latent_landmark_positions)) {
  write_image_panel(
    smoothed_observed_images[[index]], gray_palette, c(0, 1),
    sprintf("Coordinate %.2f", latent_landmark_positions[[index]])
  )
}
dev.off()

png(
  continuum_output_file("observed_mri_continuum", "png"),
  width = 390L * args$continuum_landmarks,
  height = 410L,
  res = 120L, bg = "white"
)
par(
  mfrow = c(1L, args$continuum_landmarks),
  mar = c(1, 1, 2.8, 1), cex.main = 1.05
)
for (index in seq_len(nrow(latent_landmarks))) {
  landmark <- latent_landmarks[index]
  observed <- matrix(
    image_matrix[landmark$analysis_row, ], nrow = image_dim[[1L]]
  )
  write_image_panel(
    observed, gray_palette, c(0, 1),
    sprintf("Coordinate %.2f", landmark$subject_coordinate)
  )
}
dev.off()

png(
  continuum_output_file("latent_coordinate_vs_age", "png"),
  width = 1050L, height = 780L, res = 140L
)
par(mar = c(5, 4.8, 4.8, 1))
plot(
  score_table$latent_coordinate,
  score_table$baseline_age,
  pch = 19, cex = 0.72,
  col = grDevices::adjustcolor(subject_continuum_colors, alpha.f = 0.60),
  xlab = "Continuous morphology coordinate",
  ylab = "Baseline age (years)",
  main = "Age varies along, but does not define, the MRI morphology continuum"
)
age_smoother <- smooth.spline(
  score_table$latent_coordinate,
  score_table$baseline_age,
  df = 4
)
coordinate_grid <- seq(0, 1, length.out = 300L)
lines(
  coordinate_grid,
  predict(age_smoother, coordinate_grid)$y,
  col = "#333333", lwd = 3
)
mtext(
  sprintf(
    "Spearman rho = %.3f; coordinate learned from MRI scores without age",
    continuum_diagnostics$coordinate_vs_age_spearman
  ),
  side = 3L, line = 0.5, cex = 0.85
)
dev.off()

analysis_manifest <- copy(meta)
analysis_manifest[, image_file := basename(image_file)]
score_metadata_columns <- c(
  "PTID", "baseline_age", "latent_coordinate", "linear_latent_coordinate",
  "latent_curve_residual",
  "PCA1", "PCA2", "UMAP1", "UMAP2", "TSNE1", "TSNE2"
)
setcolorder(score_table, c(
  score_metadata_columns, setdiff(names(score_table), score_metadata_columns)
))
fwrite(analysis_manifest, file.path(args$output_dir, "analysis_manifest.csv"))
fwrite(prepared$exclusions, file.path(args$output_dir, "exclusion_summary.csv"))
fwrite(reconstruction_metrics, file.path(args$output_dir, "reconstruction_metrics.csv"))
fwrite(reconstruction_summary, file.path(args$output_dir, "reconstruction_summary.csv"))
fwrite(product_index, file.path(args$output_dir, "product_fpc_eigenvalues.csv"))
fwrite(
  temporal_diagnostics,
  file.path(args$output_dir, "temporal_eigencomponent_diagnostics.csv")
)
fwrite(score_table, continuum_output_file("fpc_scores_embeddings", "csv"))
fwrite(embedding_fit$metrics,
  continuum_output_file("embedding_quality_metrics", "csv"))
fwrite(continuum_diagnostics,
  continuum_output_file("latent_coordinate_diagnostics", "csv"))
fwrite(latent_landmarks,
  continuum_output_file("latent_coordinate_landmarks", "csv"))
fwrite(latent_curve_table,
  continuum_output_file("latent_principal_curve_scores", "csv"))

compact_model <- list(
  covariance = covariance_fit$C,
  eigen = eig,
  mean_model = mean_model,
  brain_mask = brain_mask,
  foreground_frequency = mask_fit$foreground_frequency,
  image_dim = image_dim,
  t_grid_years = t_grid,
  ranks = ranks,
  temporal_eigencomponent_diagnostics = temporal_diagnostics,
  subject_ids = subject_ids,
  embeddings = embedding_fit$embeddings,
  embedding_metrics = embedding_fit$metrics,
  morphology_continuum = list(
    input_fpcs = ncol(score_input),
    score_scaling = "none",
    coordinate = setNames(continuum_fit$coordinate, subject_ids),
    residual_distance = setNames(continuum_fit$residual_distance, subject_ids),
    curve_coordinate = continuum_fit$curve_coordinate,
    curve_scores = continuum_fit$curve_scores,
    diagnostics = continuum_diagnostics,
    landmarks = latent_landmarks
  ),
  arguments = args
)
saveRDS(compact_model, continuum_output_file("adni_matfpca_model", "rds"))
suppressWarnings(capture.output(
  sessionInfo(), file = continuum_output_file("session_info", "txt")
))

run_summary <- c(
  sprintf("Manifest: %s", args$manifest),
  sprintf("Eligible subjects: %d", length(subject_ids)),
  sprintf("Eligible images: %d", nrow(meta)),
  sprintf("Training images: %d", sum(meta$split == "train")),
  sprintf("Held-out images/subjects: %d", sum(meta$split == "test")),
  sprintf("Image dimension after pooling: %s", paste(image_dim, collapse = " x ")),
  sprintf("Training-inferred brain-mask coverage: %.2f%%", 100 * mean(brain_mask)),
  sprintf("Time grid in years: %.6f to %.6f by %.6f", min(t_grid), max(t_grid), args$time_step_months / 12),
  sprintf("MatFPCA mode ranks: %s", paste(ranks, collapse = " x ")),
  sprintf(
    "Leading temporal mode: %.4f%% marginal PVE; %.4f%% non-constant energy",
    100 * temporal_diagnostics$marginal_variance_fraction[[1L]],
    100 * temporal_diagnostics$nonconstant_fraction[[1L]]
  ),
  if (any(temporal_diagnostics$dynamic_screen)) {
    paste(
      "Material non-constant temporal mode(s):",
      paste(vapply(
        which(temporal_diagnostics$dynamic_screen),
        function(row) sprintf(
          "%d (%.4f%% marginal PVE; %.2f%% non-constant energy; leading full product rank %d)",
          temporal_diagnostics$time_component[[row]],
          100 * temporal_diagnostics$marginal_variance_fraction[[row]],
          100 * temporal_diagnostics$nonconstant_fraction[[row]],
          temporal_diagnostics$leading_full_product_rank[[row]]
        ),
        character(1L)
      ), collapse = "; ")
    )
  } else {
    "Material non-constant temporal modes: none under the prespecified diagnostic screen."
  },
  if (any(temporal_diagnostics$dynamic_screen &
          temporal_diagnostics$time_component > ranks[["time"]])) {
    sprintf(
      "Temporal truncation implication: screened dynamic mode(s) lie beyond retained K=%d; the leading retained product FPCs therefore use the near-constant temporal mode.",
      ranks[["time"]]
    )
  } else {
    "Temporal truncation implication: all screened dynamic modes are retained."
  },
  sprintf("Estimated noise variance: %.8g", covariance_fit$C$sig2),
  sprintf("MatFPCA elapsed seconds: %.3f", analysis_fit$elapsed_seconds),
  sprintf(
    "Continuum input: %d leading raw FPC scores (%.2f%% of retained product variance)",
    ncol(score_input), 100 * product_index$pve[[n_score_pcs]]
  ),
  sprintf(
    "Morphology continuum variance explained (principal curve / linear axis): %.2f%% / %.2f%%",
    100 * continuum_fit$variance_explained,
    100 * continuum_fit$linear_variance_explained
  ),
  sprintf(
    "Principal-curve vs linear-axis rank correlation: %.4f",
    continuum_diagnostics$curve_vs_linear_spearman
  ),
  sprintf(
    "Morphology-coordinate vs age Spearman correlation: %.4f",
    continuum_diagnostics$coordinate_vs_age_spearman
  ),
  "Within this script, clinical outcome labels were not used for subject sampling, holdout selection, continuum fitting, embedding, plotting, or diagnostics.",
  "Primary interpretation: dominant MRI-score variation follows an almost linear continuum associated with, but not defined by, age.",
  "PCA, UMAP, and t-SNE are visualization checks of the continuous score geometry.",
  "",
  "Embedding quality:",
  paste(capture.output(print(embedding_fit$metrics)), collapse = "\n"),
  "",
  "Reconstruction summary (metrics evaluated inside the inferred brain mask):",
  paste(capture.output(print(reconstruction_summary)), collapse = "\n")
)
writeLines(run_summary, continuum_output_file("run_summary", "txt"))

log_message("Analysis complete. Results: ", args$output_dir)
cat(paste(run_summary, collapse = "\n"), "\n")
