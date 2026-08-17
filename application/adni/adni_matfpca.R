#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(png)
  library(rTensor)
})
if (!nzchar(Sys.getenv("TZ"))) Sys.setenv(TZ = "America/Vancouver")

required_packages <- c(
  "mgcv", "pracma", "foreach", "splines", "uwot", "FNN", "class", "cluster"
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
  umap_pcs = 30L,
  umap_neighbors = 30L,
  umap_min_dist = 0.20,
  probe_neighbors = 15L,
  probe_folds = 5L,
  top_fpcs = 5L,
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
    "  --holdout-fraction=0.20     Stratified fraction of subjects tested\n",
    "  --time-step-months=VALUE    Visit-time resolution (default: one day)\n",
    "  --umap-pcs=30               Leading standardized FPC scores used by UMAP\n",
    "  --umap-neighbors=30         UMAP neighborhood size\n",
    "  --max-subjects=Inf          Stratified subject cap for smoke tests\n",
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
if (args$pool < 1L || args$mean_df < 3L || args$umap_pcs < 2L) {
  stop("pool >= 1, mean-df >= 3, and umap-pcs >= 2 are required.")
}
if (args$mask_intensity < 0 || args$mask_min_fraction <= 0 ||
    args$mask_min_fraction >= 1 || args$mask_dilation < 0L) {
  stop("Invalid brain-mask settings.")
}
if (args$umap_neighbors < 2L || args$umap_min_dist < 0 ||
    args$probe_neighbors < 1L || args$probe_folds < 2L) {
  stop("Invalid UMAP or representation-probe settings.")
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

# These files belong to the superseded mixture-clustering analysis. Remove only
# the exact generated artifacts so stale results cannot be mistaken for current
# output.
legacy_cluster_outputs <- file.path(args$output_dir, c(
  "fpc_score_clustering.png", "fpc_scores_and_clusters.csv",
  "clustering_metrics.csv", "clustering_contingency.csv"
))
invisible(file.remove(legacy_cluster_outputs[file.exists(legacy_cluster_outputs)]))

log_message <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
  flush.console()
}

safe_label <- function(x) {
  result <- rep(NA_character_, length(x))
  result[x == "prevalent_ad_at_baseline"] <- "AD at baseline"
  result[x == "incident_ad_conversion"] <- "MCI-to-AD converter"
  result[x == "right_censored_non_ad"] <- "MCI non-converter"
  result[is.na(result)] <- x[is.na(result)]
  result
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

stratified_subject_sample <- function(subject_table, size, seed) {
  if (!is.finite(size) || size >= nrow(subject_table)) return(subject_table$PTID)
  size <- max(2L, as.integer(size))
  set.seed(seed)
  allocation <- subject_table[, .N, by = patient_group]
  allocation[, raw := size * N / sum(N)]
  allocation[, take := pmin(N, pmax(1L, floor(raw)))]
  while (sum(allocation$take) < size) {
    candidates <- which(allocation$take < allocation$N)
    if (!length(candidates)) break
    target <- candidates[which.max(allocation$raw[candidates] - allocation$take[candidates])]
    allocation$take[target] <- allocation$take[target] + 1L
  }
  while (sum(allocation$take) > size) {
    candidates <- which(allocation$take > 1L)
    if (!length(candidates)) break
    target <- candidates[which.min(allocation$raw[candidates] - allocation$take[candidates])]
    allocation$take[target] <- allocation$take[target] - 1L
  }
  unlist(lapply(seq_len(nrow(allocation)), function(i) {
    pool <- subject_table[patient_group == allocation$patient_group[[i]], PTID]
    sample(pool, allocation$take[[i]])
  }), use.names = FALSE)
}

prepare_manifest <- function() {
  manifest <- data.table::fread(args$manifest, na.strings = c("", "NA"))
  required <- c(
    "PTID", "image_id", "image_time_month", "mri_VISCODE2",
    "baseline_dx", "dx_at_or_near_mri", "outcome_type",
    "image_retained_for_analysis"
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

  retained[, patient_group := safe_label(outcome_type)]
  retained[, patient_group := factor(
    patient_group,
    levels = c("AD at baseline", "MCI-to-AD converter", "MCI non-converter")
  )]
  if (anyNA(retained$patient_group)) {
    stop("Unrecognized outcome_type values remain after filtering.")
  }

  subjects <- unique(retained[, .(PTID, patient_group)])
  if (anyDuplicated(subjects$PTID)) stop("A subject has inconsistent patient groups.")
  selected_subjects <- stratified_subject_sample(
    subjects, args$max_subjects, args$seed + 10L
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

# Select subjects for one-visit holdout, stratified by clinical outcome group.
subject_table <- unique(meta[, .(PTID, patient_group)])
set.seed(args$seed)
holdout_subjects <- unlist(lapply(split(subject_table, subject_table$patient_group), function(x) {
  number <- max(1L, as.integer(round(args$holdout_fraction * nrow(x))))
  sample(x$PTID, min(number, nrow(x)))
}), use.names = FALSE)
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
      patient_group = as.character(meta$patient_group[[row]]),
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
), by = .(split, patient_group)]

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
  patient_group = as.character(subject_table$patient_group[match(subject_ids, subject_table$PTID)])
)
score_table <- cbind(score_table, as.data.table(score_matrix_ordered))

n_umap_pcs <- min(args$umap_pcs, ncol(score_matrix_ordered), nrow(score_matrix_ordered) - 2L)
umap_input <- scale(score_matrix_ordered[, seq_len(n_umap_pcs), drop = FALSE])
valid_dimensions <- apply(
  umap_input, 2L,
  function(x) all(is.finite(x)) && stats::sd(x) > 1e-10
)
umap_input <- umap_input[, valid_dimensions, drop = FALSE]
if (ncol(umap_input) < 2L) stop("Too few non-degenerate FPC score dimensions.")

umap_neighbors <- min(args$umap_neighbors, nrow(umap_input) - 1L)
log_message(
  "Embedding ", ncol(umap_input), " standardized FPC scores with UMAP (",
  umap_neighbors, " neighbors, min_dist = ", args$umap_min_dist, ")..."
)
set.seed(args$seed + 3L)
umap_embedding <- uwot::umap(
  umap_input,
  n_neighbors = umap_neighbors,
  n_components = 2L,
  metric = "euclidean",
  min_dist = args$umap_min_dist,
  init = "spectral",
  n_threads = 1L,
  n_sgd_threads = 1L,
  fast_sgd = FALSE,
  ret_model = FALSE,
  verbose = FALSE,
  seed = args$seed + 3L
)
colnames(umap_embedding) <- c("UMAP1", "UMAP2")
score_table[, `:=`(UMAP1 = umap_embedding[, 1L], UMAP2 = umap_embedding[, 2L])]
true_groups <- score_table$patient_group

label_silhouette <- function(representation, labels) {
  group_id <- as.integer(factor(labels))
  mean(cluster::silhouette(group_id, stats::dist(representation))[, "sil_width"])
}

neighbor_label_diagnostics <- function(representation, labels, k, space) {
  k <- min(k, nrow(representation) - 1L)
  neighbors <- FNN::get.knn(representation, k = k)$nn.index
  neighbor_labels <- matrix(labels[neighbors], nrow = nrow(representation))
  same_group_fraction <- rowMeans(neighbor_labels == labels)
  group_sizes <- table(labels)
  result <- data.table(
    patient_group = labels,
    same_group_neighbor_fraction = same_group_fraction
  )[, .(
    subjects = .N,
    same_group_neighbor_fraction = mean(same_group_neighbor_fraction)
  ), by = patient_group]
  result[, chance_fraction :=
    (as.numeric(group_sizes[patient_group]) - 1) / (length(labels) - 1)]
  result[, enrichment_over_chance := same_group_neighbor_fraction / chance_fraction]
  result[, `:=`(space = space, neighbors = k)]
  setcolorder(result, c(
    "space", "neighbors", "patient_group", "subjects",
    "same_group_neighbor_fraction", "chance_fraction", "enrichment_over_chance"
  ))
  result
}

cross_validated_knn_probe <- function(
  representation, labels, folds, neighbors, seed
) {
  label_levels <- unique(labels)
  fold_id <- integer(length(labels))
  set.seed(seed)
  for (label in label_levels) {
    indices <- which(labels == label)
    fold_id[indices] <- sample(rep(seq_len(folds), length.out = length(indices)))
  }
  predictions <- rep(NA_character_, length(labels))
  for (fold in seq_len(folds)) {
    test <- which(fold_id == fold)
    train <- which(fold_id != fold)
    k <- min(neighbors, length(train))
    predictions[test] <- as.character(class::knn(
      train = representation[train, , drop = FALSE],
      test = representation[test, , drop = FALSE],
      cl = factor(labels[train], levels = label_levels),
      k = k
    ))
  }
  confusion <- table(
    truth = factor(labels, levels = label_levels),
    prediction = factor(predictions, levels = label_levels)
  )
  recall <- diag(confusion) / rowSums(confusion)
  list(
    accuracy = mean(predictions == labels),
    balanced_accuracy = mean(recall),
    recall = data.table(
      patient_group = label_levels,
      subjects = as.integer(rowSums(confusion)),
      recall = as.numeric(recall)
    ),
    confusion = as.data.table(confusion)
  )
}

probe_folds <- min(
  args$probe_folds,
  max(2L, min(as.integer(table(true_groups))))
)
probe <- cross_validated_knn_probe(
  umap_input,
  true_groups,
  folds = probe_folds,
  neighbors = args$probe_neighbors,
  seed = args$seed + 4L
)
mci_indices <- which(true_groups != "AD at baseline")
mci_probe <- cross_validated_knn_probe(
  umap_input[mci_indices, , drop = FALSE],
  true_groups[mci_indices],
  folds = min(args$probe_folds, min(table(true_groups[mci_indices]))),
  neighbors = args$probe_neighbors,
  seed = args$seed + 5L
)
probe_recall <- rbindlist(list(
  cbind(task = "three_group", probe$recall),
  cbind(task = "MCI_conversion_binary", mci_probe$recall)
))
probe_confusion <- rbindlist(list(
  cbind(task = "three_group", probe$confusion),
  cbind(task = "MCI_conversion_binary", mci_probe$confusion)
))
representation_metrics <- data.table(
  method = "Unsupervised UMAP of standardized product-FPC scores",
  fpc_dimensions = ncol(umap_input),
  umap_neighbors = umap_neighbors,
  umap_min_dist = args$umap_min_dist,
  known_group_silhouette_fpc = label_silhouette(umap_input, true_groups),
  known_group_silhouette_umap = label_silhouette(umap_embedding, true_groups),
  knn_probe_neighbors = args$probe_neighbors,
  knn_probe_folds = probe_folds,
  knn_probe_accuracy = probe$accuracy,
  knn_probe_balanced_accuracy = probe$balanced_accuracy,
  mci_binary_knn_probe_accuracy = mci_probe$accuracy,
  mci_binary_knn_probe_balanced_accuracy = mci_probe$balanced_accuracy,
  majority_class_accuracy = max(table(true_groups)) / length(true_groups),
  balanced_chance_accuracy = 1 / length(unique(true_groups))
)
neighbor_diagnostics <- rbindlist(list(
  neighbor_label_diagnostics(
    umap_input, true_groups, args$probe_neighbors, "standardized_FPC"
  ),
  neighbor_label_diagnostics(
    umap_embedding, true_groups, args$probe_neighbors, "UMAP_2D"
  )
))

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
diverging_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(255L)

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

top_count <- min(args$top_fpcs, nrow(product_index))
png(
  file.path(args$output_dir, "top_fpc_spatial_panels.png"),
  width = 360L * top_count, height = 390L, res = 120L
)
par(mfrow = c(1L, top_count), mar = c(1.5, 1.5, 3.5, 1))
for (component in seq_len(top_count)) {
  index <- product_index[component]
  theta <- eig$funcs$Theta[, index$time_component]
  time_position <- which.max(abs(theta))
  effect <- sqrt(index$eigenvalue) *
    outer(
      eig$funcs$Phi[, index$row_component],
      eig$funcs$Psi[, index$column_component]
    ) * theta[[time_position]]
  effect[!brain_mask] <- 0
  limit <- max(abs(effect))
  write_image_panel(
    effect, diverging_palette, c(-limit, limit),
    sprintf("FPC %d\npeak |effect| at t = %.2f y", component, t_grid[[time_position]])
  )
}
dev.off()

png(
  file.path(args$output_dir, "top_fpc_temporal_panels.png"),
  width = 360L * top_count, height = 330L, res = 120L
)
par(mfrow = c(1L, top_count), mar = c(3.5, 3.5, 3, 1))
for (component in seq_len(top_count)) {
  index <- product_index[component]
  profile <- sqrt(index$eigenvalue) * eig$funcs$Theta[, index$time_component]
  plot(
    t_grid, profile, type = "l", lwd = 2, col = "#2C7FB8",
    xlab = "Years from baseline", ylab = "1-SD loading",
    main = sprintf("FPC %d: (%d, %d, %d)", component,
      index$row_component, index$column_component, index$time_component)
  )
  abline(v = 0, h = 0, lty = 3, col = "grey60")
}
dev.off()

group_colors <- c(
  "AD at baseline" = "#D73027",
  "MCI-to-AD converter" = "#FC8D59",
  "MCI non-converter" = "#4575B4"
)
group_shapes <- c(
  "AD at baseline" = 17,
  "MCI-to-AD converter" = 19,
  "MCI non-converter" = 19
)
png(
  file.path(args$output_dir, "fpc_score_umap_by_patient_group.png"),
  width = 1050L, height = 850L, res = 140L
)
par(mar = c(4.5, 4.5, 4.5, 10), xpd = TRUE)
plot(
  umap_embedding,
  type = "n",
  xlab = "UMAP 1", ylab = "UMAP 2",
  main = "UMAP of longitudinal MatFPCA scores"
)
mtext(
  sprintf(
    "%d standardized FPC scores; patient labels used only for color",
    ncol(umap_input)
  ),
  side = 3, line = 0.5, cex = 0.82
)
draw_order <- order(factor(
  true_groups,
  levels = c("MCI non-converter", "MCI-to-AD converter", "AD at baseline")
))
points(
  umap_embedding[draw_order, , drop = FALSE],
  pch = group_shapes[true_groups[draw_order]],
  col = grDevices::adjustcolor(group_colors[true_groups[draw_order]], alpha.f = 0.78),
  cex = ifelse(true_groups[draw_order] == "AD at baseline", 1.15, 0.9)
)
legend(
  "right", inset = c(-0.38, 0),
  legend = names(group_colors),
  col = group_colors,
  pch = group_shapes[names(group_colors)],
  pt.cex = 1.1,
  bty = "n"
)
dev.off()

# Pick the subject nearest the median held-out error in every patient group.
example_subjects <- reconstruction_metrics[split == "test", {
  group_median <- median(rmse)
  .SD[which.min(abs(rmse - group_median)), .(PTID)]
}, by = patient_group]

png(
  file.path(args$output_dir, "reconstruction_examples_by_group.png"),
  width = 1600L, height = 350L * nrow(example_subjects), res = 120L
)
par(mfrow = c(nrow(example_subjects), 5L), mar = c(1, 1, 3.6, 1))
for (example_index in seq_len(nrow(example_subjects))) {
  subject <- example_subjects$PTID[[example_index]]
  test_row <- meta[PTID == subject & split == "test", analysis_row]
  training_candidates <- meta[PTID == subject & split == "train"]
  train_row <- training_candidates$analysis_row[[which.min(abs(training_candidates$time_year))]]
  observed_train <- matrix(image_matrix[train_row, ], nrow = image_dim[[1L]])
  reconstructed_train <- reconstruction_store[[subject]][[as.character(train_row)]]
  observed_test <- matrix(image_matrix[test_row, ], nrow = image_dim[[1L]])
  reconstructed_test <- reconstruction_store[[subject]][[as.character(test_row)]]
  difference <- reconstructed_test - observed_test
  difference_limit <- max(abs(difference))

  write_image_panel(
    observed_train, gray_palette, c(0, 1),
    sprintf("%s\n%s\nTrain original, t = %.2f y",
      example_subjects$patient_group[[example_index]], subject, meta$time_year[[train_row]])
  )
  write_image_panel(
    reconstructed_train, gray_palette, c(0, 1),
    sprintf("Train reconstruction\nRMSE %.3f",
      reconstruction_metrics[analysis_row == train_row, rmse])
  )
  write_image_panel(
    observed_test, gray_palette, c(0, 1),
    sprintf("Held-out original\nt = %.2f y", meta$time_year[[test_row]])
  )
  write_image_panel(
    reconstructed_test, gray_palette, c(0, 1),
    sprintf("Held-out prediction\nRMSE %.3f",
      reconstruction_metrics[analysis_row == test_row, rmse])
  )
  write_image_panel(
    difference, diverging_palette, c(-difference_limit, difference_limit),
    "Prediction - original"
  )
}
dev.off()

analysis_manifest <- copy(meta)
analysis_manifest[, image_file := basename(image_file)]
fwrite(analysis_manifest, file.path(args$output_dir, "analysis_manifest.csv"))
fwrite(prepared$exclusions, file.path(args$output_dir, "exclusion_summary.csv"))
fwrite(reconstruction_metrics, file.path(args$output_dir, "reconstruction_metrics.csv"))
fwrite(reconstruction_summary, file.path(args$output_dir, "reconstruction_summary.csv"))
fwrite(product_index, file.path(args$output_dir, "product_fpc_eigenvalues.csv"))
fwrite(score_table, file.path(args$output_dir, "fpc_scores_and_umap.csv"))
fwrite(representation_metrics,
  file.path(args$output_dir, "umap_representation_diagnostics.csv"))
fwrite(neighbor_diagnostics,
  file.path(args$output_dir, "umap_neighbor_diagnostics.csv"))
fwrite(probe_recall, file.path(args$output_dir, "knn_probe_group_recall.csv"))
fwrite(probe_confusion, file.path(args$output_dir, "knn_probe_confusion.csv"))

compact_model <- list(
  covariance = covariance_fit$C,
  eigen = eig,
  mean_model = mean_model,
  brain_mask = brain_mask,
  foreground_frequency = mask_fit$foreground_frequency,
  image_dim = image_dim,
  t_grid_years = t_grid,
  ranks = ranks,
  subject_ids = subject_ids,
  umap_embedding = umap_embedding,
  representation_metrics = representation_metrics,
  arguments = args
)
saveRDS(compact_model, file.path(args$output_dir, "adni_matfpca_model.rds"))
suppressWarnings(capture.output(
  sessionInfo(), file = file.path(args$output_dir, "session_info.txt")
))

run_summary <- c(
  sprintf("Manifest: %s", args$manifest),
  sprintf("Eligible subjects: %d", length(subject_ids)),
  sprintf("Eligible images: %d", nrow(meta)),
  sprintf("Training images: %d", sum(meta$split == "train")),
  sprintf("Held-out images/subjects: %d", sum(meta$split == "test")),
  sprintf(
    "Patient groups (subjects): %s",
    paste(
      names(table(subject_table$patient_group)), table(subject_table$patient_group),
      sep = "=", collapse = "; "
    )
  ),
  sprintf("Image dimension after pooling: %s", paste(image_dim, collapse = " x ")),
  sprintf("Training-inferred brain-mask coverage: %.2f%%", 100 * mean(brain_mask)),
  sprintf("Time grid in years: %.6f to %.6f by %.6f", min(t_grid), max(t_grid), args$time_step_months / 12),
  sprintf("MatFPCA mode ranks: %s", paste(ranks, collapse = " x ")),
  sprintf("Estimated noise variance: %.8g", covariance_fit$C$sig2),
  sprintf("MatFPCA elapsed seconds: %.3f", analysis_fit$elapsed_seconds),
  sprintf(
    "UMAP input: %d standardized FPC scores; neighbors=%d; min_dist=%.3f",
    ncol(umap_input), umap_neighbors, args$umap_min_dist
  ),
  sprintf(
    "Known-group silhouette (FPC / UMAP): %.4f / %.4f",
    representation_metrics$known_group_silhouette_fpc,
    representation_metrics$known_group_silhouette_umap
  ),
  sprintf(
    "%d-fold CV %d-NN probe accuracy / balanced accuracy: %.4f / %.4f",
    representation_metrics$knn_probe_folds,
    representation_metrics$knn_probe_neighbors,
    representation_metrics$knn_probe_accuracy,
    representation_metrics$knn_probe_balanced_accuracy
  ),
  sprintf(
    "MCI converter-vs-non-converter probe accuracy / balanced accuracy: %.4f / %.4f",
    representation_metrics$mci_binary_knn_probe_accuracy,
    representation_metrics$mci_binary_knn_probe_balanced_accuracy
  ),
  sprintf(
    "Majority / balanced-chance baselines: %.4f / %.4f",
    representation_metrics$majority_class_accuracy,
    representation_metrics$balanced_chance_accuracy
  ),
  "Interpretation: UMAP is exploratory; apparent islands are not evidence of clinical separation without support from the FPC-space diagnostics.",
  "",
  "Reconstruction summary (metrics evaluated inside the inferred brain mask):",
  paste(capture.output(print(reconstruction_summary)), collapse = "\n")
)
writeLines(run_summary, file.path(args$output_dir, "run_summary.txt"))

log_message("Analysis complete. Results: ", args$output_dir)
cat(paste(run_summary, collapse = "\n"), "\n")
