#!/usr/bin/env Rscript

# Interpret already-fitted BrainIAC MatrixFPCA scores and eigenmodes.
#
# This analysis deliberately does not reconstruct MR images or score latent
# reconstruction error.  It asks whether product-FPC scores form continuous
# morphology coordinates and whether those coordinates correspond to regional
# deformation, temporal change, or demographic/clinical features.

suppressPackageStartupMessages({
  library(data.table)
  library(RNifti)
})

required_packages <- c("jsonlite", "png", "uwot", "Rtsne")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]))
} else {
  normalizePath(file.path("application", "adni-latent", "interpret_brainiac_matfpca.R"))
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
setwd(project_root)
source(file.path("application", "adni-latent", "adni_latent_utils.R"))
source(file.path("application", "adni", "adni_continuum_utils.R"))

defaults <- list(
  output_dir = file.path(
    "application", "adni-latent", "results", "brainiac_matfpca_interpretation_full"
  ),
  atlas = file.path(
    "application", "adni-latent", "resources", "templateflow",
    "tpl-MNI152NLin2009cAsym_res-01_seg-aseg_dseg.nii.gz"
  ),
  nonlinear_dir = file.path(
    "/home", "muye", "Projects", "SurvADNI", "ADNI_raw_data",
    "downloaded_images", "preprocessed", "nonlinear"
  ),
  top_components = 20L,
  min_reliability = 0.20,
  minimum_subjects = 250L,
  bootstrap_replicates = 2000L,
  seed = 20260819L,
  cache_targets = TRUE
)

layer_specs <- data.table(
  layer = c("block_3", "block_6", "final"),
  result_dir = file.path(
    "application", "adni-latent", "results",
    c(
      "brainiac_patch_full_layer03",
      "brainiac_patch_full_layer06",
      "brainiac_patch_full"
    )
  )
)

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
  result <- default_values
  for (argument in commandArgs(trailingOnly = TRUE)) {
    if (argument %in% c("--help", "-h")) {
      cat(
        "Usage: Rscript application/adni-latent/interpret_brainiac_matfpca.R ",
        "[--option=value ...]\n",
        "  --top-components=20       Leading product FPCs tested per layer\n",
        "  --min-reliability=0.20     Minimum model-based score reliability\n",
        "  --minimum-subjects=250     Refuse undersized analysis inputs\n",
        "  --bootstrap-replicates=2000\n",
        "  --cache-targets=true       Reuse extracted regional log-Jacobians\n",
        sep = ""
      )
      quit(save = "no", status = 0)
    }
    if (!grepl("^--[^=]+=", argument)) stop("Options require --name=value: ", argument)
    pair <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", pair[[1L]], fixed = TRUE)
    if (!key %in% names(result)) stop("Unknown option: --", pair[[1L]])
    result[[key]] <- coerce_option(paste(pair[-1L], collapse = "="), default_values[[key]])
  }
  result
}

args <- parse_options(defaults)
if (args$top_components < 1L || args$minimum_subjects < 1L ||
    args$bootstrap_replicates < 100L ||
    args$min_reliability < 0 || args$min_reliability > 1) {
  stop("Invalid component, reliability, or bootstrap setting.")
}
args$output_dir <- normalizePath(args$output_dir, mustWork = FALSE)
args$atlas <- normalizePath(args$atlas)
args$nonlinear_dir <- normalizePath(args$nonlinear_dir)
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

log_message <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
  flush.console()
}

roi_labels <- list(
  whole_brain = 1:299,
  cortex = c(3L, 42L),
  white_matter = c(2L, 41L),
  hippocampus = c(17L, 53L),
  lateral_ventricles = c(4L, 5L, 43L, 44L),
  amygdala = c(18L, 54L),
  thalamus = c(10L, 49L)
)
local_regions <- c("hippocampus", "lateral_ventricles", "amygdala", "thalamus")
foreground_roi_start <- c(16L, 16L, 0L)
foreground_roi_end <- c(177L, 217L, 177L)
model_image_size <- c(96L, 96L, 96L)
patch_size <- c(16L, 16L, 16L)
patch_grid <- model_image_size %/% patch_size

nearest_resize_indices <- function(input_size, output_size) {
  floor(seq.int(0L, output_size - 1L) * input_size / output_size) + 1L
}

crop_and_resize_atlas <- function(atlas) {
  crop_indices <- Map(
    function(start, end) seq.int(start + 1L, end),
    foreground_roi_start, foreground_roi_end
  )
  cropped <- atlas[crop_indices[[1L]], crop_indices[[2L]], crop_indices[[3L]]]
  resize_indices <- Map(
    nearest_resize_indices, dim(cropped), model_image_size
  )
  cropped[resize_indices[[1L]], resize_indices[[2L]], resize_indices[[3L]]]
}

atlas_patch_counts <- function(mask) {
  counts <- numeric(prod(patch_grid))
  for (i in 0:(patch_grid[[1L]] - 1L)) {
    for (j in 0:(patch_grid[[2L]] - 1L)) {
      for (k in 0:(patch_grid[[3L]] - 1L)) {
        token <- i * patch_grid[[2L]] * patch_grid[[3L]] +
          j * patch_grid[[3L]] + k + 1L
        counts[[token]] <- sum(mask[
          i * patch_size[[1L]] + seq_len(patch_size[[1L]]),
          j * patch_size[[2L]] + seq_len(patch_size[[2L]]),
          k * patch_size[[3L]] + seq_len(patch_size[[3L]])
        ])
      }
    }
  }
  counts
}

build_token_definitions <- function(atlas_96, maximum_tokens = 8L) {
  brain_counts <- atlas_patch_counts(atlas_96 > 0)
  coordinates <- data.table(
    grid_i = rep(0:(patch_grid[[1L]] - 1L), each = patch_grid[[2L]] * patch_grid[[3L]]),
    grid_j = rep(rep(0:(patch_grid[[2L]] - 1L), each = patch_grid[[3L]]),
                 times = patch_grid[[1L]]),
    grid_k = rep(0:(patch_grid[[3L]] - 1L), times = patch_grid[[1L]] * patch_grid[[2L]])
  )
  rows <- list()
  position <- 0L
  for (region in local_regions) {
    roi_mask <- atlas_96 %in% roi_labels[[region]]
    dim(roi_mask) <- dim(atlas_96)
    roi_counts <- atlas_patch_counts(roi_mask)
    candidates <- which(roi_counts > 0)
    selected <- sort(head(candidates[order(-roi_counts[candidates], candidates)], maximum_tokens))
    centroid <- colSums(
      as.matrix(coordinates[candidates]) * roi_counts[candidates]
    ) / sum(roi_counts[candidates])
    remote_candidates <- setdiff(which(brain_counts > 0), candidates)
    distances <- rowSums((
      as.matrix(coordinates[remote_candidates]) -
        matrix(centroid, nrow = length(remote_candidates), ncol = 3L, byrow = TRUE)
    )^2)
    remote <- sort(head(remote_candidates[order(-distances, remote_candidates)], length(selected)))
    for (group in c("roi", "remote")) {
      tokens <- if (group == "roi") selected else remote
      position <- position + 1L
      rows[[position]] <- data.table(
        target = region,
        group = group,
        token_index_zero_based = tokens - 1L,
        coordinates[tokens],
        roi_voxels_in_token = roi_counts[tokens],
        brain_voxels_in_token = brain_counts[tokens]
      )
    }
  }
  rbindlist(rows)
}

required_result_files <- c(
  "matfpca_model.rds", "analysis_manifest.csv", "subject_fpc_scores.csv",
  "product_fpc_eigenvalues.csv"
)
for (path in unlist(lapply(
  layer_specs$result_dir,
  function(directory) file.path(directory, required_result_files)
))) {
  if (!file.exists(path)) stop("Missing fitted MatrixFPCA output: ", path)
}

reference_manifest <- fread(
  file.path(layer_specs$result_dir[[1L]], "analysis_manifest.csv"),
  na.strings = c("", "NA")
)
reference_manifest[, subject_id := as.character(subject_id)]
if (uniqueN(reference_manifest$subject_id) < args$minimum_subjects) {
  stop(
    "Interpretation requires at least ", args$minimum_subjects,
    " subjects; refusing to analyze an undersized cohort."
  )
}
reference_keys <- reference_manifest[, paste(subject_id, image_id, sep = "|")]
for (directory in layer_specs$result_dir[-1L]) {
  candidate <- fread(file.path(directory, "analysis_manifest.csv"))
  keys <- candidate[, paste(subject_id, image_id, sep = "|")]
  if (!identical(reference_keys, keys)) stop("Layer manifests do not contain the same visits.")
}

first_nonmissing <- function(values) {
  available <- which(!is.na(values) & as.character(values) != "")
  if (length(available)) {
    values[[available[[1L]]]]
  } else {
    # Preserve the source storage type. data.table requires a grouping result
    # to keep the same type even when one subject has an entirely missing
    # integer, numeric, character, or logical covariate.
    values[NA_integer_][[1L]]
  }
}

derive_deformation_paths <- function(nifti_path) {
  name <- basename(nifti_path)
  if (!grepl("_Warped[.]nii[.]gz$", name)) {
    stop("Cannot derive log-Jacobian filename from ", name)
  }
  prefix <- sub("_Warped[.]nii[.]gz$", "", name)
  list(
    value = file.path(args$nonlinear_dir, paste0(prefix, "_LogJacobian_masked.nii.gz")),
    valid = file.path(
      args$nonlinear_dir,
      paste0(prefix, "_analysisValidMask_templateSpace.nii.gz")
    )
  )
}

target_cache <- file.path(args$output_dir, "visit_morphology_targets.csv")
atlas <- as.array(readNifti(args$atlas))
if (!identical(dim(atlas), c(193L, 229L, 193L))) {
  stop("Atlas must use the registered 193 x 229 x 193 grid.")
}
token_definitions <- build_token_definitions(crop_and_resize_atlas(atlas))
fwrite(token_definitions, file.path(args$output_dir, "token_definitions.csv"))
physical_targets <- NULL
if (args$cache_targets && file.exists(target_cache)) {
  candidate <- fread(target_cache)
  candidate_keys <- candidate[, paste(subject_id, image_id, sep = "|")]
  if (identical(reference_keys, candidate_keys) &&
      all(names(roi_labels) %in% names(candidate))) {
    physical_targets <- candidate
    log_message("Reusing regional morphology target cache.")
  }
}

if (is.null(physical_targets)) {
  log_message("Extracting atlas-defined regional log-Jacobian targets...")
  roi_masks <- lapply(roi_labels, function(labels) atlas %in% labels)
  physical_targets <- reference_manifest[, .(
    subject_id, image_id, visit_id, time_year, time_month, split
  )]
  for (region in names(roi_labels)) physical_targets[, (region) := NA_real_]
  coverage_rows <- vector("list", nrow(reference_manifest))
  for (row in seq_len(nrow(reference_manifest))) {
    paths <- derive_deformation_paths(reference_manifest$nifti_path[[row]])
    if (!file.exists(paths$value) || !file.exists(paths$valid)) {
      stop("Missing deformation input for image ", reference_manifest$image_id[[row]])
    }
    values <- as.array(readNifti(paths$value))
    valid <- as.array(readNifti(paths$valid)) > 0.5 & is.finite(values)
    if (!identical(dim(values), dim(atlas)) || !identical(dim(valid), dim(atlas))) {
      stop("Atlas/deformation shape mismatch for ", paths$value)
    }
    coverage <- numeric(length(roi_masks))
    for (region_position in seq_along(roi_masks)) {
      region <- names(roi_masks)[[region_position]]
      mask <- roi_masks[[region]] & valid
      if (sum(mask) < 100L) stop("Insufficient valid voxels for ", region)
      set(
        physical_targets, i = row, j = region,
        value = mean(values[mask], na.rm = TRUE)
      )
      coverage[[region_position]] <- sum(mask) / sum(roi_masks[[region]])
    }
    coverage_rows[[row]] <- data.table(
      image_id = reference_manifest$image_id[[row]],
      region = names(roi_masks), valid_fraction = coverage
    )
    if (row %% 20L == 0L || row == nrow(reference_manifest)) {
      log_message("Loaded deformation maps ", row, "/", nrow(reference_manifest), ".")
    }
  }
  fwrite(physical_targets, target_cache)
  coverage_table <- rbindlist(coverage_rows)
  fwrite(
    coverage_table[, .(
      mean_valid_fraction = mean(valid_fraction),
      minimum_valid_fraction = min(valid_fraction)
    ), by = region],
    file.path(args$output_dir, "morphology_target_coverage.csv")
  )
}

# Each subject's physical trajectory is summarized by its fitted baseline level
# and annual slope. This is descriptive: visit-to-template maps are not a
# dedicated within-subject longitudinal registration.
physical_long <- melt(
  physical_targets,
  id.vars = c("subject_id", "image_id", "time_year"),
  measure.vars = names(roi_labels),
  variable.name = "region", value.name = "log_jacobian"
)
physical_summary_long <- physical_long[, {
  fit <- lm(log_jacobian ~ time_year)
  list(
    morphology_level = unname(coef(fit)[[1L]]),
    morphology_slope = unname(coef(fit)[[2L]]),
    observed_mean = mean(log_jacobian),
    visits = .N,
    time_span_years = diff(range(time_year))
  )
}, by = .(subject_id, region)]

subject_morphology <- unique(physical_summary_long[, .(
  subject_id, visits, time_span_years
)], by = "subject_id")
for (region_name in names(roi_labels)) {
  values <- physical_summary_long[region == region_name]
  subject_morphology[values, on = "subject_id", paste0(region_name, "_level") :=
    i.morphology_level]
  subject_morphology[values, on = "subject_id", paste0(region_name, "_slope") :=
    i.morphology_slope]
}

clinical_columns <- intersect(
  c(
    "baseline_age", "sex", "education_years", "baseline_mmse",
    "baseline_cdrsb", "baseline_faq", "apoe_e4_count",
    "incident_ad_conversion", "baseline_dx", "outcome_type"
  ),
  names(reference_manifest)
)
clinical <- reference_manifest[, lapply(.SD, first_nonmissing),
                               by = subject_id, .SDcols = clinical_columns]
subject_morphology <- clinical[subject_morphology, on = "subject_id"]
subject_morphology[, sex_female := as.numeric(tolower(as.character(sex)) == "female")]
conversion_text <- tolower(as.character(subject_morphology$incident_ad_conversion))
subject_morphology[, incident_conversion_binary := fifelse(
  conversion_text %in% c("true", "1", "yes"), 1,
  fifelse(conversion_text %in% c("false", "0", "no"), 0, NA_real_)
)]
numeric_clinical <- intersect(
  c(
    "baseline_age", "education_years", "baseline_mmse", "baseline_cdrsb",
    "baseline_faq", "apoe_e4_count"
  ),
  names(subject_morphology)
)
for (column in numeric_clinical) {
  set(subject_morphology, j = column, value = as.numeric(subject_morphology[[column]]))
}
setorder(subject_morphology, subject_id)
fwrite(subject_morphology, file.path(args$output_dir, "subject_morphology_summary.csv"))

phenotype_table <- rbindlist(list(
  data.table(
    phenotype = unlist(lapply(names(roi_labels), function(region) {
      c(paste0(region, "_level"), paste0(region, "_slope"))
    })),
    phenotype_group = "regional_morphology"
  ),
  data.table(
    phenotype = intersect(
      c(
        "baseline_age", "sex_female", "education_years", "baseline_mmse",
        "baseline_cdrsb", "baseline_faq", "apoe_e4_count",
        "incident_conversion_binary"
      ),
      names(subject_morphology)
    ),
    phenotype_group = "demographic_clinical"
  )
))

standardize <- function(values) {
  values <- as.numeric(values)
  deviation <- sd(values, na.rm = TRUE)
  if (!is.finite(deviation) || deviation <= 0) return(rep(NA_real_, length(values)))
  (values - mean(values, na.rm = TRUE)) / deviation
}

rank_standardize <- function(values) {
  standardize(rank(as.numeric(values), ties.method = "average"))
}

score_reliability <- function(model, manifest, product_index) {
  total <- numeric(nrow(product_index))
  conditioning <- manifest[split != "validation_test"]
  rows_by_subject <- split(
    conditioning$time_index,
    factor(conditioning$subject_id, levels = model$subject_ids)
  )
  for (indices in rows_by_subject) {
    covariance <- model$covariance$K[indices, indices, drop = FALSE]
    decomposition <- eigen(0.5 * (covariance + t(covariance)), symmetric = TRUE)
    values <- pmax(decomposition$values, 0)
    projected_theta <- crossprod(
      decomposition$vectors,
      model$eig$funcs$Theta[indices, , drop = FALSE]
    )
    for (component in seq_len(nrow(product_index))) {
      row_mode <- product_index$row_component[[component]]
      column_mode <- product_index$column_component[[component]]
      time_mode <- product_index$time_component[[component]]
      spatial_signal <- model$eig$values$lam_phi[[row_mode]] *
        model$eig$values$lam_psi[[column_mode]]
      denominator <- values + max(
        1e-8, model$covariance$sig2 / spatial_signal
      )
      total[[component]] <- total[[component]] +
        sum(projected_theta[, time_mode]^2 / denominator) *
        model$eig$values$lam_theta[[time_mode]]
    }
  }
  total / length(rows_by_subject)
}

partial_association <- function(y, x, covariates) {
  complete <- is.finite(y) & is.finite(x)
  if (ncol(covariates)) {
    complete <- complete & apply(covariates, 1L, function(row) all(is.finite(row)))
  }
  y <- y[complete]
  x <- x[complete]
  covariates <- covariates[complete, , drop = FALSE]
  if (length(y) < 8L || length(unique(x)) < 2L || sd(y) <= 0 || sd(x) <= 0) {
    return(list(n = length(y), beta = NA_real_, partial_r2 = NA_real_, p = NA_real_))
  }
  # Rank-transform first so the adjusted analysis targets the same monotone
  # relationship as Spearman correlation and is less driven by outliers.
  yz <- rank_standardize(y)
  xz <- rank_standardize(x)
  covariates <- apply(covariates, 2L, rank_standardize)
  if (is.null(dim(covariates))) covariates <- matrix(covariates, ncol = 1L)
  reduced <- cbind(intercept = 1, covariates)
  full <- cbind(reduced, phenotype = xz)
  if (qr(full)$rank < ncol(full)) {
    return(list(n = length(y), beta = NA_real_, partial_r2 = NA_real_, p = NA_real_))
  }
  reduced_residual <- qr.resid(qr(reduced), yz)
  full_fit <- lm.fit(full, yz)
  reduced_sse <- sum(reduced_residual^2)
  full_sse <- sum(full_fit$residuals^2)
  df <- length(yz) - ncol(full)
  partial_r2 <- max(0, (reduced_sse - full_sse) / reduced_sse)
  statistic <- ((reduced_sse - full_sse) / 1) / (full_sse / df)
  list(
    n = length(y), beta = unname(full_fit$coefficients[[ncol(full)]]),
    partial_r2 = partial_r2,
    p = pf(max(0, statistic), 1, df, lower.tail = FALSE)
  )
}

covariates_for <- function(phenotype, data) {
  common <- c("baseline_age", "sex_female", "visits", "time_span_years")
  if (phenotype == "baseline_age") common <- setdiff(common, "baseline_age")
  if (phenotype == "sex_female") common <- setdiff(common, "sex_female")
  common <- intersect(common, names(data))
  common <- setdiff(common, phenotype)
  as.matrix(data[, ..common])
}

layer_objects <- list()
association_rows <- list()
reliability_rows <- list()
localization_rows <- list()
continuum_score_rows <- list()
continuum_diagnostic_rows <- list()
continuum_association_rows <- list()
continuum_sensitivity_rows <- list()
embedding_metric_rows <- list()
association_position <- 0L
localization_position <- 0L
for (layer_position in seq_len(nrow(layer_specs))) {
  layer <- layer_specs$layer[[layer_position]]
  result_dir <- layer_specs$result_dir[[layer_position]]
  log_message("Interpreting ", layer, " MatrixFPCA scores...")
  model <- readRDS(file.path(result_dir, "matfpca_model.rds"))
  manifest <- fread(file.path(result_dir, "analysis_manifest.csv"), na.strings = c("", "NA"))
  manifest[, subject_id := as.character(subject_id)]
  score_table <- fread(file.path(result_dir, "subject_fpc_scores.csv"))
  score_table[, subject_id := as.character(subject_id)]
  product_index <- fread(file.path(result_dir, "product_fpc_eigenvalues.csv"))
  fpc_columns <- grep("^FPC[0-9]+$", names(score_table), value = TRUE)
  component_count <- min(args$top_components, length(fpc_columns), nrow(product_index))
  fpc_columns <- fpc_columns[seq_len(component_count)]
  product_index <- product_index[seq_len(component_count)]
  product_index[, theoretical_variance_fraction := eigenvalue / sum(
    fread(file.path(result_dir, "product_fpc_eigenvalues.csv"))$eigenvalue
  )]

  product_index[, mean_score_reliability :=
    score_reliability(model, manifest, product_index)]
  scores <- as.matrix(score_table[, ..fpc_columns])
  colnames(scores) <- paste0("FPC", product_index$product_pc)
  subject_data <- subject_morphology[score_table, on = "subject_id"]
  if (anyNA(subject_data$subject_id)) stop("Score/morphology subject mismatch for ", layer)

  # Match the raw-image application conceptually: learn a one-dimensional
  # principal curve from leading raw product-FPC scores, without phenotype
  # labels. Age only orients the otherwise arbitrary direction.
  continuum_fit <- fit_score_principal_curve(
    scores, spline_df = 3, max_iterations = 20L, tolerance = 1e-4
  )
  continuum_fit <- orient_principal_curve(continuum_fit, subject_data$baseline_age)
  embedding_fit <- compare_score_embeddings(
    scores,
    neighbors = min(15L, nrow(scores) - 1L),
    min_dist = 0.20,
    tsne_perplexity = 8,
    seed = args$seed + 100L + layer_position
  )
  embedding_metric_rows[[layer_position]] <- copy(embedding_fit$metrics)[
    , layer := layer
  ]
  continuum_score_rows[[layer_position]] <- data.table(
    layer = layer,
    subject_id = subject_data$subject_id,
    latent_continuum = continuum_fit$coordinate,
    curve_residual = continuum_fit$residual_distance,
    baseline_age = subject_data$baseline_age,
    incident_ad_conversion = subject_data$incident_ad_conversion
  )
  continuum_diagnostic_rows[[layer_position]] <- data.table(
    layer = layer,
    score_dimensions = ncol(scores),
    selected_product_cumulative_pve = product_index$cumulative_pve[[component_count]],
    principal_curve_variance_explained = continuum_fit$variance_explained,
    linear_axis_variance_explained = continuum_fit$linear_variance_explained,
    nonlinear_improvement = continuum_fit$variance_explained -
      continuum_fit$linear_variance_explained,
    curve_vs_linear_spearman = suppressWarnings(cor(
      continuum_fit$coordinate, continuum_fit$linear_coordinate,
      method = "spearman"
    )),
    coordinate_vs_age_spearman = suppressWarnings(cor(
      continuum_fit$coordinate, subject_data$baseline_age,
      method = "spearman", use = "complete.obs"
    )),
    iterations = continuum_fit$iterations
  )
  for (phenotype_position in seq_len(nrow(phenotype_table))) {
    phenotype <- phenotype_table$phenotype[[phenotype_position]]
    target <- as.numeric(subject_data[[phenotype]])
    complete <- is.finite(target) & is.finite(continuum_fit$coordinate)
    test <- if (sum(complete) >= 4L && length(unique(target[complete])) > 1L) {
      suppressWarnings(cor.test(
        continuum_fit$coordinate[complete], target[complete],
        method = "spearman", exact = FALSE
      ))
    } else NULL
    partial <- partial_association(
      continuum_fit$coordinate, target, covariates_for(phenotype, subject_data)
    )
    continuum_association_rows[[length(continuum_association_rows) + 1L]] <- data.table(
      layer = layer,
      phenotype = phenotype,
      phenotype_group = phenotype_table$phenotype_group[[phenotype_position]],
      n = sum(complete),
      spearman_rho = if (is.null(test)) NA_real_ else unname(test$estimate),
      spearman_p = if (is.null(test)) NA_real_ else test$p.value,
      adjusted_standardized_beta = partial$beta,
      adjusted_partial_r2 = partial$partial_r2,
      adjusted_p = partial$p
    )
  }
  sensitivity_dimensions <- sort(unique(pmin(c(3L, 5L, 10L, 20L), ncol(scores))))
  for (dimension in sensitivity_dimensions) {
    sensitivity_fit <- fit_score_principal_curve(
      scores[, seq_len(dimension), drop = FALSE],
      spline_df = 3, max_iterations = 20L, tolerance = 1e-4
    )
    sensitivity_fit <- orient_principal_curve(
      sensitivity_fit, subject_data$baseline_age
    )
    continuum_sensitivity_rows[[length(continuum_sensitivity_rows) + 1L]] <- data.table(
      layer = layer,
      score_dimensions = dimension,
      coordinate_vs_full_spearman = suppressWarnings(cor(
        sensitivity_fit$coordinate, continuum_fit$coordinate,
        method = "spearman"
      )),
      coordinate_vs_full_abs_spearman = abs(suppressWarnings(cor(
        sensitivity_fit$coordinate, continuum_fit$coordinate,
        method = "spearman"
      ))),
      variance_explained = sensitivity_fit$variance_explained,
      coordinate_vs_age_spearman = suppressWarnings(cor(
        sensitivity_fit$coordinate, subject_data$baseline_age,
        method = "spearman", use = "complete.obs"
      ))
    )
  }

  empirical_sd <- apply(scores, 2L, sd)
  product_index[, `:=`(
    layer = layer,
    empirical_score_sd = empirical_sd,
    empirical_variance_to_eigenvalue = empirical_sd^2 / eigenvalue
  )]
  reliability_rows[[layer_position]] <- copy(product_index)

  for (component_position in seq_len(component_count)) {
    score <- scores[, component_position]
    for (phenotype_position in seq_len(nrow(phenotype_table))) {
      phenotype <- phenotype_table$phenotype[[phenotype_position]]
      x <- as.numeric(subject_data[[phenotype]])
      complete <- is.finite(score) & is.finite(x)
      test <- if (sum(complete) >= 4L && length(unique(x[complete])) > 1L) {
        suppressWarnings(cor.test(score[complete], x[complete], method = "spearman", exact = FALSE))
      } else NULL
      partial <- partial_association(
        score, x, covariates_for(phenotype, subject_data)
      )
      association_position <- association_position + 1L
      association_rows[[association_position]] <- data.table(
        layer = layer,
        product_pc = product_index$product_pc[[component_position]],
        row_component = product_index$row_component[[component_position]],
        column_component = product_index$column_component[[component_position]],
        time_component = product_index$time_component[[component_position]],
        theoretical_variance_fraction =
          product_index$theoretical_variance_fraction[[component_position]],
        cumulative_pve = product_index$cumulative_pve[[component_position]],
        mean_score_reliability = product_index$mean_score_reliability[[component_position]],
        phenotype = phenotype,
        phenotype_group = phenotype_table$phenotype_group[[phenotype_position]],
        n = sum(complete),
        spearman_rho = if (is.null(test)) NA_real_ else unname(test$estimate),
        spearman_p = if (is.null(test)) NA_real_ else test$p.value,
        adjusted_standardized_beta = partial$beta,
        adjusted_partial_r2 = partial$partial_r2,
        adjusted_p = partial$p
      )
    }
  }

  for (component_position in seq_len(component_count)) {
    token_mode <- model$eig$funcs$Phi[, product_index$row_component[[component_position]]]
    for (region in local_regions) {
      roi_tokens <- token_definitions[target == region & group == "roi"]$token_index_zero_based + 1L
      remote_tokens <- token_definitions[target == region & group == "remote"]$token_index_zero_based + 1L
      localization_position <- localization_position + 1L
      localization_rows[[localization_position]] <- data.table(
        layer = layer,
        product_pc = product_index$product_pc[[component_position]],
        row_component = product_index$row_component[[component_position]],
        region = region,
        roi_mean_absolute_loading = mean(abs(token_mode[roi_tokens])),
        remote_mean_absolute_loading = mean(abs(token_mode[remote_tokens])),
        roi_to_remote_loading_ratio = mean(abs(token_mode[roi_tokens])) /
          mean(abs(token_mode[remote_tokens]))
      )
    }
  }

  layer_objects[[layer]] <- list(
    model = model,
    manifest = manifest,
    scores = scores,
    score_table = score_table,
    subject_data = subject_data,
    product_index = product_index,
    continuum_fit = continuum_fit,
    embedding_fit = embedding_fit
  )
}

associations <- rbindlist(association_rows)
associations[, spearman_q := p.adjust(spearman_p, method = "BH"),
             by = .(layer, phenotype)]
associations[, adjusted_q := p.adjust(adjusted_p, method = "BH"),
             by = .(layer, phenotype)]
associations[, spearman_global_q := p.adjust(spearman_p, method = "BH"), by = layer]
associations[, adjusted_global_q := p.adjust(adjusted_p, method = "BH"), by = layer]
fwrite(associations, file.path(args$output_dir, "eigenmode_phenotype_associations.csv"))
fwrite(rbindlist(reliability_rows), file.path(args$output_dir, "score_reliability.csv"))
localization <- rbindlist(localization_rows)
fwrite(localization, file.path(args$output_dir, "token_roi_localization.csv"))
continuum_scores <- rbindlist(continuum_score_rows)
continuum_diagnostics <- rbindlist(continuum_diagnostic_rows)
continuum_associations <- rbindlist(continuum_association_rows)
continuum_associations[, spearman_q := p.adjust(spearman_p, method = "BH"), by = layer]
continuum_associations[, adjusted_q := p.adjust(adjusted_p, method = "BH"), by = layer]
fwrite(continuum_scores, file.path(args$output_dir, "latent_continuum_scores.csv"))
fwrite(continuum_diagnostics, file.path(args$output_dir, "latent_continuum_diagnostics.csv"))
fwrite(continuum_associations, file.path(args$output_dir, "latent_continuum_associations.csv"))
fwrite(rbindlist(continuum_sensitivity_rows),
       file.path(args$output_dir, "latent_continuum_dimension_sensitivity.csv"))
fwrite(
  rbindlist(embedding_metric_rows),
  file.path(args$output_dir, "latent_score_embedding_quality.csv")
)

# Spatial localization and anatomy association answer different questions. A
# mode can order subjects by ventricular size while using a nearly global token
# pattern. Keep a separate, explicitly exploratory list of modes whose token
# loading is concentrated in an atlas ROI relative to matched distant tokens.
localization_summary <- localization[
  order(-roi_to_remote_loading_ratio), .SD[1L],
  by = .(layer, product_pc, row_component)
][, .(
  layer, product_pc, row_component,
  most_localized_region = region,
  maximum_roi_to_remote_ratio = roi_to_remote_loading_ratio,
  localized_roi_mean_absolute_loading = roi_mean_absolute_loading,
  localized_remote_mean_absolute_loading = remote_mean_absolute_loading
)]
localized_pool <- merge(
  associations[
    phenotype_group == "regional_morphology" &
      mean_score_reliability >= args$min_reliability &
      is.finite(adjusted_partial_r2)
  ],
  localization_summary,
  by = c("layer", "product_pc", "row_component")
)
localized_pool[, phenotype_region := sub("_(level|slope)$", "", phenotype)]
localized_pool[, localization_matches_phenotype :=
  phenotype_region == most_localized_region]
localized_pool[, localized_selection_score :=
  abs(spearman_rho) * sqrt(adjusted_partial_r2) *
    sqrt(theoretical_variance_fraction) *
    log1p(maximum_roi_to_remote_ratio)]
localized_candidates <- localized_pool[
  maximum_roi_to_remote_ratio >= 1.5 & is.finite(localized_selection_score)
][
  order(-localized_selection_score, spearman_global_q, phenotype),
  .SD[1L], by = .(layer, product_pc)
][
  order(-localized_selection_score, spearman_global_q, product_pc),
  head(.SD, 3L), by = layer
]
fwrite(
  localized_candidates,
  file.path(args$output_dir, "localized_eigenmode_candidates.csv")
)

# Select one explicitly post-hoc anatomy-linked coordinate per layer. Selection
# is constrained to leading, estimable product modes and is reported as
# exploratory to avoid presenting a supervised winner as an unsupervised fact.
eligible_anatomy <- associations[
  phenotype_group == "regional_morphology" &
    mean_score_reliability >= args$min_reliability &
    is.finite(adjusted_partial_r2)
]
eligible_anatomy[, selection_score :=
  abs(spearman_rho) * sqrt(adjusted_partial_r2) *
    sqrt(theoretical_variance_fraction)]
selected_coordinates <- eligible_anatomy[
  order(-selection_score, spearman_global_q, product_pc), .SD[1L], by = layer
]
if (nrow(selected_coordinates) != nrow(layer_specs)) {
  stop("Could not select one reliable morphology coordinate for every layer.")
}

bootstrap_coordinate <- function(score, target, seed) {
  complete <- is.finite(score) & is.finite(target)
  score <- score[complete]
  target <- target[complete]
  observed <- suppressWarnings(cor(score, target, method = "spearman"))
  set.seed(seed)
  draws <- replicate(args$bootstrap_replicates, {
    sample_index <- sample(seq_along(score), replace = TRUE)
    suppressWarnings(cor(score[sample_index], target[sample_index], method = "spearman"))
  })
  leave_one_out <- vapply(seq_along(score), function(index) {
    suppressWarnings(cor(score[-index], target[-index], method = "spearman"))
  }, numeric(1L))
  list(
    rho = observed,
    ci_low = unname(quantile(draws, 0.025, na.rm = TRUE)),
    ci_high = unname(quantile(draws, 0.975, na.rm = TRUE)),
    loo_min = min(leave_one_out, na.rm = TRUE),
    loo_max = max(leave_one_out, na.rm = TRUE)
  )
}

coordinate_score_rows <- list()
for (selection_position in seq_len(nrow(selected_coordinates))) {
  selection <- selected_coordinates[selection_position]
  object <- layer_objects[[selection$layer]]
  component_position <- match(selection$product_pc, object$product_index$product_pc)
  score <- object$scores[, component_position]
  target <- as.numeric(object$subject_data[[selection$phenotype]])
  stability <- bootstrap_coordinate(score, target, args$seed + selection_position)
  orientation <- ifelse(stability$rho >= 0, 1, -1)
  selected_coordinates[selection_position, `:=`(
    orientation_sign = orientation,
    oriented_spearman_rho = orientation * stability$rho,
    bootstrap_ci_low = orientation * ifelse(orientation > 0, stability$ci_low, stability$ci_high),
    bootstrap_ci_high = orientation * ifelse(orientation > 0, stability$ci_high, stability$ci_low),
    leave_one_out_abs_rho_min = min(abs(c(stability$loo_min, stability$loo_max))),
    leave_one_out_rho_min = orientation * ifelse(orientation > 0, stability$loo_min, stability$loo_max),
    leave_one_out_rho_max = orientation * ifelse(orientation > 0, stability$loo_max, stability$loo_min)
  )]
  coordinate_score_rows[[selection_position]] <- data.table(
    layer = selection$layer,
    product_pc = selection$product_pc,
    phenotype = selection$phenotype,
    subject_id = object$subject_data$subject_id,
    oriented_score_z = standardize(score * orientation),
    target_value = target,
    target_z = standardize(target),
    baseline_age = object$subject_data$baseline_age,
    sex = object$subject_data$sex,
    incident_ad_conversion = object$subject_data$incident_ad_conversion
  )
}
coordinate_scores <- rbindlist(coordinate_score_rows)
fwrite(selected_coordinates, file.path(args$output_dir, "selected_morphology_coordinates.csv"))
fwrite(coordinate_scores, file.path(args$output_dir, "morphology_coordinate_scores.csv"))

# Top demographic/clinical annotations are kept separate from morphology-based
# coordinate selection.
selected_demographic <- associations[
  phenotype_group == "demographic_clinical" &
    mean_score_reliability >= args$min_reliability &
    is.finite(adjusted_partial_r2)
][order(-adjusted_partial_r2), .SD[1L], by = .(layer, phenotype)]
fwrite(selected_demographic, file.path(args$output_dir, "top_demographic_associations.csv"))
fwrite(
  associations[phenotype_group == "regional_morphology"][
    order(layer, spearman_global_q, -abs(spearman_rho))
  ],
  file.path(args$output_dir, "ranked_anatomical_associations.csv")
)

temporal_summary <- associations[
  is.finite(adjusted_partial_r2), .(
    tested_associations = .N,
    maximum_partial_r2 = max(adjusted_partial_r2),
    median_partial_r2 = median(adjusted_partial_r2),
    phenotype_wise_adjusted_fdr_0.10 = sum(adjusted_q < 0.10, na.rm = TRUE),
    layer_wide_spearman_fdr_0.10 = sum(spearman_global_q < 0.10, na.rm = TRUE),
    layer_wide_adjusted_fdr_0.10 = sum(adjusted_global_q < 0.10, na.rm = TRUE)
  ), by = .(layer, time_component, phenotype_group)
]
fwrite(temporal_summary, file.path(args$output_dir, "temporal_mode_interpretation.csv"))

# Cross-layer score concordance: high values mean two layers recover a similar
# subject ordering even when their component numbers differ.
concordance_rows <- list()
concordance_position <- 0L
for (left_position in seq_len(length(layer_objects) - 1L)) {
  for (right_position in seq.int(left_position + 1L, length(layer_objects))) {
    left_name <- names(layer_objects)[[left_position]]
    right_name <- names(layer_objects)[[right_position]]
    left <- layer_objects[[left_name]]
    right <- layer_objects[[right_name]]
    right_order <- match(left$score_table$subject_id, right$score_table$subject_id)
    correlation <- cor(left$scores, right$scores[right_order, , drop = FALSE],
                       method = "spearman")
    for (left_component in seq_len(nrow(correlation))) {
      best <- which.max(abs(correlation[left_component, ]))
      concordance_position <- concordance_position + 1L
      concordance_rows[[concordance_position]] <- data.table(
        left_layer = left_name,
        left_product_pc = left$product_index$product_pc[[left_component]],
        right_layer = right_name,
        matched_right_product_pc = right$product_index$product_pc[[best]],
        spearman_rho = correlation[left_component, best],
        absolute_rho = abs(correlation[left_component, best])
      )
    }
  }
}
concordance <- rbindlist(concordance_rows)
fwrite(concordance, file.path(args$output_dir, "cross_layer_score_concordance.csv"))

continuum_wide <- dcast(
  continuum_scores, subject_id ~ layer, value.var = "latent_continuum"
)
continuum_layer_concordance <- rbindlist(lapply(
  combn(layer_specs$layer, 2L, simplify = FALSE),
  function(pair) data.table(
    left_layer = pair[[1L]], right_layer = pair[[2L]],
    spearman_rho = suppressWarnings(cor(
      continuum_wide[[pair[[1L]]]], continuum_wide[[pair[[2L]]]],
      method = "spearman"
    )),
    absolute_spearman_rho = abs(suppressWarnings(cor(
      continuum_wide[[pair[[1L]]]], continuum_wide[[pair[[2L]]]],
      method = "spearman"
    )))
  )
))
fwrite(
  continuum_layer_concordance,
  file.path(args$output_dir, "latent_continuum_layer_concordance.csv")
)

# Choose distinct subjects nearest fixed continuum percentiles, then retrieve
# each subject's actual observed baseline MRI. This is not decoding: these are
# the source images whose longitudinal visits generated the BrainIAC features.
continuum_landmark_percentiles <- seq(0.05, 0.95, length.out = 6L)
landmark_rows <- list()
embedding_rows <- list()
for (layer_position in seq_len(nrow(layer_specs))) {
  layer_name <- layer_specs$layer[[layer_position]]
  object <- layer_objects[[layer_name]]
  coordinate <- object$continuum_fit$coordinate
  positions <- as.numeric(quantile(
    coordinate, probs = continuum_landmark_percentiles, names = FALSE
  ))
  chosen <- integer()
  selected <- vapply(positions, function(target) {
    candidates <- setdiff(seq_along(coordinate), chosen)
    ordering <- order(
      abs(coordinate[candidates] - target),
      object$continuum_fit$residual_distance[candidates],
      object$subject_data$subject_id[candidates]
    )
    index <- candidates[[ordering[[1L]]]]
    chosen <<- c(chosen, index)
    index
  }, integer(1L))
  layer_landmarks <- rbindlist(lapply(seq_along(selected), function(index) {
    subject_position <- selected[[index]]
    subject <- object$subject_data$subject_id[[subject_position]]
    visits <- object$manifest[subject_id == subject]
    baseline_visit <- visits[which.min(abs(time_year))]
    image_name <- sub("[.]nii[.]gz$", ".png", basename(baseline_visit$nifti_path))
    data.table(
      layer = layer_name,
      target_percentile = continuum_landmark_percentiles[[index]],
      target_coordinate = positions[[index]],
      subject_id = subject,
      subject_coordinate = coordinate[[subject_position]],
      curve_residual = object$continuum_fit$residual_distance[[subject_position]],
      baseline_age = object$subject_data$baseline_age[[subject_position]],
      image_id = baseline_visit$image_id,
      image_time_year = baseline_visit$time_year,
      source_nifti_path = baseline_visit$nifti_path,
      axial_png = normalizePath(file.path(
        "data", "ADNI", "affine_axial_slice094", image_name
      )),
      coronal_png = normalizePath(file.path(
        "data", "ADNI", "affine_coronal_slice110_upright", image_name
      ))
    )
  }))
  landmark_rows[[layer_position]] <- layer_landmarks
  embeddings <- object$embedding_fit$embeddings
  embedding_rows[[layer_position]] <- data.table(
    layer = layer_name,
    subject_id = object$subject_data$subject_id,
    latent_continuum = coordinate,
    latent_continuum_percentile =
      (frank(coordinate, ties.method = "average") - 0.5) / length(coordinate),
    baseline_age = object$subject_data$baseline_age,
    PCA1 = embeddings$PCA[, 1L], PCA2 = embeddings$PCA[, 2L],
    UMAP1 = embeddings$UMAP[, 1L], UMAP2 = embeddings$UMAP[, 2L],
    TSNE1 = embeddings$`t-SNE`[, 1L], TSNE2 = embeddings$`t-SNE`[, 2L]
  )
}
continuum_landmarks <- rbindlist(landmark_rows)
continuum_embeddings <- rbindlist(embedding_rows)
fwrite(
  continuum_landmarks,
  file.path(args$output_dir, "latent_continuum_observed_mri_landmarks.csv")
)
fwrite(
  continuum_embeddings,
  file.path(args$output_dir, "latent_score_embeddings.csv")
)

# -------------------------------------------------------------------------
# Figures

group_color <- function(values) {
  text <- tolower(as.character(values))
  ifelse(text %in% c("true", "1", "yes"), "#B2182B", "#2166AC")
}

continuum_palette <- hcl.colors(256L, "Viridis")
continuum_color <- function(coordinate) {
  index <- 1L + floor(pmin(1, pmax(0, coordinate)) * 255L)
  continuum_palette[index]
}

draw_coordinate_age <- function(coordinate, age, title) {
  colors <- continuum_color(coordinate)
  rho <- suppressWarnings(cor(
    coordinate, age, method = "spearman", use = "complete.obs"
  ))
  plot(
    coordinate, age, pch = 19, cex = 0.85,
    col = grDevices::adjustcolor(colors, 0.72),
    xlab = "Label-free latent continuum coordinate",
    ylab = "Baseline age (years)", main = title
  )
  lines(lowess(coordinate, age, f = 0.8), col = "#333333", lwd = 3)
  mtext(sprintf("Spearman rho = %.3f; n = %d", rho, length(age)),
        side = 3L, line = 0.35, cex = 0.82)
  grid(col = "grey90")
}

for (layer_name in layer_specs$layer) {
  values <- continuum_embeddings[layer == layer_name]
  png(
    file.path(args$output_dir, paste0("latent_coordinate_vs_age_", layer_name, ".png")),
    width = 1050L, height = 780L, res = 140L
  )
  par(mar = c(5, 4.8, 4.8, 1))
  draw_coordinate_age(
    values$latent_continuum, values$baseline_age,
    paste("Age along the", gsub("_", " ", layer_name), "latent continuum")
  )
  dev.off()

  png(
    file.path(
      args$output_dir,
      paste0("fpc_score_embedding_comparison_by_continuum_", layer_name, ".png")
    ),
    width = 1800L, height = 620L, res = 140L
  )
  par(mfrow = c(1L, 3L), mar = c(4.2, 4.2, 4.2, 1), oma = c(0, 0, 2.5, 0))
  embedding_names <- c("PCA", "UMAP", "t-SNE")
  embedding_columns <- list(
    PCA = c("PCA1", "PCA2"),
    UMAP = c("UMAP1", "UMAP2"),
    `t-SNE` = c("TSNE1", "TSNE2")
  )
  for (method in embedding_names) {
    columns <- embedding_columns[[method]]
    plot(
      values[[columns[[1L]]]], values[[columns[[2L]]]],
      pch = 19, cex = 0.95,
      col = grDevices::adjustcolor(
        continuum_color(values$latent_continuum_percentile), 0.82
      ),
      xlab = paste(method, "1"), ylab = paste(method, "2"), main = method
    )
    legend_values <- seq(0, 1, length.out = 5L)
    legend(
      "topright", legend = sprintf("%.2f", legend_values),
      col = continuum_color(legend_values), pch = 19, bty = "n", cex = 0.72,
      title = "percentile"
    )
  }
  mtext(
    sprintf(
      "%s: principal curve through %d scores; color is continuum percentile",
      gsub("_", " ", layer_name),
      layer_objects[[layer_name]]$product_index[, .N]
    ),
    side = 3L, outer = TRUE, line = 0.6, cex = 0.95
  )
  dev.off()
}

draw_observed_mri <- function(path, title) {
  image <- png::readPNG(path)
  plot.new()
  rasterImage(image, 0, 0, 1, 1, interpolate = TRUE)
  title(main = title, cex.main = 0.90)
  box()
}

for (layer_name in layer_specs$layer) {
  landmarks <- continuum_landmarks[layer == layer_name][order(target_percentile)]
  png(
    file.path(args$output_dir, paste0("observed_mri_continuum_", layer_name, ".png")),
    width = 390L * nrow(landmarks), height = 410L,
    res = 120L, bg = "white"
  )
  par(mfrow = c(1L, nrow(landmarks)), mar = c(1, 1, 3.2, 1))
  for (index in seq_len(nrow(landmarks))) {
    draw_observed_mri(
      landmarks$axial_png[[index]],
      sprintf(
        "Percentile %.2f\narc %.2f; age %.1f",
        landmarks$target_percentile[[index]],
        landmarks$subject_coordinate[[index]], landmarks$baseline_age[[index]]
      )
    )
  }
  dev.off()
}

# Direct visual comparison with the raw-image MatrixFPCA continuum. All four
# rows show actual observed baseline scans; none is a decoded reconstruction.
raw_score_path <- file.path(
  "application", "adni", "results", "fpc_scores_embeddings_pcs30.csv"
)
raw_landmark_path <- file.path(
  "application", "adni", "results", "latent_coordinate_landmarks_pcs30.csv"
)
if (file.exists(raw_score_path) && file.exists(raw_landmark_path)) {
  raw_scores <- fread(raw_score_path, select = c(
    "PTID", "baseline_age", "latent_coordinate"
  ))
  png(
    file.path(args$output_dir, "raw_vs_brainiac_latent_coordinate_vs_age.png"),
    width = 1500L, height = 1260L, res = 150L
  )
  par(mfrow = c(2L, 2L), mar = c(4.5, 4.5, 4.4, 1.0), oma = c(0, 0, 2.2, 0))
  draw_coordinate_age(
    raw_scores$latent_coordinate, raw_scores$baseline_age,
    "Raw-image MatrixFPCA"
  )
  for (layer_name in layer_specs$layer) {
    values <- continuum_embeddings[layer == layer_name]
    draw_coordinate_age(
      values$latent_continuum, values$baseline_age,
      paste("BrainIAC", gsub("_", " ", layer_name))
    )
  }
  mtext(
    "Raw-image and BrainIAC MatrixFPCA continua: age comparison",
    side = 3L, outer = TRUE, cex = 1.12, font = 2
  )
  dev.off()

  raw_landmarks <- fread(raw_landmark_path)[order(target_percentile)]
  comparison_landmarks <- rbindlist(list(
    data.table(
      layer = "raw_image",
      target_percentile = raw_landmarks$target_percentile,
      subject_coordinate = raw_landmarks$subject_coordinate,
      baseline_age = raw_landmarks$baseline_age,
      axial_png = raw_landmarks$source_image_path
    ),
    continuum_landmarks[, .(
      layer, target_percentile, subject_coordinate, baseline_age, axial_png
    )]
  ), use.names = TRUE)
  comparison_layers <- c("raw_image", layer_specs$layer)
  png(
    file.path(args$output_dir, "raw_vs_brainiac_observed_mri_continua.png"),
    width = 390L * length(continuum_landmark_percentiles),
    height = 400L * length(comparison_layers), res = 120L, bg = "white"
  )
  par(
    mfrow = c(length(comparison_layers), length(continuum_landmark_percentiles)),
    mar = c(0.7, 0.7, 3.0, 0.7)
  )
  for (layer_name in comparison_layers) {
    landmarks <- comparison_landmarks[layer == layer_name][order(target_percentile)]
    for (index in seq_len(nrow(landmarks))) {
      draw_observed_mri(
        landmarks$axial_png[[index]],
        sprintf(
          "%s  p=%.2f\narc %.2f; age %.1f",
          ifelse(index == 1L, gsub("_", " ", layer_name), "coord"),
          landmarks$target_percentile[[index]],
          landmarks$subject_coordinate[[index]], landmarks$baseline_age[[index]]
        )
      )
    }
  }
  dev.off()
}

png(file.path(args$output_dir, "continuous_morphology_coordinates.png"),
    width = 1650, height = 570, res = 150)
par(mfrow = c(1L, 3L), mar = c(4.5, 4.5, 3.5, 1.0), oma = c(0, 0, 2.0, 0))
for (layer_name in layer_specs$layer) {
  selection <- selected_coordinates[layer == layer_name]
  values <- coordinate_scores[layer == layer_name]
  colors <- group_color(values$incident_ad_conversion)
  plot(
    values$oriented_score_z, values$target_z,
    pch = 19, col = grDevices::adjustcolor(colors, 0.78),
    xlab = sprintf("Oriented %s product FPC%d (z)", layer_name, selection$product_pc),
    ylab = paste(gsub("_", " ", selection$phenotype), "(z)"),
    main = sprintf(
      "%s: rho %.2f [%.2f, %.2f]",
      layer_name, selection$oriented_spearman_rho,
      selection$bootstrap_ci_low, selection$bootstrap_ci_high
    )
  )
  abline(lm(target_z ~ oriented_score_z, data = values), col = "grey25", lwd = 2)
  grid(col = "grey88")
}
mtext(
  "Exploratory continuous morphology coordinates (blue non-converter; red converter)",
  side = 3, outer = TRUE, cex = 1.15, font = 2
)
dev.off()

continuum_plot_targets <- continuum_associations[
  phenotype_group == "regional_morphology" & is.finite(spearman_rho)
][order(-abs(spearman_rho)), .SD[1L], by = layer]
fwrite(
  continuum_plot_targets,
  file.path(args$output_dir, "latent_continuum_top_morphology.csv")
)
png(file.path(args$output_dir, "latent_continuum_interpretation.png"),
    width = 1650, height = 1020, res = 150)
par(mfrow = c(2L, 3L), mar = c(4.4, 4.5, 3.4, 1.0), oma = c(0, 0, 2.0, 0))
for (layer_name in layer_specs$layer) {
  object <- layer_objects[[layer_name]]
  coordinate <- object$continuum_fit$coordinate
  selection <- continuum_plot_targets[layer == layer_name]
  target <- as.numeric(object$subject_data[[selection$phenotype]])
  colors <- group_color(object$subject_data$incident_ad_conversion)
  plot(
    coordinate, standardize(target), pch = 19,
    col = grDevices::adjustcolor(colors, 0.78),
    xlab = paste(layer_name, "unsupervised continuum"),
    ylab = paste(gsub("_", " ", selection$phenotype), "(z)"),
    main = sprintf("rho %.2f; q %.3f", selection$spearman_rho, selection$spearman_q)
  )
  lines(lowess(coordinate, standardize(target)), col = "grey20", lwd = 2)
  grid(col = "grey88")
}
for (layer_name in layer_specs$layer) {
  object <- layer_objects[[layer_name]]
  coordinate <- object$continuum_fit$coordinate
  age <- object$subject_data$baseline_age
  colors <- group_color(object$subject_data$incident_ad_conversion)
  rho <- suppressWarnings(cor(coordinate, age, method = "spearman"))
  plot(
    coordinate, age, pch = 19, col = grDevices::adjustcolor(colors, 0.78),
    xlab = paste(layer_name, "unsupervised continuum"),
    ylab = "Baseline age", main = sprintf("Age rho %.2f", rho)
  )
  lines(lowess(coordinate, age), col = "grey20", lwd = 2)
  grid(col = "grey88")
}
mtext(
  "Unsupervised score-space continua: morphology (top) and age orientation (bottom)",
  side = 3, outer = TRUE, cex = 1.12, font = 2
)
dev.off()

heatmap_phenotypes <- c(
  "whole_brain_level", "white_matter_level", "hippocampus_level",
  "lateral_ventricles_level", "amygdala_level", "thalamus_level",
  "hippocampus_slope", "lateral_ventricles_slope",
  "baseline_age", "baseline_mmse", "baseline_cdrsb", "baseline_faq",
  "incident_conversion_binary"
)
heatmap_phenotypes <- intersect(heatmap_phenotypes, phenotype_table$phenotype)
png(file.path(args$output_dir, "eigenmode_interpretation_heatmap.png"),
    width = 1780, height = 1450, res = 150)
par(mfrow = c(3L, 1L), mar = c(6.8, 4.8, 2.8, 1.2))
palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101L)
for (layer_name in layer_specs$layer) {
  display_components <- seq_len(min(10L, nrow(layer_objects[[layer_name]]$product_index)))
  values <- associations[
    layer == layer_name & product_pc %in% display_components &
      phenotype %in% heatmap_phenotypes
  ]
  matrix_values <- dcast(values, product_pc ~ phenotype, value.var = "spearman_rho")
  component_labels <- paste0("FPC", matrix_values$product_pc)
  matrix_values <- as.matrix(matrix_values[, -"product_pc"])
  image(
    x = seq_len(ncol(matrix_values)), y = seq_len(nrow(matrix_values)),
    z = t(matrix_values), zlim = c(-1, 1), col = palette,
    axes = FALSE, xlab = "", ylab = "", main = paste(layer_name, "leading product FPCs")
  )
  axis(1, at = seq_len(ncol(matrix_values)), labels = gsub("_", " ", colnames(matrix_values)),
       las = 2, cex.axis = 0.62)
  axis(2, at = seq_len(nrow(matrix_values)), labels = component_labels, las = 1, cex.axis = 0.75)
  box()
}
dev.off()

plot_midplane <- function(values, title) {
  cube <- token_vector_to_cube(values)
  planes <- list(
    sagittal = cube[3L, , ], coronal = cube[, 3L, ], axial = cube[, , 3L]
  )
  limit <- max(abs(unlist(planes)))
  for (name in names(planes)) {
    image(
      t(planes[[name]])[, 6:1, drop = FALSE], axes = FALSE, useRaster = FALSE,
      col = hcl.colors(80, "Blue-Red 3"), zlim = c(-limit, limit),
      main = paste(title, name)
    )
    box()
  }
}

png(file.path(args$output_dir, "selected_eigenmode_token_time_maps.png"),
    width = 1700, height = 1180, res = 150)
layout(matrix(seq_len(12L), nrow = 3L, byrow = TRUE), widths = c(1, 1, 1, 1.25))
par(mar = c(3.5, 3.5, 2.6, 0.8))
for (layer_name in layer_specs$layer) {
  selection <- selected_coordinates[layer == layer_name]
  object <- layer_objects[[layer_name]]
  orientation <- selection$orientation_sign
  token_values <- orientation * object$model$eig$funcs$Phi[, selection$row_component]
  plot_midplane(token_values, sprintf("%s FPC%d", layer_name, selection$product_pc))
  theta <- orientation * object$model$eig$funcs$Theta[, selection$time_component]
  plot(
    object$model$t_grid, theta, type = "l", lwd = 2.5, col = "#762A83",
    xlab = "Years", ylab = "Oriented temporal loading",
    main = sprintf(
      "%s FPC%d: T%d\n%s",
      layer_name, selection$product_pc, selection$time_component,
      gsub("_", " ", selection$phenotype)
    )
  )
  abline(h = 0, col = "grey70", lty = 2)
}
dev.off()

top_localized_candidates <- localized_candidates[
  localization_matches_phenotype == TRUE
][
  order(-localized_selection_score), .SD[1L], by = layer
]
if (nrow(top_localized_candidates)) {
  png(file.path(args$output_dir, "localized_eigenmode_candidates.png"),
      width = 1700, height = 520 * nrow(top_localized_candidates), res = 150)
  layout(
    matrix(seq_len(4L * nrow(top_localized_candidates)),
           nrow = nrow(top_localized_candidates), byrow = TRUE),
    widths = c(1, 1, 1, 1.35)
  )
  par(mar = c(4.5, 3.5, 3.2, 0.8))
  for (candidate_position in seq_len(nrow(top_localized_candidates))) {
    candidate <- top_localized_candidates[candidate_position]
    object <- layer_objects[[candidate$layer]]
    component_position <- match(
      candidate$product_pc, object$product_index$product_pc
    )
    score <- object$scores[, component_position]
    target <- as.numeric(object$subject_data[[candidate$phenotype]])
    orientation <- ifelse(
      suppressWarnings(cor(score, target, method = "spearman")) >= 0, 1, -1
    )
    token_values <- orientation * object$model$eig$funcs$Phi[, candidate$row_component]
    plot_midplane(
      token_values,
      sprintf("%s localized FPC%d", candidate$layer, candidate$product_pc)
    )
    theta <- orientation * object$model$eig$funcs$Theta[, candidate$time_component]
    plot(
      object$model$t_grid, theta, type = "l", lwd = 2.5, col = "#1B7837",
      xlab = "Years", ylab = "Oriented temporal loading",
      main = sprintf(
        "%s FPC%d: %s\n%s; ROI/remote %.1f",
        candidate$layer, candidate$product_pc,
        gsub("_", " ", candidate$phenotype),
        gsub("_", " ", candidate$most_localized_region),
        candidate$maximum_roi_to_remote_ratio
      )
    )
    abline(h = 0, col = "grey70", lty = 2)
  }
  dev.off()
}

session_file <- file.path(args$output_dir, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_file)
summary <- list(
  status = "complete",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  matrixfpca_refit = FALSE,
  image_reconstruction_performed = FALSE,
  observed_image_landmarks_retrieved = TRUE,
  observed_image_landmark_method = paste(
    "Actual baseline MRI from distinct subjects nearest six principal-curve",
    "percentiles; no feature inversion or decoding"
  ),
  subjects = nrow(subject_morphology),
  visits = nrow(physical_targets),
  layers = layer_specs$layer,
  top_components_tested = args$top_components,
  minimum_score_reliability = args$min_reliability,
  bootstrap_replicates = args$bootstrap_replicates,
  physical_target = paste(
    "Atlas-region mean visit-to-template log-Jacobian summarized as",
    "subject baseline level and annual slope"
  ),
  association_adjustment = paste(
    "Partial rank regression adjusted for baseline age, sex, visit count, and",
    "observed time span, omitting the tested covariate itself"
  ),
  coordinate_selection = paste(
    "Exploratory composite of absolute Spearman correlation, adjusted partial",
    "R-squared, and product-mode variance among leading FPCs with reliability at least",
    args$min_reliability
  ),
  selected_coordinates = lapply(seq_len(nrow(selected_coordinates)), function(row) {
    as.list(selected_coordinates[row])
  })
)
jsonlite::write_json(
  summary, file.path(args$output_dir, "run_summary.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null"
)
log_message("Complete. Interpretation results: ", args$output_dir)
