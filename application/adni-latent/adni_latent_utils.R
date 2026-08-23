source(file.path("application", "adni-latent", "latent_io.R"))

trapezoid_weights <- function(grid) {
  grid <- as.numeric(grid)
  if (length(grid) < 2L || is.unsorted(grid, strictly = TRUE)) {
    stop("The temporal grid must be strictly increasing and contain at least two points.")
  }
  if (length(grid) == 2L) return(rep(diff(grid) / 2, 2L))
  c(
    (grid[[2L]] - grid[[1L]]) / 2,
    (grid[3:length(grid)] - grid[1:(length(grid) - 2L)]) / 2,
    (grid[[length(grid)]] - grid[[length(grid) - 1L]]) / 2
  )
}

eigen_matfcov_quadrature <- function(covariance) {
  stopifnot(is(covariance, "matfcov"))
  eig_U <- eigen(0.5 * (covariance$U + t(covariance$U)), symmetric = TRUE)
  eig_V <- eigen(0.5 * (covariance$V + t(covariance$V)), symmetric = TRUE)
  weights <- trapezoid_weights(covariance$t_grid)
  sqrt_weights <- sqrt(weights)
  weighted_K <- sweep(
    sweep(covariance$K, 1L, sqrt_weights, "*"),
    2L, sqrt_weights, "*"
  )
  eig_K <- eigen(0.5 * (weighted_K + t(weighted_K)), symmetric = TRUE)
  Theta <- sweep(eig_K$vectors, 1L, sqrt_weights, "/")
  result <- create_matrix_eigobj(
    eig_U$vectors, eig_V$vectors, Theta,
    eig_U$values, eig_V$values, eig_K$values,
    covariance$t_grid
  )
  attr(result, "quadrature_weights") <- weights
  result
}

rank_for_pve <- function(values, pve, cap = Inf) {
  values <- pmax(as.numeric(values), 0)
  if (!length(values) || !is.finite(sum(values)) || sum(values) <= 0) {
    stop("Cannot select a rank from a non-positive eigenspectrum.")
  }
  rank <- which(cumsum(values) / sum(values) >= pve)[[1L]]
  min(rank, as.integer(min(cap, length(values))))
}

positive_mode_rank <- function(values, relative_tolerance = 1e-10) {
  values <- as.numeric(values)
  if (!length(values) || !is.finite(max(values)) || max(values) <= 0) {
    stop("Cannot select positive modes from a non-positive eigenspectrum.")
  }
  max(1L, sum(values > max(values) * relative_tolerance))
}

effective_rank <- function(values) {
  values <- pmax(as.numeric(values), 0)
  if (!length(values) || sum(values) <= 0) return(NA_real_)
  sum(values)^2 / sum(values^2)
}

fit_feature_projection <- function(
  metadata, pve = 0.98, rank_cap = Inf, subject_column = "subject_id"
) {
  if (!is.numeric(pve) || length(pve) != 1L || !is.finite(pve) ||
      pve <= 0 || pve > 1) {
    stop("Feature-PCA pve must be one finite number in (0, 1].")
  }
  if (!is.numeric(rank_cap) || length(rank_cap) != 1L ||
      is.na(rank_cap) || rank_cap < 1) {
    stop("Feature-PCA rank_cap must be a positive number or Inf.")
  }
  first <- read_npy_matrix(metadata$latent_path[[1L]])
  feature_count <- ncol(first)
  token_count <- nrow(first)
  feature_sum <- numeric(feature_count)
  feature_crossproduct <- matrix(0, feature_count, feature_count)
  row_count <- 0
  subject_counts <- table(metadata[[subject_column]])
  visit_weights <- 1 / as.numeric(subject_counts[metadata[[subject_column]]])
  total_weight <- 0
  squared_weight_sum <- 0
  for (row in seq_len(nrow(metadata))) {
    latent <- read_npy_matrix(metadata$latent_path[[row]])
    if (!identical(dim(latent), c(token_count, feature_count))) {
      stop("Inconsistent latent shape at metadata row ", row)
    }
    weight <- visit_weights[[row]]
    feature_sum <- feature_sum + weight * colSums(latent)
    feature_crossproduct <- feature_crossproduct + weight * crossprod(latent)
    row_count <- row_count + nrow(latent)
    total_weight <- total_weight + weight * nrow(latent)
    squared_weight_sum <- squared_weight_sum + nrow(latent) * weight^2
  }
  center <- feature_sum / total_weight
  covariance_denominator <- total_weight - squared_weight_sum / total_weight
  covariance <- (
    feature_crossproduct - total_weight * tcrossprod(center)
  ) / max(.Machine$double.eps, covariance_denominator)
  decomposition <- eigen(0.5 * (covariance + t(covariance)), symmetric = TRUE)
  positive <- which(decomposition$values > max(decomposition$values) * 1e-10)
  if (!length(positive)) stop("The latent feature covariance has no positive rank.")
  positive_values <- pmax(decomposition$values[positive], 0)
  target_rank <- which(cumsum(positive_values) / sum(positive_values) >= pve)[[1L]]
  retained_rank <- min(target_rank, length(positive), floor(rank_cap))
  retained_pve <- sum(positive_values[seq_len(retained_rank)]) / sum(positive_values)
  if (retained_pve + 1e-12 < pve) {
    warning(
      "Feature-PCA rank cap retained only ", sprintf("%.2f%%", 100 * retained_pve),
      " variance, below the requested ", sprintf("%.2f%%", 100 * pve), ".",
      call. = FALSE
    )
  }
  list(
    center = center,
    rotation = decomposition$vectors[, seq_len(retained_rank), drop = FALSE],
    eigenvalues = decomposition$values,
    retained_rank = retained_rank,
    requested_pve = pve,
    retained_pve = retained_pve,
    token_count = token_count,
    original_feature_count = feature_count,
    training_token_rows = row_count,
    effective_weighted_token_rows = total_weight,
    training_subjects = length(subject_counts),
    subject_weighted = TRUE
  )
}

project_latent_matrix <- function(latent, projection) {
  sweep(latent, 2L, projection$center, "-") %*% projection$rotation
}

unproject_latent_matrix <- function(latent, projection) {
  sweep(latent %*% t(projection$rotation), 2L, projection$center, "+")
}

fit_latent_mean_streaming <- function(
  metadata, projection, df, boundary, subject_column = "subject_id"
) {
  basis <- splines::ns(
    metadata$time_year,
    df = as.integer(df) - 1L,
    intercept = FALSE,
    Boundary.knots = boundary
  )
  design <- cbind(intercept = 1, basis)
  counts <- table(metadata[[subject_column]])
  weights <- 1 / as.numeric(counts[metadata[[subject_column]]])
  normal <- crossprod(design, design * weights)
  output_size <- projection$token_count * projection$retained_rank
  rhs <- matrix(0, ncol(design), output_size)
  for (row in seq_len(nrow(metadata))) {
    projected <- project_latent_matrix(
      read_npy_matrix(metadata$latent_path[[row]]), projection
    )
    rhs <- rhs + tcrossprod(design[row, ] * weights[[row]], as.vector(projected))
  }
  ridge <- 1e-8 * max(mean(diag(normal)), 1)
  list(
    coefficients = solve(normal + diag(ridge, ncol(normal)), rhs),
    knots = attr(basis, "knots"),
    Boundary.knots = attr(basis, "Boundary.knots"),
    token_count = projection$token_count,
    feature_count = projection$retained_rank
  )
}

latent_mean_design <- function(model, times) {
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

predict_latent_mean <- function(model, times) {
  values <- latent_mean_design(model, times) %*% model$coefficients
  if (length(times) == 1L) {
    return(matrix(values, nrow = model$token_count, ncol = model$feature_count))
  }
  lapply(seq_along(times), function(index) {
    matrix(values[index, ], nrow = model$token_count, ncol = model$feature_count)
  })
}

evaluate_temporal_functions <- function(eig, times) {
  values <- vapply(seq_len(ncol(eig$funcs$Theta)), function(component) {
    stats::approx(
      eig$funcs$t_grid, eig$funcs$Theta[, component],
      xout = times, rule = 2, ties = "ordered"
    )$y
  }, numeric(length(times)))
  matrix(values, nrow = length(times), ncol = ncol(eig$funcs$Theta))
}

reconstruct_subject_latent <- function(model, subject_id, times, original_features = FALSE) {
  subject_position <- match(subject_id, model$subject_ids)
  if (is.na(subject_position)) stop("Unknown subject: ", subject_id)
  times <- as.numeric(times)
  theta <- evaluate_temporal_functions(model$eig, times)
  subject_scores <- model$scores[subject_position, , , , drop = FALSE]
  dim(subject_scores) <- c(
    model$ranks[["row"]], model$ranks[["column"]], model$ranks[["time"]]
  )
  means <- predict_latent_mean(model$mean_model, times)
  if (length(times) == 1L) means <- list(means)
  reconstructions <- lapply(seq_along(times), function(time_position) {
    centered <- matrix(
      0, model$representation_shape[[1L]], model$representation_shape[[2L]]
    )
    for (time_component in seq_len(model$ranks[["time"]])) {
      score_slice <- matrix(
        subject_scores[, , time_component],
        nrow = model$ranks[["row"]],
        ncol = model$ranks[["column"]]
      )
      centered <- centered + theta[time_position, time_component] *
        model$eig$funcs$Phi %*%
        score_slice %*%
        t(model$eig$funcs$Psi)
    }
    result <- centered + means[[time_position]]
    if (original_features) {
      result <- unproject_latent_matrix(result, model$feature_projection)
    }
    result
  })
  if (length(times) == 1L) reconstructions[[1L]] else reconstructions
}

ordered_product_scores <- function(scores, eig, max_components = Inf) {
  values <- eig$values$lam_phi %o% eig$values$lam_psi %o% eig$values$lam_theta
  index <- data.table::as.data.table(arrayInd(seq_along(values), dim(values)))
  data.table::setnames(index, c("row_component", "column_component", "time_component"))
  index[, eigenvalue := as.vector(values)]
  index <- index[is.finite(eigenvalue) & eigenvalue > 0][order(-eigenvalue)]
  index[, `:=`(
    product_pc = .I,
    cumulative_pve = cumsum(eigenvalue) / sum(eigenvalue)
  )]
  score_matrix <- matrix(scores@data, nrow = dim(scores)[[1L]])
  score_index <- index[seq_len(min(nrow(index), floor(max_components)))]
  offsets <- with(score_index, {
    row_component +
      (column_component - 1L) * dim(scores)[[2L]] +
      (time_component - 1L) * dim(scores)[[2L]] * dim(scores)[[3L]]
  })
  list(index = index, scores = score_matrix[, offsets, drop = FALSE])
}
