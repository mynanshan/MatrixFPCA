# Independent SPACO- style comparator: smooth CP/PARAFAC for a tensor with
# missing longitudinal slabs and no auxiliary covariates. The authors' public
# experiment repository imports a non-shipped `spaco` module, and the PyPI
# artifact contains metadata only, so this is explicitly not claimed as their
# exact software implementation.

.spaco_tensor <- function(data) {
  longitudinal <- .as_longitudinal_arrays(data)
  n <- length(longitudinal$Ly)
  p <- dim(longitudinal$Ly[[1L]])[2L]
  q <- dim(longitudinal$Ly[[1L]])[3L]
  nt <- length(data$t_grid)
  X <- array(NA_real_, c(n, nt, p * q))
  observed <- matrix(FALSE, n, nt)
  for (i in seq_len(n)) {
    visit_index <- vapply(longitudinal$Lt[[i]], function(time) {
      which.min(abs(data$t_grid - time))
    }, integer(1L))
    discrepancy <- abs(data$t_grid[visit_index] - longitudinal$Lt[[i]])
    if (any(discrepancy > 1e-8 * max(1, diff(range(data$t_grid))))) {
      stop("spaco_minus requires observation times on data$t_grid.")
    }
    X[i, visit_index, ] <- matrix(
      longitudinal$Ly[[i]], nrow = length(visit_index), ncol = p * q
    )
    observed[i, visit_index] <- TRUE
  }
  list(X = X, observed = observed, p = p, q = q)
}

.spaco_solve <- function(design, response, penalty = NULL, ridge = 1e-8) {
  H <- crossprod(design)
  if (!is.null(penalty)) H <- H + penalty
  scale_H <- max(1, max(abs(diag(H))))
  H <- H + ridge * scale_H * diag(nrow(H))
  tryCatch(
    solve(H, crossprod(design, response)),
    error = function(error) qr.solve(H, crossprod(design, response), tol = 1e-9)
  )
}

.spaco_score_fit <- function(design, response, ridge_multiplier = 0,
                             numerical_ridge = 1e-8) {
  H <- crossprod(design)
  scale_H <- mean(diag(H))
  if (!is.finite(scale_H) || scale_H <= 0) scale_H <- 1
  numerical_scale <- max(1, max(abs(diag(H))))
  regularized <- 0.5 * (H + t(H)) +
    (ridge_multiplier * scale_H + numerical_ridge * numerical_scale) *
      diag(nrow(H))
  coefficients <- tryCatch(
    solve(regularized, crossprod(design, response)),
    error = function(error) qr.solve(
      regularized, crossprod(design, response), tol = 1e-9
    )
  )
  list(
    coefficients = as.numeric(coefficients), gram = H,
    regularized_gram = regularized,
    penalty = ridge_multiplier * scale_H
  )
}

.spaco_select_score_ridge <- function(X, observed, Phi, V,
                                      candidates = c(0, 10^seq(-4, 0, by = 0.25)),
                                      numerical_ridge = 1e-8) {
  candidates <- sort(unique(as.numeric(candidates)))
  if (!length(candidates) || any(!is.finite(candidates)) || any(candidates < 0)) {
    stop("score-ridge candidates must be finite and nonnegative.")
  }

  # Phi and V have already been estimated by the expensive ALS fit. Tuning
  # therefore requires only small rank-by-rank score solves. Sufficient
  # statistics for each grid time avoid repeatedly processing all pixels.
  time_design <- lapply(seq_len(nrow(Phi)), function(time) {
    sweep(V, 2L, Phi[time, ], "*")
  })
  time_gram <- lapply(time_design, crossprod)
  prediction_sse <- matrix(0, dim(X)[1L], length(candidates))
  prediction_n <- matrix(0, dim(X)[1L], length(candidates))

  for (i in seq_len(dim(X)[1L])) {
    times <- which(observed[i, ])
    if (length(times) < 2L) next
    response <- lapply(times, function(time) X[i, time, ])
    rhs <- Map(crossprod, time_design[times], response)
    total_gram <- Reduce(`+`, time_gram[times])
    total_rhs <- Reduce(`+`, rhs)
    for (held_out in seq_along(times)) {
      train_gram <- total_gram - time_gram[[times[held_out]]]
      train_rhs <- total_rhs - rhs[[held_out]]
      scale_H <- mean(diag(train_gram))
      if (!is.finite(scale_H) || scale_H <= 0) scale_H <- 1
      numerical_scale <- max(1, max(abs(diag(train_gram))))
      for (candidate_index in seq_along(candidates)) {
        regularized <- 0.5 * (train_gram + t(train_gram)) +
          (candidates[candidate_index] * scale_H +
             numerical_ridge * numerical_scale) * diag(nrow(train_gram))
        score <- tryCatch(
          solve(regularized, train_rhs),
          error = function(error) qr.solve(regularized, train_rhs, tol = 1e-9)
        )
        residual <- response[[held_out]] -
          time_design[[times[held_out]]] %*% score
        prediction_sse[i, candidate_index] <-
          prediction_sse[i, candidate_index] + sum(residual^2)
        prediction_n[i, candidate_index] <-
          prediction_n[i, candidate_index] + length(residual)
      }
    }
  }
  subject_error <- prediction_sse / prediction_n
  cv_error <- colMeans(subject_error, na.rm = TRUE)
  cv_se <- vapply(seq_along(candidates), function(candidate_index) {
    values <- subject_error[, candidate_index]
    values <- values[is.finite(values)]
    if (length(values) < 2L) 0 else stats::sd(values) / sqrt(length(values))
  }, numeric(1L))
  valid <- is.finite(cv_error)
  minimum <- if (any(valid)) which.min(ifelse(valid, cv_error, Inf)) else NA_integer_
  selected <- if (is.finite(minimum)) {
    # The conventional one-standard-error rule avoids selecting an unstable,
    # nearly unregularized score solve for a negligible prediction-error gain.
    eligible <- valid & cv_error <= cv_error[minimum] + cv_se[minimum]
    max(which(eligible))
  } else {
    which.min(candidates)
  }
  list(
    ridge_multiplier = candidates[selected], candidates = candidates,
    cv_error = cv_error, cv_se = cv_se, minimum_cv_ridge = if (
      is.finite(minimum)
    ) candidates[minimum] else candidates[selected]
  )
}

fit_spaco_minus <- function(data, ncomp = 5L, nbasis = 7L,
                            smoothness = 1e-2, center = TRUE,
                            maxit = 100L, tolerance = 1e-5,
                            ridge = 1e-8, score_ridge = "auto",
                            score_ridge_candidates = c(
                              0, 10^seq(-4, 0, by = 0.25)
                            ),
                            trace = FALSE) {
  tensor <- .spaco_tensor(data)
  X <- tensor$X
  observed <- tensor$observed
  n <- dim(X)[1L]
  nt <- dim(X)[2L]
  npix <- dim(X)[3L]
  rank <- min(max(1L, as.integer(ncomp)), n - 1L, nt, npix)
  automatic_score_ridge <- is.character(score_ridge) &&
    length(score_ridge) == 1L && identical(tolower(score_ridge), "auto")
  if (!automatic_score_ridge) {
    score_ridge <- as.numeric(score_ridge)[1L]
    if (!is.finite(score_ridge) || score_ridge < 0) {
      stop("score_ridge must be 'auto' or a finite nonnegative number.")
    }
  }
  nbasis <- min(nt, max(4L, as.integer(nbasis)))
  basis <- fda::create.bspline.basis(
    range(data$t_grid), nbasis = nbasis, norder = min(4L, nbasis)
  )
  B <- fda::eval.basis(data$t_grid, basis)
  penalty_basis <- fda::eval.penalty(basis, Lfdobj = 2L)
  time_weights <- .trapz_weights(data$t_grid)

  mean_surface <- matrix(0, nt, npix)
  if (center) {
    for (pixel in seq_len(npix)) {
      row <- which(observed, arr.ind = TRUE)
      response <- X[cbind(row[, 1L], row[, 2L], pixel)]
      mean_coef <- .spaco_solve(
        B[row[, 2L], , drop = FALSE], response,
        penalty = 1e-4 * penalty_basis, ridge = ridge
      )
      mean_surface[, pixel] <- B %*% mean_coef
    }
    for (i in seq_len(n)) for (time in which(observed[i, ])) {
      X[i, time, ] <- X[i, time, ] - mean_surface[time, ]
    }
  }

  # Missing slabs are zero after centering. The SVD supplies only an ALS start.
  unfolded <- matrix(0, n, nt * npix)
  for (i in seq_len(n)) for (time in which(observed[i, ])) {
    block <- (time - 1L) * npix + seq_len(npix)
    unfolded[i, block] <- X[i, time, ]
  }
  unfolded <- scale(unfolded, center = TRUE, scale = FALSE)
  initial <- svd(unfolded, nu = rank, nv = rank)
  U <- initial$u[, seq_len(rank), drop = FALSE] %*%
    diag(initial$d[seq_len(rank)], rank)
  Phi <- matrix(0, nt, rank)
  V <- matrix(0, npix, rank)
  for (component in seq_len(rank)) {
    surface <- matrix(initial$v[, component], nrow = npix, ncol = nt)
    surface_svd <- svd(t(surface), nu = 1L, nv = 1L)
    Phi[, component] <- surface_svd$u[, 1L] * sqrt(surface_svd$d[1L])
    V[, component] <- surface_svd$v[, 1L] * sqrt(surface_svd$d[1L])
  }

  normalize_factors <- function(U, Phi, V) {
    for (component in seq_len(rank)) {
      norm_time <- sqrt(sum(time_weights * Phi[, component]^2))
      norm_pixel <- sqrt(sum(V[, component]^2))
      if (!is.finite(norm_time) || norm_time <= 1e-12) norm_time <- 1
      if (!is.finite(norm_pixel) || norm_pixel <= 1e-12) norm_pixel <- 1
      Phi[, component] <- Phi[, component] / norm_time
      V[, component] <- V[, component] / norm_pixel
      U[, component] <- U[, component] * norm_time * norm_pixel
    }
    list(U = U, Phi = Phi, V = V)
  }
  normalized <- normalize_factors(U, Phi, V)
  U <- normalized$U; Phi <- normalized$Phi; V <- normalized$V
  loss_previous <- Inf
  converged <- FALSE
  for (iteration in seq_len(max(1L, as.integer(maxit)))) {
    # Conditional subject scores.
    for (i in seq_len(n)) {
      times <- which(observed[i, ])
      design_i <- do.call(rbind, lapply(times, function(time) {
        sweep(V, 2L, Phi[time, ], "*")
      }))
      response_i <- matrix(
        X[i, times, , drop = FALSE], nrow = length(times), ncol = npix
      )
      U[i, ] <- .spaco_solve(
        design_i, as.vector(t(response_i)), ridge = ridge
      )
    }

    # Pixel loadings.
    observation_index <- which(observed, arr.ind = TRUE)
    common_design <- U[observation_index[, 1L], , drop = FALSE] *
      Phi[observation_index[, 2L], , drop = FALSE]
    for (pixel in seq_len(npix)) {
      response <- X[cbind(observation_index[, 1L], observation_index[, 2L], pixel)]
      V[pixel, ] <- .spaco_solve(common_design, response, ridge = ridge)
    }

    # Joint smooth update of all temporal factors in a B-spline basis.
    nobs <- nrow(observation_index) * npix
    temporal_design <- matrix(0, nobs, nbasis * rank)
    temporal_response <- numeric(nobs)
    position <- 0L
    for (row in seq_len(nrow(observation_index))) {
      i <- observation_index[row, 1L]
      time <- observation_index[row, 2L]
      rows <- position + seq_len(npix)
      for (component in seq_len(rank)) {
        cols <- (component - 1L) * nbasis + seq_len(nbasis)
        temporal_design[rows, cols] <- (U[i, component] * V[, component]) %o% B[time, ]
      }
      temporal_response[rows] <- X[i, time, ]
      position <- position + npix
    }
    temporal_penalty <- smoothness * kronecker(diag(rank), penalty_basis)
    temporal_coef <- .spaco_solve(
      temporal_design, temporal_response, temporal_penalty, ridge
    )
    Phi <- B %*% matrix(temporal_coef, nbasis, rank)
    normalized <- normalize_factors(U, Phi, V)
    U <- normalized$U; Phi <- normalized$Phi; V <- normalized$V

    residual_sse <- 0
    total_observations <- 0L
    for (i in seq_len(n)) for (time in which(observed[i, ])) {
      prediction <- as.numeric((U[i, ] * Phi[time, ]) %*% t(V))
      residual_sse <- residual_sse + sum((X[i, time, ] - prediction)^2)
      total_observations <- total_observations + npix
    }
    loss <- residual_sse / total_observations
    relative_change <- abs(loss_previous - loss) / max(1, abs(loss_previous), abs(loss))
    if (trace) message("SPACO- iteration ", iteration, ", loss = ", signif(loss, 5L))
    if (is.finite(relative_change) && relative_change <= tolerance) {
      converged <- TRUE
      break
    }
    loss_previous <- loss
  }

  score_ridge_selection <- if (automatic_score_ridge) {
    .spaco_select_score_ridge(
      X, observed, Phi, V, candidates = score_ridge_candidates,
      numerical_ridge = ridge
    )
  } else {
    list(
      ridge_multiplier = score_ridge, candidates = score_ridge,
      cv_error = NA_real_, cv_se = NA_real_, minimum_cv_ridge = score_ridge
    )
  }
  selected_score_ridge <- score_ridge_selection$ridge_multiplier
  score_fits <- lapply(seq_len(n), function(i) {
    times <- which(observed[i, ])
    design_i <- do.call(rbind, lapply(times, function(time) {
      sweep(V, 2L, Phi[time, ], "*")
    }))
    response_i <- matrix(
      X[i, times, , drop = FALSE], nrow = length(times), ncol = npix
    )
    .spaco_score_fit(
      design_i, as.vector(t(response_i)), selected_score_ridge, ridge
    )
  })
  U <- t(vapply(score_fits, `[[`, numeric(rank), "coefficients"))

  score_covariance <- stats::cov(U)
  score_ee <- eigen(0.5 * (score_covariance + t(score_covariance)), symmetric = TRUE)
  positive <- score_ee$values > 1e-10 * max(score_ee$values, 0)
  base_functions <- array(0, c(nt, tensor$p, tensor$q, rank))
  for (component in seq_len(rank)) {
    base_functions[, , , component] <- array(
      Phi[, component] %o% V[, component],
      c(nt, tensor$p, tensor$q)
    )
  }
  base_matrix <- matrix(base_functions, nt * npix, rank)
  psi <- array(
    base_matrix %*% score_ee$vectors[, positive, drop = FALSE],
    c(nt, tensor$p, tensor$q, sum(positive))
  )
  pixel_noise <- numeric(npix)
  for (pixel in seq_len(npix)) {
    residual <- numeric()
    for (i in seq_len(n)) for (time in which(observed[i, ])) {
      prediction <- sum(U[i, ] * Phi[time, ] * V[pixel, ])
      residual <- c(residual, X[i, time, pixel] - prediction)
    }
    pixel_noise[pixel] <- mean(residual^2)
  }
  .truncate_eigenfit(list(
    psi = psi, nu = score_ee$values[positive], sig2 = mean(pixel_noise),
    tgrid = data$t_grid,
    diagnostics = list(
      implementation = "independent SPACO-minus style smooth CP",
      converged = converged, iterations = iteration, rank = rank,
      nbasis = nbasis, smoothness = smoothness, numerical_ridge = ridge,
      selected_score_ridge = selected_score_ridge,
      minimum_cv_score_ridge = score_ridge_selection$minimum_cv_ridge,
      score_ridge_candidates = score_ridge_selection$candidates,
      score_ridge_cv_error = score_ridge_selection$cv_error,
      score_ridge_cv_se = score_ridge_selection$cv_se,
      subject_scores = U, temporal_factors = Phi, pixel_loadings = V,
      pixel_noise = pixel_noise, mean_surface = mean_surface
    )
  ), ncomp)
}
