# Paper-faithful independent implementation of the alternating-regression,
# tensor-product B-spline MD-FPCA estimator in Shi et al. (2024). This is not
# copied from author software: no public longitudinal-image code was located.

.shi_tensor_setup <- function(data, nbasis_time, nbasis_row, nbasis_col,
                              center, mean_ridge) {
  longitudinal <- .as_longitudinal_arrays(data)
  Ly <- longitudinal$Ly
  Lt <- longitudinal$Lt
  p <- dim(Ly[[1L]])[2L]
  q <- dim(Ly[[1L]])[3L]
  row_grid <- if (p == 1L) 0.5 else seq(0, 1, length.out = p)
  col_grid <- if (q == 1L) 0.5 else seq(0, 1, length.out = q)
  nbasis_time <- min(length(data$t_grid), max(4L, as.integer(nbasis_time)))
  nbasis_row <- min(p, max(if (p < 4L) p else 4L, as.integer(nbasis_row)))
  nbasis_col <- min(q, max(if (q < 4L) q else 4L, as.integer(nbasis_col)))
  make_basis <- function(rangeval, nbasis, ngrid) {
    if (ngrid == 1L) return(NULL)
    fda::create.bspline.basis(rangeval, nbasis = nbasis,
                              norder = min(4L, nbasis))
  }
  basis_time <- make_basis(range(data$t_grid), nbasis_time, length(data$t_grid))
  basis_row <- make_basis(c(0, 1), nbasis_row, p)
  basis_col <- make_basis(c(0, 1), nbasis_col, q)
  eval_basis <- function(x, basis) {
    if (is.null(basis)) matrix(1, length(x), 1L) else fda::eval.basis(x, basis)
  }
  basis_gram <- function(basis) {
    if (is.null(basis)) matrix(1, 1L, 1L) else fda::inprod(basis, basis)
  }
  basis_penalty <- function(basis) {
    if (is.null(basis) || basis$nbasis < 3L) return(matrix(0, nrow(basis_gram(basis)), nrow(basis_gram(basis))))
    fda::eval.penalty(basis, Lfdobj = 2L)
  }

  Br <- eval_basis(row_grid, basis_row)
  Bc <- eval_basis(col_grid, basis_col)
  Bs <- kronecker(Bc, Br)
  Gr <- basis_gram(basis_row)
  Gc <- basis_gram(basis_col)
  Gt <- basis_gram(basis_time)
  Pr <- basis_penalty(basis_row)
  Pc <- basis_penalty(basis_col)
  Pt <- basis_penalty(basis_time)
  Gs <- kronecker(Gc, Gr)
  gram <- kronecker(Gt, Gs)
  penalty <- kronecker(Pt, Gs) + kronecker(
    Gt, kronecker(Pc, Gr) + kronecker(Gc, Pr)
  )
  design <- lapply(seq_along(Ly), function(i) {
    Bt <- eval_basis(Lt[[i]], basis_time)
    kronecker(Bt, Bs)
  })
  response <- lapply(Ly, function(Y) {
    as.vector(t(matrix(Y, nrow = dim(Y)[1L], ncol = p * q)))
  })
  subject_weight <- 1 / lengths(response)

  mean_coef <- rep(0, ncol(design[[1L]]))
  if (center) {
    H <- Reduce(`+`, Map(function(X, weight) weight * crossprod(X),
                         design, subject_weight))
    b <- Reduce(`+`, Map(function(X, y, weight) weight * crossprod(X, y),
                         design, response, subject_weight))
    scale_H <- max(1, max(abs(diag(H))))
    mean_coef <- tryCatch(
      solve(H + mean_ridge * scale_H * diag(nrow(H)), b),
      error = function(error) qr.solve(
        H + mean_ridge * scale_H * diag(nrow(H)), b, tol = 1e-10
      )
    )
  }
  centered <- Map(function(X, y) as.numeric(y - X %*% mean_coef), design, response)
  list(
    design = design, response = centered, raw_response = response,
    subject_weight = subject_weight, gram = gram, penalty = penalty,
    mean_coef = mean_coef, Bs = Bs, basis_time = basis_time,
    basis_row = basis_row, basis_col = basis_col,
    eval_basis = eval_basis, p = p, q = q,
    nbasis = c(time = ncol(Gt), row = ncol(Gr), col = ncol(Gc))
  )
}

.shi_orthonormalize <- function(beta, gram, previous = NULL, tolerance = 1e-10) {
  if (!is.null(previous) && ncol(previous)) {
    beta <- beta - previous %*% as.numeric(crossprod(previous, gram %*% beta))
  }
  norm_beta <- sqrt(max(0, as.numeric(crossprod(beta, gram %*% beta))))
  if (!is.finite(norm_beta) || norm_beta <= tolerance) return(NULL)
  as.numeric(beta / norm_beta)
}

.shi_solve_constrained <- function(H, b, gram, previous, ridge) {
  scale_H <- max(1, max(abs(diag(H))))
  H <- 0.5 * (H + t(H)) + ridge * scale_H * diag(nrow(H))
  if (is.null(previous) || !ncol(previous)) {
    return(tryCatch(solve(H, b), error = function(error) qr.solve(H, b, tol = 1e-10)))
  }
  constraints <- crossprod(previous, gram)
  KKT <- rbind(
    cbind(H, t(constraints)),
    cbind(constraints, matrix(0, nrow(constraints), nrow(constraints)))
  )
  rhs <- c(b, rep(0, nrow(constraints)))
  solution <- tryCatch(
    solve(KKT, rhs),
    error = function(error) qr.solve(KKT, rhs, tol = 1e-10)
  )
  solution[seq_len(nrow(H))]
}

fit_shi_mdfpca <- function(data, ncomp = 5L, nbasis_time = 7L,
                           nbasis_row = 5L, nbasis_col = 5L,
                           center = TRUE, smoothness = 0,
                           maxit = 100L, tolerance = 1e-5,
                           ridge = 1e-8, mean_ridge = 1e-6,
                           noise_correct = TRUE, trace = FALSE) {
  setup <- .shi_tensor_setup(
    data, nbasis_time, nbasis_row, nbasis_col, center, mean_ridge
  )
  X <- setup$design
  y <- setup$response
  n <- length(y)
  if (n < 2L) stop("shi_mdfpca requires at least two subjects.")
  ncomp <- min(max(1L, as.integer(ncomp)), n - 1L, ncol(X[[1L]]))
  maxit <- max(1L, as.integer(maxit))
  if (!is.finite(smoothness) || smoothness < 0) stop("smoothness must be nonnegative.")

  # A stable, data-driven initialization. It affects only the starting point
  # of the published alternating least-squares criterion.
  initial_coef <- t(vapply(seq_len(n), function(i) {
    H_i <- crossprod(X[[i]])
    scale_i <- max(1, max(abs(diag(H_i))))
    as.numeric(tryCatch(
      solve(H_i + 1e-4 * scale_i * diag(nrow(H_i)), crossprod(X[[i]], y[[i]])),
      error = function(error) qr.solve(
        H_i + 1e-4 * scale_i * diag(nrow(H_i)), crossprod(X[[i]], y[[i]]),
        tol = 1e-9
      )
    ))
  }, numeric(ncol(X[[1L]]))))
  initial_coef <- scale(initial_coef, center = TRUE, scale = FALSE)
  gram_root <- chol(0.5 * (setup$gram + t(setup$gram)))
  initial_svd <- svd(initial_coef %*% t(gram_root), nu = 0L)

  coefficients <- matrix(0, nrow(setup$gram), 0L)
  iterations <- integer(ncomp)
  converged <- logical(ncomp)
  for (component in seq_len(ncomp)) {
    if (component <= ncol(initial_svd$v)) {
      beta <- backsolve(gram_root, initial_svd$v[, component])
    } else {
      beta <- cos(seq_len(nrow(setup$gram)) * component)
    }
    beta <- .shi_orthonormalize(beta, setup$gram, coefficients)
    if (is.null(beta)) stop("Could not initialize Shi MD-FPCA component ", component, ".")

    for (iteration in seq_len(maxit)) {
      current <- cbind(coefficients, beta)
      scores <- lapply(seq_len(n), function(i) {
        stats::coef(stats::lm.fit(X[[i]] %*% current, y[[i]]))
      })
      alpha <- vapply(scores, function(score) score[component], numeric(1L))
      residual <- lapply(seq_len(n), function(i) {
        if (component == 1L) return(y[[i]])
        as.numeric(y[[i]] - X[[i]] %*% coefficients %*% scores[[i]][seq_len(component - 1L)])
      })
      H <- Reduce(`+`, lapply(seq_len(n), function(i) {
        setup$subject_weight[i] * alpha[i]^2 * crossprod(X[[i]])
      })) + smoothness * setup$penalty
      b <- Reduce(`+`, lapply(seq_len(n), function(i) {
        setup$subject_weight[i] * alpha[i] * crossprod(X[[i]], residual[[i]])
      }))
      beta_new <- .shi_solve_constrained(
        H, b, setup$gram, coefficients, ridge
      )
      beta_new <- .shi_orthonormalize(beta_new, setup$gram, coefficients)
      if (is.null(beta_new)) stop("Shi MD-FPCA component update became degenerate.")
      if (as.numeric(crossprod(beta, setup$gram %*% beta_new)) < 0) beta_new <- -beta_new
      change <- sqrt(as.numeric(crossprod(
        beta_new - beta, setup$gram %*% (beta_new - beta)
      )))
      beta <- beta_new
      if (trace) message("Shi component ", component, ", iteration ", iteration,
                         ", change = ", signif(change, 4L))
      if (is.finite(change) && change <= tolerance) {
        converged[component] <- TRUE
        break
      }
    }
    iterations[component] <- iteration
    coefficients <- cbind(coefficients, beta)
  }

  score_matrix <- t(vapply(seq_len(n), function(i) {
    as.numeric(stats::coef(stats::lm.fit(X[[i]] %*% coefficients, y[[i]])))
  }, numeric(ncomp)))
  fitted <- lapply(seq_len(n), function(i) as.numeric(
    X[[i]] %*% coefficients %*% score_matrix[i, ]
  ))
  residual_sse <- sum(vapply(seq_len(n), function(i) {
    sum((y[[i]] - fitted[[i]])^2)
  }, numeric(1L)))
  residual_df <- max(1, sum(lengths(y)) - n * ncomp)
  sig2 <- residual_sse / residual_df
  score_covariance <- stats::cov(score_matrix)
  noise_covariance <- matrix(0, ncomp, ncomp)
  if (noise_correct && is.finite(sig2) && sig2 > 0) {
    noise_covariance <- Reduce(`+`, lapply(seq_len(n), function(i) {
      score_design <- X[[i]] %*% coefficients
      gram_i <- crossprod(score_design)
      scale_i <- max(1, max(abs(diag(gram_i))))
      tryCatch(
        solve(gram_i + 1e-8 * scale_i * diag(ncomp)),
        error = function(error) matrix(0, ncomp, ncomp)
      )
    })) * sig2 / n
    score_covariance <- score_covariance - noise_covariance
  }
  score_ee <- eigen(0.5 * (score_covariance + t(score_covariance)), symmetric = TRUE)
  positive <- score_ee$values > 1e-10 * max(score_ee$values, 0)
  if (!any(positive)) stop("Shi MD-FPCA produced no positive score variance.")
  coefficients <- coefficients %*% score_ee$vectors[, positive, drop = FALSE]

  Bt <- setup$eval_basis(data$t_grid, setup$basis_time)
  grid_design <- kronecker(Bt, setup$Bs)
  values <- grid_design %*% coefficients
  psi <- array(0, c(length(data$t_grid), setup$p, setup$q, ncol(values)))
  for (component in seq_len(ncol(values))) {
    psi[, , , component] <- aperm(
      array(values[, component], c(setup$p, setup$q, length(data$t_grid))),
      c(3L, 1L, 2L)
    )
  }
  .truncate_eigenfit(list(
    psi = psi,
    nu = score_ee$values[positive],
    sig2 = sig2,
    tgrid = data$t_grid,
    diagnostics = list(
      implementation = "independent paper-faithful reimplementation",
      nbasis = setup$nbasis,
      iterations = iterations,
      converged = converged,
      smoothness = smoothness,
      noise_correct = noise_correct,
      noise_score_covariance = noise_covariance,
      coefficients = coefficients,
      mean_coefficients = setup$mean_coef
    )
  ), ncomp)
}
