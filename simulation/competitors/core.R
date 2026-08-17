# Competing estimators and common evaluation utilities for simu-methods.R.
#
# All returned eigenfunctions use the array convention
#   time x row x column x component.
# The Hilbert-space inner product integrates over time (trapezoid rule) and
# sums over pixels.  Consequently, signs, component labels, and rotations
# within a leading eigenspace do not affect the subspace metrics below.

.trapz_weights <- function(t) {
  t <- as.numeric(t)
  if (length(t) < 2L || is.unsorted(t) || anyDuplicated(t)) {
    stop("t must contain at least two sorted, distinct points.")
  }
  dt <- diff(t)
  w <- numeric(length(t))
  w[1L] <- dt[1L] / 2
  w[length(t)] <- dt[length(dt)] / 2
  if (length(t) > 2L) {
    w[2L:(length(t) - 1L)] <- (dt[-length(dt)] + dt[-1L]) / 2
  }
  w
}

.as_longitudinal_arrays <- function(data) {
  if (!inherits(data, "matfdata.irreg")) {
    stop("Competitors currently require a matfdata.irreg object.")
  }
  list(
    Ly = lapply(seq_len(count_subj(data)), function(i) {
      Xi <- get_subj_tnsr(data, i)@data
      aperm(Xi, c(3L, 1L, 2L))
    }),
    Lt = get_t_list(data)
  )
}

.dct_bank <- function(t, nfun) {
  x <- (t - min(t)) / diff(range(t))
  vapply(0L:(nfun - 1L), function(k) {
    if (k == 0L) rep(1, length(x)) else sqrt(2) * cos(pi * k * x)
  }, numeric(length(x)))
}

.weighted_completion <- function(Phi, tgrid, w, n_extra) {
  if (n_extra <= 0L) return(matrix(0, nrow(Phi), 0L))
  candidates <- .dct_bank(tgrid, max(3L * n_extra, n_extra + ncol(Phi) + 2L))
  if (ncol(Phi)) {
    candidates <- candidates - Phi %*% crossprod(Phi * w, candidates)
  }
  Q <- qr.Q(qr(sqrt(w) * candidates))
  Q <- Q[, seq_len(min(n_extra, ncol(Q))), drop = FALSE] / sqrt(w)
  if (ncol(Q) < n_extra) {
    Q <- cbind(Q, matrix(0, nrow(Q), n_extra - ncol(Q)))
  }
  Q
}

.zero_pixel_fit <- function(j, k, nt, n, npc, tgrid, w, status) {
  list(
    j = j,
    k = k,
    phi = .weighted_completion(matrix(0, nt, 0L), tgrid, w, npc),
    lambda = rep(0, npc),
    sig2 = 0,
    scores = matrix(0, n, npc),
    status = status
  )
}

# Happ & Greven (2018), with a PACE decomposition for every vectorized pixel.
# This is adapted from /home/muye/Downloads/matfpca-happ.R and made independent
# of the other project from which that script sourced R/matfunc.R.
fit_happ_pace <- function(data, npc = 3L, ncomp = NULL, ncores = 1L,
                          backend = c("multisession", "multicore"),
                          center_scores = TRUE) {
  if (!requireNamespace("fdapace", quietly = TRUE)) {
    stop("Package 'fdapace' is required for happ_pace.")
  }
  backend <- match.arg(backend)
  longitudinal <- .as_longitudinal_arrays(data)
  Ly <- longitudinal$Ly
  Lt <- longitudinal$Lt
  n <- length(Ly)
  if (n < 2L) stop("happ_pace requires at least two subjects.")
  p <- dim(Ly[[1L]])[2L]
  q <- dim(Ly[[1L]])[3L]
  tgrid <- data$t_grid
  nt <- length(tgrid)
  w <- .trapz_weights(tgrid)
  npc <- max(1L, as.integer(npc))

  pixel_index <- expand.grid(j = seq_len(p), k = seq_len(q))
  pixel_task <- function(index) {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1L)
    }
    j <- pixel_index$j[index]
    k <- pixel_index$k[index]
    zero_fit <- function(status) {
      .zero_pixel_fit(j, k, nt, n, npc, tgrid, w, status)
    }
    y <- lapply(Ly, function(Y) as.numeric(Y[, j, k]))
    y_all <- unlist(y, use.names = FALSE)
    if (length(y_all) < 2L || !all(is.finite(y_all)) ||
        max(abs(y_all - mean(y_all))) <= 1e-12 + 1e-8 * max(1, abs(y_all))) {
      return(zero_fit("degenerate"))
    }

    tryCatch({
      fpca <- suppressWarnings(fdapace::FPCA(
        Ly = y,
        Lt = Lt,
        optns = list(maxK = npc, usergrid = TRUE, FVEthreshold = 0.998)
      ))
      if (is.null(fpca$phi) || is.null(fpca$lambda)) {
        return(zero_fit("empty_fpca"))
      }
      phi_raw <- as.matrix(fpca$phi)
      lambda_raw <- pmax(as.numeric(fpca$lambda), 0)
      n_available <- min(ncol(phi_raw), length(lambda_raw), npc)
      if (n_available < 1L || !any(lambda_raw[seq_len(n_available)] > 0)) {
        return(zero_fit("zero_eigen"))
      }
      phi_raw <- fdapace::ConvertSupport(
        fromGrid = fpca$workGrid,
        toGrid = tgrid,
        phi = phi_raw[, seq_len(n_available), drop = FALSE]
      )

      qr_fit <- qr(sqrt(w) * phi_raw)
      phi <- qr.Q(qr_fit) / sqrt(w)
      rotation <- qr.R(qr_fit)
      flip <- colSums(phi * phi_raw) < 0
      if (any(flip)) {
        phi[, flip] <- -phi[, flip, drop = FALSE]
        rotation[flip, ] <- -rotation[flip, , drop = FALSE]
      }

      scores <- if (!is.null(fpca$xiEst)) {
        as.matrix(fpca$xiEst)[, seq_len(n_available), drop = FALSE] %*%
          t(rotation[seq_len(n_available), seq_len(n_available), drop = FALSE])
      } else {
        curves <- fdapace::FPCApred(
          fpca, newLy = y, newLt = Lt, K = n_available
        )$ypred
        t(vapply(seq_len(n), function(i) {
          as.numeric(crossprod(phi[, seq_len(n_available), drop = FALSE], w * curves[i, ]))
        }, numeric(n_available)))
      }

      if (n_available < npc) {
        phi <- cbind(
          phi,
          .weighted_completion(phi, tgrid, w, npc - n_available)
        )
        scores <- cbind(scores, matrix(0, n, npc - n_available))
      }
      sig2 <- as.numeric(fpca$sigma2)[1L]
      if (!is.finite(sig2) || sig2 < 0) sig2 <- 0
      if (any(!is.finite(phi)) || any(!is.finite(scores))) {
        return(zero_fit("nonfinite_fpca"))
      }
      list(
        j = j,
        k = k,
        phi = phi[, seq_len(npc), drop = FALSE],
        lambda = c(lambda_raw[seq_len(n_available)], rep(0, npc - n_available)),
        sig2 = sig2,
        scores = scores[, seq_len(npc), drop = FALSE],
        status = "ok"
      )
    }, error = function(error) {
      ans <- zero_fit("fpca_error")
      ans$error <- conditionMessage(error)
      ans
    })
  }

  indices <- seq_len(nrow(pixel_index))
  pixel_fits <- if (ncores <= 1L) {
    lapply(indices, pixel_task)
  } else {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("Packages 'future' and 'future.apply' are required when ncores > 1.")
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    strategy <- if (backend == "multicore") future::multicore else future::multisession
    future::plan(strategy, workers = as.integer(ncores))
    future.apply::future_lapply(
      indices, pixel_task, future.seed = TRUE, future.packages = "fdapace"
    )
  }

  local_phi <- array(0, c(nt, p, q, npc))
  score_matrix <- matrix(0, n, p * q * npc)
  for (pixel_fit in pixel_fits) {
    j <- pixel_fit$j
    k <- pixel_fit$k
    local_phi[, j, k, ] <- pixel_fit$phi
    block <- ((j - 1L) * q + k - 1L) * npc + seq_len(npc)
    score_matrix[, block] <- pixel_fit$scores
  }
  if (center_scores) {
    score_matrix <- scale(score_matrix, center = TRUE, scale = FALSE)
  }
  score_svd <- svd(score_matrix / sqrt(n - 1L), nu = 0L)
  tolerance <- max(dim(score_matrix)) * .Machine$double.eps * max(score_svd$d, 0)
  rank_effective <- sum(score_svd$d > tolerance)
  if (rank_effective < 1L) stop("The HAPP score covariance has zero numerical rank.")
  if (is.null(ncomp)) ncomp <- rank_effective
  ncomp <- min(max(1L, as.integer(ncomp)), rank_effective)
  coefficients <- score_svd$v[, seq_len(ncomp), drop = FALSE]
  eigenvalues <- score_svd$d[seq_len(ncomp)]^2
  eigenfunctions <- array(0, c(nt, p, q, ncomp))
  for (j in seq_len(p)) {
    for (k in seq_len(q)) {
      block <- ((j - 1L) * q + k - 1L) * npc + seq_len(npc)
      eigenfunctions[, j, k, ] <- local_phi[, j, k, ] %*%
        coefficients[block, , drop = FALSE]
    }
  }
  status <- vapply(pixel_fits, `[[`, character(1L), "status")
  valid_sig2 <- vapply(pixel_fits, function(x) {
    identical(x$status, "ok") && is.finite(x$sig2) && x$sig2 >= 0
  }, logical(1L))
  sig2 <- if (any(valid_sig2)) {
    mean(vapply(pixel_fits[valid_sig2], `[[`, numeric(1L), "sig2"))
  } else {
    NA_real_
  }
  canonicalize_eigenfit(list(
    psi = eigenfunctions,
    nu = eigenvalues,
    sig2 = sig2,
    diagnostics = list(pixel_status = table(status)),
    tgrid = tgrid
  ))
}

# A deliberately simple nonseparable baseline: ridge-smooth every subject's
# pixel trajectories in a common B-spline basis, then perform PCA in the
# resulting Hilbert-space coefficient coordinates.  Unlike MatFPCA, it makes
# no product-covariance assumption.  Unlike HAPP-PACE, it does not estimate a
# separate covariance smoother for every pixel.
fit_spline_pca <- function(data, nbasis = 7L, ncomp = NULL,
                           smoothness = 1e-4, ridge = 1e-8) {
  longitudinal <- .as_longitudinal_arrays(data)
  Ly <- longitudinal$Ly
  Lt <- longitudinal$Lt
  n <- length(Ly)
  if (n < 2L) stop("spline_pca requires at least two subjects.")
  p <- dim(Ly[[1L]])[2L]
  q <- dim(Ly[[1L]])[3L]
  npix <- p * q
  nbasis <- max(4L, as.integer(nbasis))
  if (!is.finite(smoothness) || smoothness < 0 ||
      !is.finite(ridge) || ridge < 0) {
    stop("smoothness and ridge must be finite and nonnegative.")
  }
  basis <- fda::create.bspline.basis(range(data$t_grid), nbasis = nbasis, norder = 4L)
  basis_grid <- fda::eval.basis(data$t_grid, basis)
  gram <- fda::inprod(basis, basis)
  penalty <- fda::eval.penalty(basis, Lfdobj = 2L)
  gram_root <- chol(0.5 * (gram + t(gram)))

  whitened <- matrix(0, n, nbasis * npix)
  for (i in seq_len(n)) {
    B <- fda::eval.basis(Lt[[i]], basis)
    Yi <- matrix(Ly[[i]], nrow = length(Lt[[i]]), ncol = npix)
    H <- crossprod(B) + smoothness * penalty
    scale_H <- max(1, max(abs(H)))
    H <- H + ridge * scale_H * diag(nbasis)
    coef_i <- tryCatch(
      solve(H, crossprod(B, Yi)),
      error = function(error) qr.solve(H, crossprod(B, Yi), tol = 1e-10)
    )
    whitened[i, ] <- as.vector(gram_root %*% coef_i)
  }
  whitened <- scale(whitened, center = TRUE, scale = FALSE)
  coef_svd <- svd(whitened / sqrt(n - 1L), nu = 0L)
  tolerance <- max(dim(whitened)) * .Machine$double.eps * max(coef_svd$d, 0)
  rank_effective <- sum(coef_svd$d > tolerance)
  if (rank_effective < 1L) stop("The spline coefficient covariance has zero numerical rank.")
  if (is.null(ncomp)) ncomp <- rank_effective
  ncomp <- min(max(1L, as.integer(ncomp)), rank_effective)

  eigenfunctions <- array(0, c(length(data$t_grid), p, q, ncomp))
  for (component in seq_len(ncomp)) {
    white_coef <- matrix(coef_svd$v[, component], nbasis, npix)
    function_coef <- backsolve(gram_root, white_coef)
    eigenfunctions[, , , component] <- array(
      basis_grid %*% function_coef,
      c(length(data$t_grid), p, q)
    )
  }
  canonicalize_eigenfit(list(
    psi = eigenfunctions,
    nu = coef_svd$d[seq_len(ncomp)]^2,
    sig2 = NA_real_,
    diagnostics = list(nbasis = nbasis, smoothness = smoothness),
    tgrid = data$t_grid
  ))
}

canonicalize_eigenfit <- function(fit, tolerance = 1e-10) {
  stopifnot(length(dim(fit$psi)) == 4L, length(fit$nu) == dim(fit$psi)[4L])
  dims <- dim(fit$psi)
  weights <- rep(.trapz_weights(fit$tgrid), dims[2L] * dims[3L])
  P <- matrix(fit$psi, nrow = prod(dims[1L:3L]), ncol = dims[4L])
  eigenvalues <- pmax(as.numeric(fit$nu), 0)
  keep <- is.finite(eigenvalues) & eigenvalues > 0 &
    colSums(is.finite(P)) == nrow(P)
  if (!any(keep)) stop("The eigenfit contains no finite positive component.")
  P <- P[, keep, drop = FALSE]
  eigenvalues <- eigenvalues[keep]

  # Eigen-decompose diag(sqrt(w)) C_hat diag(sqrt(w)) through its small Gram
  # matrix. This preserves the estimated covariance while enforcing numerical
  # weighted orthonormality and descending eigenvalue order.
  M <- sweep(sqrt(weights) * P, 2L, sqrt(eigenvalues), "*")
  gram <- crossprod(M)
  ee <- eigen(0.5 * (gram + t(gram)), symmetric = TRUE)
  threshold <- tolerance * max(ee$values, 0)
  keep <- ee$values > threshold
  values <- ee$values[keep]
  vectors_w <- M %*% sweep(
    ee$vectors[, keep, drop = FALSE], 2L, sqrt(values), "/"
  )
  P <- vectors_w / sqrt(weights)
  fit$psi <- array(P, c(dims[1L:3L], length(values)))
  fit$nu <- values
  fit
}

product_covariance_eigenfit <- function(covariance, ncomp) {
  stopifnot(inherits(covariance, "matfcov"))
  weights <- .trapz_weights(covariance$t_grid)
  sqrt_w <- sqrt(weights)
  time_ee <- eigen(
    (sqrt_w * covariance$K) * rep(sqrt_w, each = length(sqrt_w)),
    symmetric = TRUE
  )
  row_ee <- eigen(0.5 * (covariance$U + t(covariance$U)), symmetric = TRUE)
  col_ee <- eigen(0.5 * (covariance$V + t(covariance$V)), symmetric = TRUE)
  time_values <- pmax(time_ee$values, 0)
  row_values <- pmax(row_ee$values, 0)
  col_values <- pmax(col_ee$values, 0)
  products <- as.vector(time_values %o% row_values %o% col_values)
  positive <- which(products > 10 * .Machine$double.eps * max(products))
  ordering <- positive[order(products[positive], decreasing = TRUE)]
  ordering <- ordering[seq_len(min(as.integer(ncomp), length(ordering)))]
  index <- arrayInd(ordering, c(length(time_values), length(row_values), length(col_values)))
  psi <- array(0, c(
    length(weights), nrow(covariance$U), nrow(covariance$V), length(ordering)
  ))
  for (component in seq_along(ordering)) {
    psi[, , , component] <- (time_ee$vectors[, index[component, 1L]] / sqrt_w) %o%
      row_ee$vectors[, index[component, 2L]] %o%
      col_ee$vectors[, index[component, 3L]]
  }
  canonicalize_eigenfit(list(
    psi = psi,
    nu = products[ordering],
    sig2 = covariance$sig2,
    tgrid = covariance$t_grid
  ))
}

.identifiable_evaluation_rank <- function(covariance, requested, tie_tolerance = 1e-8) {
  time_values <- pmax(eigen(
    (sqrt(.trapz_weights(covariance$t_grid)) * covariance$K) *
      rep(sqrt(.trapz_weights(covariance$t_grid)), each = nrow(covariance$K)),
    symmetric = TRUE,
    only.values = TRUE
  )$values, 0)
  row_values <- pmax(eigen(covariance$U, symmetric = TRUE, only.values = TRUE)$values, 0)
  col_values <- pmax(eigen(covariance$V, symmetric = TRUE, only.values = TRUE)$values, 0)
  spectrum <- sort(
    as.vector(time_values %o% row_values %o% col_values),
    decreasing = TRUE
  )
  spectrum <- spectrum[spectrum > 10 * .Machine$double.eps * max(spectrum)]
  rank <- min(max(1L, as.integer(requested)), length(spectrum))
  # If the requested boundary cuts an equal-eigenvalue block, include the
  # entire block. The associated individual eigenvectors/order are not
  # identifiable, whereas their joint eigenspace is.
  while (rank < length(spectrum) &&
         abs(spectrum[rank + 1L] - spectrum[rank]) <=
           tie_tolerance * max(1, spectrum[rank])) {
    rank <- rank + 1L
  }
  rank
}

.apply_product_covariance <- function(covariance, vector) {
  nt <- length(covariance$t_grid)
  p <- nrow(covariance$U)
  q <- nrow(covariance$V)
  X <- rTensor::as.tensor(array(vector, c(nt, p, q)))
  as.vector(rTensor::ttl(
    X,
    list(covariance$K, covariance$U, covariance$V),
    ms = 1L:3L
  )@data)
}

evaluate_full_eigenfit <- function(fit, truth, eval_ncomp = 5L) {
  fit <- canonicalize_eigenfit(fit)
  identifiable_rank <- .identifiable_evaluation_rank(truth, eval_ncomp)
  truth_eigenfit <- product_covariance_eigenfit(
    truth, ncomp = identifiable_rank
  )
  dims <- dim(fit$psi)
  weights <- rep(.trapz_weights(fit$tgrid), dims[2L] * dims[3L])
  P_hat <- matrix(fit$psi, nrow = length(weights))
  P_true <- matrix(truth_eigenfit$psi, nrow = length(weights))
  r <- min(identifiable_rank, ncol(P_true))
  r_hat <- min(r, ncol(P_hat))
  P_hat <- P_hat[, seq_len(r_hat), drop = FALSE]
  P_true <- P_true[, seq_len(r), drop = FALSE]
  overlap <- crossprod(sqrt(weights) * P_true, sqrt(weights) * P_hat)
  singular_values <- pmin(svd(overlap, nu = 0L, nv = 0L)$d, 1)

  lambda_hat <- c(fit$nu[seq_len(r_hat)], rep(0, r - r_hat))
  lambda_true <- truth_eigenfit$nu[seq_len(r)]
  fit_norm2 <- sum(fit$nu^2)
  cross_term <- 0
  P_all <- matrix(fit$psi, nrow = length(weights))
  for (component in seq_along(fit$nu)) {
    weighted_component <- weights * P_all[, component]
    cross_term <- cross_term + fit$nu[component] * sum(
      weighted_component * .apply_product_covariance(truth, weighted_component)
    )
  }
  truth_norm2 <- matfcovNorm2(truth)
  covariance_error2 <- max(0, fit_norm2 + truth_norm2 - 2 * cross_term)
  time_weights <- .trapz_weights(truth$t_grid)
  time_values <- pmax(eigen(
    (sqrt(time_weights) * truth$K) * rep(sqrt(time_weights), each = length(time_weights)),
    symmetric = TRUE, only.values = TRUE
  )$values, 0)
  row_values <- pmax(eigen(truth$U, symmetric = TRUE, only.values = TRUE)$values, 0)
  col_values <- pmax(eigen(truth$V, symmetric = TRUE, only.values = TRUE)$values, 0)
  truth_spectrum <- sort(
    as.vector(time_values %o% row_values %o% col_values), decreasing = TRUE
  )
  oracle_rank <- min(length(fit$nu), length(truth_spectrum))
  oracle_rank_error <- max(
    0, (truth_norm2 - sum(truth_spectrum[seq_len(oracle_rank)]^2)) / truth_norm2
  )
  relative_covariance_error <- covariance_error2 / truth_norm2

  c(
    rel_hs_cov = relative_covariance_error,
    oracle_rank_rel_hs = oracle_rank_error,
    excess_rel_hs_cov = max(0, relative_covariance_error - oracle_rank_error),
    mean_sin2 = max(0, 1 - sum(singular_values^2) / r),
    max_sin2 = if (r_hat < r) 1 else max(0, 1 - min(singular_values)^2),
    rel_eigenvalue_sse = sum((lambda_hat - lambda_true)^2) / sum(lambda_true^2),
    evaluated_components = r,
    estimated_rank = length(fit$nu),
    sig2_hat = fit$sig2
  )
}

evaluate_product_fit <- function(fit, truth, eval_ncomp = 5L,
                                 covariance_ncomp = NULL) {
  if (!is.null(covariance_ncomp)) {
    return(evaluate_full_eigenfit(
      product_covariance_eigenfit(fit$C, covariance_ncomp),
      truth,
      eval_ncomp = eval_ncomp
    ))
  }
  covariance_error2 <- max(0, matfcovDiffNorm2(fit$C, truth))
  identifiable_rank <- .identifiable_evaluation_rank(truth, eval_ncomp)
  eigen_metrics <- evaluate_full_eigenfit(
    product_covariance_eigenfit(fit$C, identifiable_rank),
    truth,
    eval_ncomp = eval_ncomp
  )
  eigen_metrics["rel_hs_cov"] <- covariance_error2 / matfcovNorm2(truth)
  eigen_metrics["oracle_rank_rel_hs"] <- 0
  eigen_metrics["excess_rel_hs_cov"] <- eigen_metrics["rel_hs_cov"]
  eigen_metrics["sig2_hat"] <- fit$C$sig2
  time_weights <- .trapz_weights(fit$C$t_grid)
  time_values <- pmax(eigen(
    (sqrt(time_weights) * fit$C$K) * rep(sqrt(time_weights), each = length(time_weights)),
    symmetric = TRUE, only.values = TRUE
  )$values, 0)
  row_values <- pmax(eigen(fit$C$U, symmetric = TRUE, only.values = TRUE)$values, 0)
  col_values <- pmax(eigen(fit$C$V, symmetric = TRUE, only.values = TRUE)$values, 0)
  product_values <- as.vector(time_values %o% row_values %o% col_values)
  eigen_metrics["estimated_rank"] <- sum(
    product_values > 10 * .Machine$double.eps * max(product_values, 0)
  )
  eigen_metrics
}

# Keep package-specific and research-prototype adapters isolated from the
# common metrics above. Source order matters: common.R uses these utilities.
.competitor_directory <- file.path("simulation", "competitors")
for (.competitor_file in c(
  "common.R", "happ-face.R", "mfaces.R", "mgsfpca.R",
  "product-fpca.R", "shi-mdfpca.R", "spaco-minus.R"
)) {
  source(file.path(.competitor_directory, .competitor_file))
}
rm(.competitor_directory, .competitor_file)
