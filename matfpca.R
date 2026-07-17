library(foreach)
library(rTensor)

source("matfdata.R")

# TODO: add support for mean functions


# Helpers ----------------------------------------------------------------

offdiag <- function(X, lower.half = TRUE) {
  stopifnot(is.matrix(X))
  stopifnot(nrow(X) == ncol(X))
  n <- nrow(X)
  idx <- lower.tri(X)
  if (!lower.half) idx <- idx | upper.tri(X)
  X[idx]
}

sym <- function(X) {
  0.5 * (X + t(X))
}

tr <- function(X, Y = NULL) {
  if (is.null(Y)) return(sum(diag(X)))
  sum(t(X) * Y)
}

psd_project <- function(A) {
  A <- sym(A)
  ee <- eigen(A, symmetric = TRUE)
  ee$vectors %*% (pmax(ee$values, 0) * t(ee$vectors))
}

normalize_psd_trace <- function(A, target_trace) {
  A <- psd_project(A)
  A_trace <- tr(A)
  if (!is.finite(A_trace) || A_trace <= 10 * .Machine$double.eps) {
    stop("The de-noised covariance estimate has no positive trace.")
  }
  A * (target_trace / A_trace)
}

trapz_average <- function(x, y) {
  stopifnot(length(x) == length(y), length(x) >= 2L)
  if (is.unsorted(x, strictly = FALSE) || anyDuplicated(x)) {
    stop("x must be sorted and contain unique values.")
  }
  sum(diff(x) * (y[-1L] + y[-length(y)]) / 2) /
    (x[length(x)] - x[1L])
}

inv_spd <- function(A, ridge = 1e-8) {
  # question: what's the best choice of tol and ridge?
  A <- 0.5 * (A + t(A))
  if (ridge > 0) {
    A <- A + ridge * diag(1, nrow(A))
  }
  chol2inv(chol(A))
}

ginv_spd <- function(A, tol = 1e-3, ridge = 1e-8) {
  A <- 0.5 * (A + t(A))
  if (ridge > 0) {
    A <- A + ridge * diag(1, nrow(A))
  }
  ee <- eigen(A, symmetric = TRUE)
  val <- pmax(ee$values, 0)
  inv <- ifelse(val > tol * max(val), 1 / val, 0)
  ee$vectors %*% (inv * t(ee$vectors))
}

rescale_trace <- function(W, target_tr) {
  s <- sum(diag(W))
  if (!is.finite(s) || abs(s) < 1e-14) return(W)
  W * (target_tr / s)
}

trace_normalize <- function(A, target_trace) {
  A_trace <- tr(A)
  if (!is.finite(A_trace) || abs(A_trace) < 10 * .Machine$double.eps) {
    return(diag(target_trace / nrow(A), nrow(A)))
  }
  A * (target_trace / A_trace)
}

# Cov Estimation ----------------------------------------------------------

###### Formulate row covariances ######

get_working_K <- function(X, W1, W2) {
  # unfold_3(X) (W_2 \otimes W_1) unfold_3(X)^T
  # Input:
  #   X: one observation tensor, p x q x Ni
  #   W1, W2: row/col weight matrices, p x p and q x q
  # Return:
  #   working K, Ni x Ni matrix
  XW <- rTensor::ttl(X, list(W1, W2), ms = c(1, 2))
  rTensor::unfold(XW, row_idx = 3, col_idx = c(1, 2))@data %*%
    t(rTensor::unfold(X, row_idx = 3, col_idx = c(1, 2))@data)
}

get_working_V <- function(X, W1, W3) {
  # unfold_2(X) (W_3 \otimes W_1) unfold_2(X)^T
  # Input:
  #   X: one observation tensor, p x q x Ni
  #   W1, W3: row/time weight matrices, p x p and Ni x Ni
  # Return:
  #   working V, q x q matrix
  XW <- rTensor::ttl(X, list(W1, W3), ms = c(1, 3))
  rTensor::unfold(XW, row_idx = 2, col_idx = c(1, 3))@data %*%
    t(rTensor::unfold(X, row_idx = 2, col_idx = c(1, 3))@data)
}

get_working_U <- function(X, W2, W3) {
  # unfold_1(X) (W_3 \otimes W_2) unfold_1(X)^T
  # Input:
  #   X: one observation tensor, p x q x Ni
  #   W2, W3: col/time weight matrices, q x q and Ni x Ni
  # Return:
  #   working U, p x p matrix
  XW <- rTensor::ttl(X, list(W2, W3), ms = c(2, 3))
  rTensor::unfold(XW, row_idx = 1, col_idx = c(2, 3))@data %*%
    t(rTensor::unfold(X, row_idx = 1, col_idx = c(2, 3))@data)
}

###### Weight matrices ######

WEIGHT_SCHEMES <- c(
  "identity", "invcov", "invridcov", "gaussian_refined"
)

match_weight_scheme <- function(wScheme, choices = WEIGHT_SCHEMES) {
  if (length(wScheme) == 1L && identical(wScheme, "optimal")) {
    warning(
      "wScheme = 'optimal' is deprecated; use 'gaussian_refined'.",
      call. = FALSE
    )
    wScheme <- "gaussian_refined"
  }
  match.arg(wScheme, choices)
}

compute_W <- function(
  U, V, LK,
  sig2 = 0,
  wScheme = WEIGHT_SCHEMES,
  nIter.opt = 1,
  rescale.opt = TRUE,
  ridge.opt = 1e-8
) {
  wScheme <- match_weight_scheme(wScheme)
  eps <- 10 * .Machine$double.eps

  if (sig2 <= eps && wScheme %in% c("gaussian_refined", "invridcov")) {
    wScheme <- "invcov"
  }

  if (wScheme == "gaussian_refined") {
    W0 <- compute_W.shared(U, V, LK, sig2, "invridcov", eps)
    return(compute_W.gaussian_refined(
      U, V, LK, sig2,
      W0$W1, W0$W2, W0$W3,
      nIter = nIter.opt,
      ridge = ridge.opt,
      rescale = rescale.opt
    ))
  }

  compute_W.shared(U, V, LK, sig2, wScheme, eps)
}

compute_W.shared <- function(U, V, LK, sig2, wScheme, eps = 1e-10) {
  p <- nrow(U)
  q <- nrow(V)

  W <- rep(list(NULL), 3)
  names(W) <- paste0("W", 1:3)
  class(W) <- c("mode_weights", "shared")

  if (wScheme == "identity") {
    W$W1 <- diag(1, p)
    W$W2 <- diag(1, q)
    W$W3 <- lapply(LK, \(Ki) diag(1, nrow(Ki)))
    return(W)
  }

  if (wScheme == "invcov") {
    W$W1 <- ginv_spd(U)
    W$W2 <- ginv_spd(V)
    W$W3 <- lapply(LK, \(Ki) ginv_spd(Ki))
    return(W)
  }

  if (wScheme == "invridcov") {
    W$W1 <- inv_spd(U + sig2 * diag(1, p))
    W$W2 <- inv_spd(V + sig2 * diag(1, q))
    W$W3 <- lapply(LK, \(Ki) inv_spd(Ki + sig2 * diag(1, nrow(Ki))))
    return(W)
  }

  stop("Unknown wScheme.")
}

# scalar shorthands
.beta_mode <- function(S, W) {  # tr(S W)
  sum(diag(S %*% W))
}

.gamma_mode <- function(S, W) {  # tr(S W S W)
  sum(diag(S %*% W %*% S %*% W))
}
.zeta_mode <- function(S, W) {   # tr(S W^2)
  sum(diag(S %*% (W %*% W)))
}
.eta_mode <- function(W) {       # tr(W^2)
  sum(W * W)
}

.gaussian_block_weight <- function(
  S,
  A,
  B,
  C,
  ridge = 1e-8,
  tol = 1e-10
) {
  # Singular-safe form of S (A S^2 + B S + C I)^(-1). The inverse-covariance
  # rewrite is equivalent only when S is well-conditioned and positive definite.
  S <- sym(S)
  ee <- eigen(S, symmetric = TRUE)
  values <- pmax(ee$values, 0)
  if (length(values) && max(values) > 0) {
    values[values <= tol * max(values)] <- 0
  }
  denom <- A * values^2 + B * values + C
  ridge_scale <- max(abs(denom), 1)
  denom <- denom + ridge * ridge_scale
  weight_values <- ifelse(denom > 0, values / denom, 0)
  sym(ee$vectors %*% (weight_values * t(ee$vectors)))
}

compute_W.gaussian_refined <- function(
  U, V, LK, sig2,
  W1_init, W2_init, W3_init,
  nIter = 1,
  ridge = 1e-8,
  rescale = TRUE
) {
  p <- nrow(U); q <- nrow(V)
  n <- length(LK)

  stopifnot(length(W3_init) == n)

  # constants from your theorem (Gaussian case)
  a1 <- sum(U * U) + p^2
  a2 <- sum(V * V) + q^2
  a3 <- sum(vapply(LK, \(Ki) sum(Ki * Ki) + (sum(diag(Ki)))^2, numeric(1)))

  b1 <- 2 * sig2 * (p^2 + p)
  b2 <- 2 * sig2 * (q^2 + q)
  Nis <- vapply(LK, nrow, integer(1))
  Nbar <- mean(Nis)
  b3 <- 2 * sig2 * sum(vapply(LK, \(Ki) sum(diag(Ki)), numeric(1)) * (Nis + 1))

  c1 <- sig2^2 * (p^2 + p)
  c2 <- sig2^2 * (q^2 + q)
  c3 <- sig2^2 * sum(Nis^2 + Nis)

  # initialize per-target weights from the shared init
  W2U <- W2_init
  W3U <- W3_init

  W1V <- W1_init
  W3V <- W3_init

  W1K <- W1_init
  W2K <- W2_init

  # ----- Target U: alternate (W2, W3i) -----
  for (it in seq_len(nIter)) {
    beta2  <- .beta_mode(V, W2U)
    gamma2 <- .gamma_mode(V, W2U)
    zeta2  <- .zeta_mode(V, W2U)
    eta2   <- .eta_mode(W2U)

    beta3i  <- mapply(\(Ki, W3i) .beta_mode(Ki, W3i), LK, W3U)
    gamma3i <- mapply(\(Ki, W3i) .gamma_mode(Ki, W3i), LK, W3U)
    zeta3i  <- mapply(\(Ki, W3i) .zeta_mode(Ki, W3i),  LK, W3U)
    eta3i   <- vapply(W3U, .eta_mode, numeric(1))

    beta3bar <- mean(beta3i)
    gamma3bar <- mean(gamma3i)
    zeta3bar <- mean(zeta3i)
    eta3bar <- mean(eta3i)

    W2U <- beta3bar * .gaussian_block_weight(
      V,
      A = a1 * gamma3bar,
      B = b1 * zeta3bar,
      C = c1 * eta3bar,
      ridge = ridge
    )

    W3U <- mapply(
      function(Ki) {
        beta2 * .gaussian_block_weight(
          Ki,
          A = a1 * gamma2,
          B = b1 * zeta2,
          C = c1 * eta2,
          ridge = ridge
        )
      },
      LK,
      SIMPLIFY = FALSE
    )

  }

  if (rescale) {
    beta2  <- .beta_mode(V, W2U)
    beta3bar <- mean(mapply(\(Ki, W3i) .beta_mode(Ki, W3i), LK, W3U))
    W2U <- W2U * (q / beta2)
    W3U <- lapply(W3U, function(W3i) W3i * (Nbar / beta3bar))
  }

  # ----- Target V: alternate (W1, W3i) -----
  for (it in seq_len(nIter)) {
    beta1  <- .beta_mode(U, W1V)
    gamma1 <- .gamma_mode(U, W1V)
    zeta1  <- .zeta_mode(U, W1V)
    eta1   <- .eta_mode(W1V)

    beta3i  <- mapply(\(Ki, W3i) .beta_mode(Ki, W3i), LK, W3V)
    gamma3i <- mapply(\(Ki, W3i) .gamma_mode(Ki, W3i), LK, W3V)
    zeta3i  <- mapply(\(Ki, W3i) .zeta_mode(Ki, W3i),  LK, W3V)
    eta3i   <- vapply(W3V, .eta_mode, numeric(1))

    beta3bar <- mean(beta3i)
    gamma3bar <- mean(gamma3i)
    zeta3bar <- mean(zeta3i)
    eta3bar <- mean(eta3i)

    W1V <- beta3bar * .gaussian_block_weight(
      U,
      A = a2 * gamma3bar,
      B = b2 * zeta3bar,
      C = c2 * eta3bar,
      ridge = ridge
    )

    W3V <- mapply(
      function(Ki) {
        beta1 * .gaussian_block_weight(
          Ki,
          A = a2 * gamma1,
          B = b2 * zeta1,
          C = c2 * eta1,
          ridge = ridge
        )
      },
      LK,
      SIMPLIFY = FALSE
    )

  }

  if (rescale) {
    beta1  <- .beta_mode(U, W1V)
    beta3bar <- mean(mapply(\(Ki, W3i) .beta_mode(Ki, W3i), LK, W3V))
    W1V <- W1V * (p / beta1)
    W3V <- lapply(W3V, function(W3i) W3i * (Nbar / beta3bar))
  }

  # ----- Target K: alternate (W1, W2) -----
  for (it in seq_len(nIter)) {
    beta2 <- .beta_mode(V, W2K)
    gamma2 <- .gamma_mode(V, W2K)
    zeta2  <- .zeta_mode(V, W2K)
    eta2   <- .eta_mode(W2K)

    beta1 <- .beta_mode(U, W1K)
    gamma1 <- .gamma_mode(U, W1K)
    zeta1  <- .zeta_mode(U, W1K)
    eta1   <- .eta_mode(W1K)

    W1K <- beta2 * .gaussian_block_weight(
      U,
      A = a3 * gamma2,
      B = b3 * zeta2,
      C = c3 * eta2,
      ridge = ridge
    )
    if (rescale) W1K <- rescale_trace(W1K, p)

    W2K <- beta1 * .gaussian_block_weight(
      V,
      A = a3 * gamma1,
      B = b3 * zeta1,
      C = c3 * eta1,
      ridge = ridge
    )
    if (rescale) W2K <- rescale_trace(W2K, q)
  }

  W <- list(
    U = list(W2 = W2U, W3 = W3U),
    V = list(W1 = W1V, W3 = W3V),
    K = list(W1 = W1K, W2 = W2K)
  )
  class(W) <- c("mode_weights", "gaussian_refined")
  W
}

# Backward-compatible entry point for code written before the method was
# renamed. The weights use Gaussian fourth-moment identities, but are not
# generally optimal after the covariance factors and noise are estimated.
compute_W.optimal <- function(...) {
  .Deprecated("compute_W.gaussian_refined")
  compute_W.gaussian_refined(...)
}

get_W <- function(W, m, target = c("U", "V", "K")) {
  stopifnot(inherits(W, "mode_weights"))
  m <- as.integer(m)[1]
  stopifnot(m %in% 1:3)
  target <- match.arg(target)

  if (inherits(W, "shared")) return(W[[paste0("W", m)]])
  W[[target]][[paste0("W", m)]]
}


###### Update U,V,K ######

.spatial_replication_factor <- function(
  p, q,
  U_tilde = NULL,
  V_tilde = NULL,
  c_delta = 0.5
) {
  p <- as.numeric(p)[1L]
  q <- as.numeric(q)[1L]
  c_delta <- as.numeric(c_delta)[1L]
  if (!is.finite(p) || !is.finite(q) || p < 1 || q < 1) {
    stop("p and q must be positive finite dimensions.")
  }
  fallback <- c_delta * p * q

  if (!is.null(U_tilde) && !is.null(V_tilde)) {
    norm_U2 <- sum(U_tilde * U_tilde)
    norm_V2 <- sum(V_tilde * V_tilde)
    if (is.finite(norm_U2) && is.finite(norm_V2) &&
        norm_U2 > 0 && norm_V2 > 0) {
      fallback <- c_delta * p^2 / norm_U2 * q^2 / norm_V2
    }
  }

  if (!is.finite(fallback) || fallback <= 0) fallback <- 1
  min(max(fallback, 1), p * q)
}

choose.nbasis <- function(
  nbasis = "theory",
  nbasis.min = NULL,
  nbasis.min.K = NULL,
  nbasis.min.G = NULL,
  ndeg = 3,
  ndiag,
  noffdiag,
  nloc,
  Nis = NULL,
  p = NULL,
  q = NULL,
  U_tilde = NULL,
  V_tilde = NULL,
  c_delta = 0.5,
  smooth_order = 2,
  basis_multiplier.K = 3,
  basis_multiplier.G = 1,
  obs_per_basis.K = 8,
  obs_per_basis.G = 8
) {
  # Return separate marginal basis dimensions for the bivariate K smoother
  # and the univariate G smoother. The default rule follows the theoretical
  # rate orders with practical constants, then caps by observed design rows.
  if (is.numeric(nbasis)) {
    if (!length(nbasis) %in% 1:2 || any(!is.finite(nbasis))) {
      stop("Numeric nbasis must have length one or two.")
    }
    out <- rep(as.integer(nbasis), length.out = 2L)
    out <- pmax(as.integer(ndeg) + 2L, pmin(out, as.integer(nloc)))
    return(stats::setNames(out, c("K", "G")))
  }

  nbasis <- match.arg(nbasis, c("theory", "auto", "sparse", "dense"))
  if (nbasis %in% c("auto", "sparse")) nbasis <- "theory"

  nbasis.min <- if (is.null(nbasis.min)) max(5L, as.integer(ndeg) + 2L) else
    as.integer(nbasis.min)
  nbasis.min.K <- if (is.null(nbasis.min.K)) nbasis.min else
    as.integer(nbasis.min.K)
  nbasis.min.G <- if (is.null(nbasis.min.G)) nbasis.min else
    as.integer(nbasis.min.G)
  if (nloc < max(nbasis.min.K, nbasis.min.G)) {
    stop("Too few unique time locations for the requested spline degree.")
  }

  if (nbasis == "theory") {
    if (is.null(Nis)) {
      stop("Automatic theory basis selection requires subject time counts.")
    }
    if (is.null(p) || is.null(q)) {
      stop("Automatic theory basis selection requires p and q.")
    }

    Nis <- as.numeric(Nis)
    if (any(!is.finite(Nis)) || any(Nis < 1)) {
      stop("Subject time counts must be positive and finite.")
    }
    mo <- Nis * (Nis - 1L) / 2
    md <- Nis
    Mo <- sum(mo)
    Md <- sum(md)
    n_o_eff <- if (Mo > 0) Mo^2 / sum(mo^2) else 0
    n_d_eff <- if (Md > 0) Md^2 / sum(md^2) else 0

    r_space <- .spatial_replication_factor(
      p, q,
      U_tilde = U_tilde,
      V_tilde = V_tilde,
      c_delta = c_delta
    )
    N_K_eff <- n_o_eff * r_space
    N_G_eff <- n_d_eff * r_space
    if (min(
      basis_multiplier.K, basis_multiplier.G,
      obs_per_basis.K, obs_per_basis.G
    ) <= 0) {
      stop("Basis multipliers and observation-per-basis caps must be positive.")
    }

    J <- ceiling(basis_multiplier.K * N_K_eff^(1 / (2 * smooth_order + 2)))
    L <- ceiling(basis_multiplier.G * N_G_eff^(1 / (2 * smooth_order)))

    J_cap <- floor(sqrt(noffdiag / obs_per_basis.K))
    L_cap <- floor(ndiag / obs_per_basis.G)
  } else {
    # Dense designs can support a richer basis, while REML still controls
    # the effective degrees of freedom.
    J <- nbasis.min - 1L + ceiling(noffdiag^(1 / (2 * ndeg)))
    L <- nbasis.min - 1L + ceiling(ndiag^(1 / (2 * ndeg - 1)))

    # Keep at least about four observations per tensor-product coefficient.
    J_cap <- max(nbasis.min.K, floor(sqrt(noffdiag / 4)))
    L_cap <- max(nbasis.min.G, floor(ndiag / 4))
  }

  c(
    K = max(
      nbasis.min.K,
      min(as.integer(J), as.integer(J_cap), as.integer(nloc))
    ),
    G = max(
      nbasis.min.G,
      min(as.integer(L), as.integer(L_cap), as.integer(nloc))
    )
  )
}

.fit_temporal_two_step <- function(
  offdiag_data,
  diag_data,
  t_grid,
  nbasis,
  ndeg,
  noise,
  sig2.min
) {
  nt <- length(t_grid)
  covEstFunc <- if (nrow(offdiag_data) > 1e5) mgcv::bam else mgcv::gam
  formula_K <- stats::as.formula(sprintf(
    paste0(
      "Kval ~ te(t1, t2, k = c(%d, %d), ",
      "m = c(%d, %d), bs = c('ps', 'ps'))"
    ),
    nbasis[["K"]], nbasis[["K"]], ndeg, ndeg
  ))

  fit_K <- covEstFunc(
    formula_K,
    data = offdiag_data,
    weights = offdiag_data$w,
    method = "REML"
  )

  grid <- data.table::data.table(
    t1 = rep(t_grid, each = nt),
    t2 = rep(t_grid, times = nt)
  )
  K_hat <- matrix(predict(fit_K, newdata = grid), nt, nt)
  # Estimate sigma^2 against the spline surface before PSD projection. The
  # projection can inflate the fitted diagonal and bias the residual average
  # downward. The outer covariance wrappers project the final K estimate.
  K_hat <- sym(K_hat)

  if (!noise) {
    return(list(
      K = K_hat,
      G = NULL,
      sig2 = 0,
      fit = list(K = fit_K, G = NULL),
      nbasis = nbasis
    ))
  }

  # The diagonal is needed only for the scalar noise variance. Estimate that
  # coefficient directly from diagonal residuals instead of fitting a second
  # univariate smoother. For refined weights, rho is the known fold-specific
  # alpha / beta nugget multiplier; for identity weights, rho = 1.
  nugget_multiplier <- if ("nugget_multiplier" %in% names(diag_data)) {
    diag_data$nugget_multiplier
  } else {
    rep(1, nrow(diag_data))
  }
  diagonal_weight <- diag_data$w
  if (any(!is.finite(nugget_multiplier)) || any(nugget_multiplier <= 0)) {
    stop("Diagonal nugget multipliers must be positive and finite.")
  }
  if (any(!is.finite(diagonal_weight)) || any(diagonal_weight <= 0)) {
    stop("Diagonal averaging weights must be positive and finite.")
  }

  diag_residual <- diag_data$Kval - as.numeric(predict(
    fit_K,
    newdata = data.table::data.table(t1 = diag_data$t, t2 = diag_data$t)
  ))
  denominator <- sum(diagonal_weight * nugget_multiplier^2)
  sig2_unconstrained <- sum(
    diagonal_weight * nugget_multiplier * diag_residual
  ) / denominator
  sig2.min <- if (is.null(sig2.min)) 0 else as.numeric(sig2.min)[1L]
  sig2_hat <- max(sig2.min, sig2_unconstrained)
  if (!is.finite(sig2_hat)) {
    stop("Weighted diagonal residuals did not produce a finite noise estimate.")
  }
  G_hat <- diag(K_hat) + sig2_hat

  list(
    K = K_hat,
    G = G_hat,
    sig2 = sig2_hat,
    fit = list(K = fit_K, G = NULL),
    nbasis = nbasis
  )
}

.fit_temporal_joint <- function(
  scatter_data,
  t_grid,
  nbasis,
  ndeg,
  noise,
  sig2.min
) {
  nt <- length(t_grid)
  covEstFunc <- if (nrow(scatter_data) > 1e5) mgcv::bam else mgcv::gam
  scatter_data <- data.table::copy(scatter_data)
  scatter_data[, diag_ind := as.integer(is_diag)]
  calibrated_nugget <- "nugget_multiplier" %in% names(scatter_data)
  nugget_term <- if (calibrated_nugget) {
    if (any(!is.finite(scatter_data$nugget_multiplier)) ||
        any(scatter_data$nugget_multiplier < 0) ||
        any(scatter_data$nugget_multiplier[scatter_data$is_diag] <= 0)) {
      stop("Nugget multipliers must be positive on the diagonal and finite.")
    }
    "nugget_multiplier"
  } else {
    "diag_ind"
  }

  formula_K <- stats::as.formula(sprintf(
    paste0(
      "Kval ~ te(t1, t2, k = c(%d, %d), ",
      "m = c(%d, %d), bs = c('ps', 'ps'))%s"
    ),
    nbasis[["K"]], nbasis[["K"]], ndeg, ndeg,
    if (noise) paste0(" + ", nugget_term) else ""
  ))

  fit_K <- covEstFunc(
    formula_K,
    data = scatter_data,
    weights = scatter_data$w,
    method = "REML"
  )

  grid <- data.table::data.table(
    t1 = rep(t_grid, each = nt),
    t2 = rep(t_grid, times = nt),
    diag_ind = 0L
  )
  if (calibrated_nugget) grid$nugget_multiplier <- 0
  K_hat <- matrix(predict(fit_K, newdata = grid), nt, nt)
  K_hat <- sym(K_hat)

  if (!noise) {
    return(list(
      K = K_hat,
      G = NULL,
      sig2 = 0,
      fit = list(K = fit_K, G = NULL),
      nbasis = nbasis,
      temporal_method = "joint"
    ))
  }

  sig2.min <- if (is.null(sig2.min)) 0 else as.numeric(sig2.min)[1L]
  delta_hat <- unname(stats::coef(fit_K)[[nugget_term]])
  if (!is.finite(delta_hat)) {
    stop("Joint temporal fit did not produce a finite diagonal noise coefficient.")
  }
  sig2_hat <- max(sig2.min, delta_hat)

  list(
    K = K_hat,
    G = diag(K_hat) + sig2_hat,
    sig2 = sig2_hat,
    fit = list(K = fit_K, G = NULL),
    nbasis = nbasis,
    temporal_method = "joint"
  )
}

estimate_K <- function(
  data,
  W1, W2,
  noise = FALSE,
  nbasis = "theory",
  ndeg = 3,
  sig2.min = NULL,
  temporal_method = c("two_step", "joint"),
  ...
) {
  UseMethod("estimate_K")
}

estimate_K.matfdata.grid <- function(
  data,
  W1, W2,
  noise = FALSE,
  nbasis = "theory",
  ndeg = 3,
  sig2.min = NULL,
  temporal_method = c("two_step", "joint"),
  use_subject_weights = FALSE,
  U_tilde = NULL,
  V_tilde = NULL,
  ...
) {
  if (!is(data, "matfdata.grid")) {
    stop("data needs to be a matfdata.grid")
  }
  temporal_method <- match.arg(temporal_method)

  n <- dim(data$X)[1]
  nt <- length(data$t_grid)
  p <- data$vdim[1]
  q <- data$vdim[2]
  P <- p * q
  ndiag <- n * nt
  noffdiag <- n * nt * (nt - 1L)
  nb <- choose.nbasis(
    nbasis = nbasis,
    ndeg = ndeg,
    ndiag = ndiag,
    noffdiag = noffdiag,
    nloc = nt,
    Nis = rep(nt, n),
    p = p,
    q = q,
    U_tilde = U_tilde,
    V_tilde = V_tilde
  )

  scatter <- foreach::foreach(i = seq_len(n), .combine = rbind) %do% {
    Ki_wk <- get_working_K(get_subj_tnsr(data, i), W1, W2) / P
    mesh <- pracma::meshgrid(data$t_grid, data$t_grid)
    t1 <- as.vector(mesh$X)
    t2 <- as.vector(mesh$Y)
    is_diag <- t1 == t2
    w_off <- if (use_subject_weights) 1 / (nt * (nt - 1L)) else 1
    w_diag <- if (use_subject_weights) 1 / nt else 1
    data.table::data.table(
      t1 = t1,
      t2 = t2,
      Kval = as.vector(Ki_wk),
      is_diag = is_diag,
      w = ifelse(is_diag, w_diag, w_off)
    )
  }

  if (temporal_method == "joint") {
    return(.fit_temporal_joint(
      scatter_data = scatter,
      t_grid = data$t_grid,
      nbasis = nb,
      ndeg = ndeg,
      noise = noise,
      sig2.min = sig2.min
    ))
  }

  out <- .fit_temporal_two_step(
    offdiag_data = scatter[!scatter$is_diag, ],
    diag_data = data.table::data.table(
      t = scatter$t1[scatter$is_diag],
      Kval = scatter$Kval[scatter$is_diag],
      w = scatter$w[scatter$is_diag]
    ),
    t_grid = data$t_grid,
    nbasis = nb,
    ndeg = ndeg,
    noise = noise,
    sig2.min = sig2.min
  )
  out$temporal_method <- "two_step"
  out
}

estimate_K.matfdata.irreg <- function(
  data, W1, W2,
  noise = FALSE,
  nbasis = "theory", ndeg = 3,
  sig2.min = NULL,
  temporal_method = c("two_step", "joint"),
  use_subject_weights = FALSE,
  U_tilde = NULL,
  V_tilde = NULL,
  ...
) {
  stopifnot(inherits(data, "matfdata.irreg"))
  temporal_method <- match.arg(temporal_method)
  n  <- length(data$X_list)
  nt <- length(data$t_grid)
  P  <- prod(data$vdim)
  Nis <- vapply(data$t_idx_list, length, integer(1))

  ndiag <- sum(Nis)
  noffdiag <- sum(Nis * (Nis - 1L))
  nb <- choose.nbasis(
    nbasis = nbasis,
    ndeg = ndeg,
    ndiag = ndiag,
    noffdiag = noffdiag,
    nloc = length(unique(unlist(data$t_idx_list))),
    Nis = Nis,
    p = data$vdim[1],
    q = data$vdim[2],
    U_tilde = U_tilde,
    V_tilde = V_tilde
  )

  # Build separate off-diagonal and diagonal pseudo-observations.
  K_wk_all <- foreach::foreach(i = seq_len(n), .combine = rbind) %do% {
    Ki_wk <- get_working_K(get_subj_tnsr(data, i), W1, W2) / P
    tid <- data$t_idx_list[[i]]
    t_cov_mesh <- pracma::meshgrid(data$t_grid[tid])
    t1 <- as.vector(t_cov_mesh$X)
    t2 <- as.vector(t_cov_mesh$Y)

    dt <- data.table::data.table(
      t1 = t1,
      t2 = t2,
      Kval = as.vector(Ki_wk),
      is_diag = t1 == t2
    )

    if (use_subject_weights) {
      Ni <- length(tid)
      dt$w <- ifelse(dt$is_diag, 1 / Ni, 1 / (Ni * (Ni - 1L)))
    } else {
      dt$w <- 1
    }
    dt
  }

  if (temporal_method == "joint") {
    return(.fit_temporal_joint(
      scatter_data = K_wk_all,
      t_grid = data$t_grid,
      nbasis = nb,
      ndeg = ndeg,
      noise = noise,
      sig2.min = sig2.min
    ))
  }

  out <- .fit_temporal_two_step(
    offdiag_data = K_wk_all[!K_wk_all$is_diag, ],
    diag_data = data.table::data.table(
      t = K_wk_all$t1[K_wk_all$is_diag],
      Kval = K_wk_all$Kval[K_wk_all$is_diag],
      w = K_wk_all$w[K_wk_all$is_diag]
    ),
    t_grid = data$t_grid,
    nbasis = nb,
    ndeg = ndeg,
    noise = noise,
    sig2.min = sig2.min
  )
  out$temporal_method <- "two_step"
  out
}

.estimate_K_crossfit_irreg <- function(
  data,
  W1,
  W2,
  beta,
  alpha = NULL,
  noise = FALSE,
  nbasis = "theory",
  ndeg = 3,
  sig2.min = NULL,
  temporal_method = c("two_step", "joint"),
  crossfit_diagonal = c("refined", "identity"),
  use_subject_weights = FALSE,
  U_tilde = NULL,
  V_tilde = NULL
) {
  stopifnot(inherits(data, "matfdata.irreg"))
  temporal_method <- match.arg(temporal_method)
  crossfit_diagonal <- match.arg(crossfit_diagonal)
  n <- count_subj(data)
  if (length(W1) != n || length(W2) != n || length(beta) != n) {
    stop("Cross-fitted temporal weights and scales must have one entry per subject.")
  }
  beta <- as.numeric(beta)
  if (any(!is.finite(beta)) || any(beta <= 0)) {
    stop("Cross-fitted temporal signal scales must be positive and finite.")
  }
  if (crossfit_diagonal == "refined") {
    alpha <- as.numeric(alpha)
    if (length(alpha) != n || any(!is.finite(alpha)) || any(alpha <= 0)) {
      stop("Cross-fitted temporal noise scales must be positive and finite.")
    }
  }

  nt <- length(data$t_grid)
  P <- prod(data$vdim)
  Nis <- vapply(data$t_idx_list, length, integer(1))
  ndiag <- sum(Nis)
  noffdiag <- sum(Nis * (Nis - 1L))
  nb <- choose.nbasis(
    nbasis = nbasis,
    ndeg = ndeg,
    ndiag = ndiag,
    noffdiag = noffdiag,
    nloc = length(unique(unlist(data$t_idx_list))),
    Nis = Nis,
    p = data$vdim[1],
    q = data$vdim[2],
    U_tilde = U_tilde,
    V_tilde = V_tilde
  )

  # Refined working covariances are divided by beta, which puts their signal
  # on a common K scale. Their diagonal noise contribution is then multiplied
  # by alpha / beta. The legacy strategy replaces the diagonal with identity
  # scatter; the refined strategy retains it and supplies its known multiplier
  # to the temporal noise fit.
  scatter <- foreach::foreach(i = seq_len(n), .combine = rbind) %do% {
    Xi <- get_subj_tnsr(data, i)
    Ki_wk <- get_working_K(Xi, W1[[i]], W2[[i]]) / (P * beta[[i]])
    if (crossfit_diagonal == "identity") {
      diag(Ki_wk) <- apply(Xi@data^2, 3L, sum) / P
    }

    tid <- data$t_idx_list[[i]]
    t_cov_mesh <- pracma::meshgrid(data$t_grid[tid])
    t1 <- as.vector(t_cov_mesh$X)
    t2 <- as.vector(t_cov_mesh$Y)
    is_diag <- t1 == t2
    Ni <- length(tid)
    w_off <- if (use_subject_weights) 1 / (Ni * (Ni - 1L)) else 1
    w_diag <- if (use_subject_weights) 1 / Ni else 1

    dt <- data.table::data.table(
      t1 = t1,
      t2 = t2,
      Kval = as.vector(Ki_wk),
      is_diag = is_diag,
      w = ifelse(is_diag, w_diag, w_off)
    )
    if (crossfit_diagonal == "refined") {
      dt$nugget_multiplier <- ifelse(is_diag, alpha[[i]] / beta[[i]], 0)
    }
    dt
  }

  if (temporal_method == "joint") {
    out <- .fit_temporal_joint(
      scatter_data = scatter,
      t_grid = data$t_grid,
      nbasis = nb,
      ndeg = ndeg,
      noise = noise,
      sig2.min = sig2.min
    )
    out$crossfit_diagonal <- crossfit_diagonal
    return(out)
  }

  diag_data <- data.table::data.table(
    t = scatter$t1[scatter$is_diag],
    Kval = scatter$Kval[scatter$is_diag],
    w = scatter$w[scatter$is_diag]
  )
  if (crossfit_diagonal == "refined") {
    diag_data$nugget_multiplier <-
      scatter$nugget_multiplier[scatter$is_diag]
  }
  out <- .fit_temporal_two_step(
    offdiag_data = scatter[!scatter$is_diag, ],
    diag_data = diag_data,
    t_grid = data$t_grid,
    nbasis = nb,
    ndeg = ndeg,
    noise = noise,
    sig2.min = sig2.min
  )
  out$temporal_method <- "two_step"
  out$crossfit_diagonal <- crossfit_diagonal
  out
}

estimate_U <- function(data, W2, W3, ...) {
  UseMethod("estimate_U")
}

estimate_U.matfdata.grid <- function(data, W2, W3, ...) {
  if (!is(data, "matfdata.grid")) {
    stop("data needs to be a matfdata.grid")
  }

  n <- dim(data$X)[1]
  nt <- length(data$t_grid)
  q <- data$vdim[2]

  U_hat <- foreach::foreach(
    i = seq_len(n),
    .combine = "+"
  ) %do% {
    U_wk <- get_working_U(get_subj_tnsr(data, i), W2, W3) / (q * nt)
    U_wk / n
  }
  U_hat
}

estimate_U.matfdata.irreg <- function(data, W2, W3, ...) {
  if (!is(data, "matfdata.irreg")) {
    stop("data needs to be a matfdata.irreg")
  }

  n <- length(data$X_list)
  nt <- length(data$t_grid)
  q <- data$vdim[2]
  Nis <- sapply(data$t_idx_list, length)
  Nbar <- mean(Nis)

  # average all working U
  U_hat <- foreach::foreach(
    i = seq_len(n),
    .combine = "+"
  ) %do% {
    W2i <- if (is.list(W2)) W2[[i]] else W2
    U_wk <- get_working_U(get_subj_tnsr(data, i), W2i, W3[[i]]) / (q * Nbar)
    U_wk / n
  }
  U_hat
}

estimate_V <- function(data, W1, W3, ...) {
  UseMethod("estimate_V")
}

estimate_V.matfdata.grid <- function(data, W1, W3, ...) {
  if (!is(data, "matfdata.grid")) {
    stop("data needs to be a matfdata.grid")
  }

  n <- dim(data$X)[1]
  nt <- length(data$t_grid)
  p <- data$vdim[1]

  V_hat <- foreach::foreach(
    i = seq_len(n),
    .combine = "+"
  ) %do% {
    V_wk <- get_working_V(get_subj_tnsr(data, i), W1, W3) / (p * nt)
    V_wk / n
  }
  V_hat
}

estimate_V.matfdata.irreg <- function(data, W1, W3, ...) {
  if (!is(data, "matfdata.irreg")) {
    stop("data needs to be a matfdata.irreg")
  }

  n <- length(data$X_list)
  nt <- length(data$t_grid)
  p <- data$vdim[1]
  Nis <- sapply(data$t_idx_list, length)
  Nbar <- mean(Nis)

  # average all working V
  V_hat <- foreach::foreach(
    i = seq_len(n),
    .combine = "+"
  ) %do% {
    W1i <- if (is.list(W1)) W1[[i]] else W1
    V_wk <- get_working_V(get_subj_tnsr(data, i), W1i, W3[[i]]) / (p * Nbar)
    V_wk / n
  }
  V_hat
}

### Main Estimation Functions

.fit_matfcov_identity <- function(
  data,
  noise = FALSE,
  nbasis = "theory",
  ndeg = 3,
  temporal_method = c("two_step", "joint")
) {
  stopifnot(inherits(data, "matfdata.irreg"))
  temporal_method <- match.arg(temporal_method)

  p <- data$vdim[1]
  q <- data$vdim[2]
  t_grid <- data$t_grid
  Nis <- vapply(data$t_idx_list, length, integer(1))

  W1 <- diag(1, p); W2 <- diag(1, q)
  W3 <- lapply(Nis, \(Ni) diag(1, Ni))

  Ksig2 <- estimate_K(
    data, W1, W2,
    noise = noise,
    nbasis = nbasis,
    ndeg = ndeg,
    temporal_method = temporal_method
  )
  K <- psd_project(sym(Ksig2$K))
  sig2 <- Ksig2$sig2

  U_hat <- sym(estimate_U(data, W2, W3))
  Gamma_U <- U_hat - diag(sig2, p)
  U_unprojected <- trace_normalize(Gamma_U, p)
  U <- normalize_psd_trace(Gamma_U, p)
  V_hat <- sym(estimate_V(data, W1, W3))
  Gamma_V <- V_hat - diag(sig2, q)
  V_unprojected <- trace_normalize(Gamma_V, q)
  V <- normalize_psd_trace(Gamma_V, q)

  list(
    C = create_matrix_covfunc(U, V, K, t_grid, sig2),
    W = NULL,
    temporal_fit = Ksig2,
    raw = list(
      Gamma_U = Gamma_U,
      Gamma_V = Gamma_V,
      U_unprojected = U_unprojected,
      V_unprojected = V_unprojected
    )
  )
}

.make_weight_folds <- function(n, nfolds, fold_id = NULL) {
  nfolds_numeric <- as.numeric(nfolds)
  if (length(nfolds_numeric) != 1L || !is.finite(nfolds_numeric) ||
      nfolds_numeric != as.integer(nfolds_numeric)) {
    stop("nfolds must be one finite integer.")
  }
  nfolds <- as.integer(nfolds_numeric)
  if (nfolds < 1L || nfolds > n) {
    stop("nfolds must lie between 1 and the number of subjects.")
  }

  if (is.null(fold_id)) {
    return(rep(seq_len(nfolds), length.out = n))
  }
  if (length(fold_id) != n || anyNA(fold_id)) {
    stop("fold_id must contain one non-missing value per subject.")
  }
  fold_id <- match(fold_id, unique(fold_id))
  if (length(unique(fold_id)) != nfolds) {
    stop("fold_id must contain exactly nfolds distinct folds.")
  }
  as.integer(fold_id)
}

.compute_crossfit_weights <- function(
  data,
  fold_id,
  wScheme,
  noise = FALSE,
  nbasis = "theory",
  ndeg = 3,
  opwIter = 3,
  temporal_method = c("two_step", "joint")
) {
  stopifnot(inherits(data, "matfdata.irreg"))
  temporal_method <- match.arg(temporal_method)
  wScheme <- match_weight_scheme(
    wScheme,
    setdiff(WEIGHT_SCHEMES, "identity")
  )
  n <- count_subj(data)
  p <- data$vdim[1]
  q <- data$vdim[2]
  fold_levels <- sort(unique(fold_id))

  W_subject <- list(
    U = list(W2 = vector("list", n), W3 = vector("list", n)),
    V = list(W1 = vector("list", n), W3 = vector("list", n)),
    K = list(W1 = vector("list", n), W2 = vector("list", n))
  )
  class(W_subject) <- c(
    "crossfit_mode_weights", "mode_weights", wScheme
  )
  attr(W_subject, "fold_id") <- fold_id

  beta_K <- numeric(n)
  alpha_K <- numeric(n)
  pilot_covariances <- vector("list", length(fold_levels))
  training_subjects <- vector("list", length(fold_levels))
  names(pilot_covariances) <- names(training_subjects) <- paste0(
    "fold", fold_levels
  )
  U_tilde <- matrix(0, p, p)
  V_tilde <- matrix(0, q, q)
  K_tilde <- matrix(0, length(data$t_grid), length(data$t_grid))
  sig2_tilde <- 0

  for (fold_position in seq_along(fold_levels)) {
    fold <- fold_levels[[fold_position]]
    test_index <- which(fold_id == fold)
    train_index <- which(fold_id != fold)
    if (!length(train_index)) {
      stop("Each cross-fitting fold must leave at least one training subject.")
    }

    pilot_fit <- .fit_matfcov_identity(
      subset_subj(data, train_index),
      noise = noise,
      nbasis = nbasis,
      ndeg = ndeg,
      temporal_method = temporal_method
    )
    pilot <- pilot_fit$C
    pilot_raw <- pilot_fit$raw
    # Evaluation designs are safe to use: only their time locations enter.
    # No held-out response contributes to these fold-specific weights.
    LK_pilot <- lapply(
      data$t_idx_list,
      function(tid) pilot$K[tid, tid, drop = FALSE]
    )
    fold_weights <- compute_W(
      pilot$U, pilot$V, LK_pilot, pilot$sig2,
      wScheme = wScheme,
      nIter.opt = opwIter
    )

    pilot_U_cal <- if (!is.null(pilot_raw$U_unprojected)) {
      pilot_raw$U_unprojected
    } else {
      pilot$U
    }
    pilot_V_cal <- if (!is.null(pilot_raw$V_unprojected)) {
      pilot_raw$V_unprojected
    } else {
      pilot$V
    }

    W1K <- get_W(fold_weights, 1L, "K")
    W2K <- get_W(fold_weights, 2L, "K")
    fold_beta_K <-
      (tr(pilot_U_cal, W1K) / p) * (tr(pilot_V_cal, W2K) / q)
    if (!is.finite(fold_beta_K) || fold_beta_K <= 0) {
      stop("Cross-fitted pilot produced a non-positive temporal signal scale.")
    }
    fold_alpha_K <- (tr(W1K) / p) * (tr(W2K) / q)
    if (!is.finite(fold_alpha_K) || fold_alpha_K <= 0) {
      stop("Cross-fitted pilot produced a non-positive temporal noise scale.")
    }

    W2U <- get_W(fold_weights, 2L, "U")
    W3U <- get_W(fold_weights, 3L, "U")
    W1V <- get_W(fold_weights, 1L, "V")
    W3V <- get_W(fold_weights, 3L, "V")
    for (i in test_index) {
      W_subject$K$W1[[i]] <- W1K
      W_subject$K$W2[[i]] <- W2K
      beta_K[[i]] <- fold_beta_K
      alpha_K[[i]] <- fold_alpha_K
      W_subject$U$W2[[i]] <- W2U
      W_subject$U$W3[[i]] <- W3U[[i]]
      W_subject$V$W1[[i]] <- W1V
      W_subject$V$W3[[i]] <- W3V[[i]]
    }

    pilot_covariances[[fold_position]] <- pilot
    training_subjects[[fold_position]] <- train_index
    fold_fraction <- length(test_index) / n
    U_tilde <- U_tilde + fold_fraction * pilot$U
    V_tilde <- V_tilde + fold_fraction * pilot$V
    K_tilde <- K_tilde + fold_fraction * pilot$K
    sig2_tilde <- sig2_tilde + fold_fraction * pilot$sig2
  }

  list(
    W = W_subject,
    beta_K = beta_K,
    alpha_K = alpha_K,
    pilot = pilot_covariances,
    training_subjects = training_subjects,
    U_tilde = U_tilde,
    V_tilde = V_tilde,
    K_tilde = K_tilde,
    sig2_tilde = sig2_tilde
  )
}

.fit_matfcov_crossfit <- function(
  data,
  fold_id,
  wScheme,
  noise = FALSE,
  nbasis = "theory",
  ndeg = 3,
  opwIter = 3,
  temporal_method = c("two_step", "joint"),
  crossfit_diagonal = c("refined", "identity")
) {
  stopifnot(inherits(data, "matfdata.irreg"))
  temporal_method <- match.arg(temporal_method)
  crossfit_diagonal <- match.arg(crossfit_diagonal)
  p <- data$vdim[1]
  q <- data$vdim[2]
  t_grid <- data$t_grid
  n <- count_subj(data)
  Nis <- vapply(data$t_idx_list, length, integer(1))
  Nbar <- mean(Nis)

  crossfit <- .compute_crossfit_weights(
    data,
    fold_id = fold_id,
    wScheme = wScheme,
    noise = noise,
    nbasis = nbasis,
    ndeg = ndeg,
    opwIter = opwIter,
    temporal_method = temporal_method
  )
  W <- crossfit$W

  Ksig2 <- .estimate_K_crossfit_irreg(
    data,
    W1 = get_W(W, 1L, "K"),
    W2 = get_W(W, 2L, "K"),
    beta = crossfit$beta_K,
    alpha = crossfit$alpha_K,
    noise = noise,
    nbasis = nbasis,
    ndeg = ndeg,
    temporal_method = temporal_method,
    crossfit_diagonal = crossfit_diagonal,
    U_tilde = crossfit$U_tilde,
    V_tilde = crossfit$V_tilde
  )
  K_uncalibrated <- psd_project(sym(Ksig2$K))
  K_scale <- fun2dInprod(
    K_uncalibrated, crossfit$K_tilde, t_grid
  ) / fun2dNorm2(K_uncalibrated, t_grid)
  if (!is.finite(K_scale) || K_scale <= 0) {
    stop("Fold-pilot temporal calibration produced a non-positive scale.")
  }
  K <- K_scale * K_uncalibrated
  # The legacy mixed strategy keeps its stable identity-pilot noise estimate.
  # With refined diagonals, the temporal fit accounts for each fold's known
  # alpha / beta nugget multiplier and directly estimates the common sig2.
  sig2 <- if (!noise) {
    0
  } else if (crossfit_diagonal == "refined") {
    max(0, Ksig2$sig2)
  } else {
    max(0, crossfit$sig2_tilde)
  }
  Ksig2$crossfit_scale <- K_scale
  Ksig2$sig2_uncalibrated <- Ksig2$sig2
  Ksig2$sig2 <- sig2
  Ksig2$K_calibrated <- K
  if (crossfit_diagonal == "refined" || temporal_method == "joint") {
    Ksig2$G <- diag(K) + sig2
  }

  W2U <- get_W(W, 2L, "U")
  W3U <- get_W(W, 3L, "U")
  U_hat <- sym(estimate_U(data, W2U, W3U))
  noise_scale_U <- mean(vapply(seq_len(n), function(i) {
    tr(W2U[[i]]) * tr(W3U[[i]]) / (q * Nbar)
  }, numeric(1)))
  Gamma_U <- U_hat - diag(sig2 * noise_scale_U, p)
  U_unprojected <- trace_normalize(Gamma_U, p)
  U <- normalize_psd_trace(Gamma_U, p)

  W1V <- get_W(W, 1L, "V")
  W3V <- get_W(W, 3L, "V")
  V_hat <- sym(estimate_V(data, W1V, W3V))
  noise_scale_V <- mean(vapply(seq_len(n), function(i) {
    tr(W1V[[i]]) * tr(W3V[[i]]) / (p * Nbar)
  }, numeric(1)))
  Gamma_V <- V_hat - diag(sig2 * noise_scale_V, q)
  V_unprojected <- trace_normalize(Gamma_V, q)
  V <- normalize_psd_trace(Gamma_V, q)

  list(
    C = create_matrix_covfunc(U, V, K, t_grid, sig2),
    W = W,
    pilot = crossfit$pilot,
    beta_K = crossfit$beta_K,
    alpha_K = crossfit$alpha_K,
    fold_id = fold_id,
    training_subjects = crossfit$training_subjects,
    temporal_fit = Ksig2,
    raw = list(
      Gamma_U = Gamma_U,
      Gamma_V = Gamma_V,
      U_unprojected = U_unprojected,
      V_unprojected = V_unprojected
    )
  )
}

.fit_matfcov_from_pilot <- function(
  data,
  pilot,
  pilot_raw = NULL,
  wScheme,
  noise = FALSE,
  nbasis = "theory",
  ndeg = 3,
  opwIter = 3,
  temporal_method = c("two_step", "joint")
) {
  stopifnot(inherits(data, "matfdata.irreg"), is(pilot, "matfcov"))
  temporal_method <- match.arg(temporal_method)
  wScheme <- match_weight_scheme(
    wScheme,
    setdiff(WEIGHT_SCHEMES, "identity")
  )

  p <- data$vdim[1]
  q <- data$vdim[2]
  t_grid <- data$t_grid
  Nis <- vapply(data$t_idx_list, length, integer(1))
  Nbar <- mean(Nis)
  LK_pilot <- lapply(
    data$t_idx_list,
    \(tid) pilot$K[tid, tid, drop = FALSE]
  )
  pilot_U_cal <- if (!is.null(pilot_raw$U_unprojected)) {
    pilot_raw$U_unprojected
  } else {
    pilot$U
  }
  pilot_V_cal <- if (!is.null(pilot_raw$V_unprojected)) {
    pilot_raw$V_unprojected
  } else {
    pilot$V
  }

  W <- compute_W(
    pilot$U, pilot$V, LK_pilot, pilot$sig2,
    wScheme = wScheme,
    nIter.opt = opwIter
  )

  # Temporal update: pilot covariances calibrate the unavailable beta factors.
  W1K <- get_W(W, 1, "K")
  W2K <- get_W(W, 2, "K")
  Ksig2 <- estimate_K(
    data, W1K, W2K,
    noise = noise,
    nbasis = nbasis,
    ndeg = ndeg,
    temporal_method = temporal_method,
    U_tilde = pilot$U,
    V_tilde = pilot$V
  )
  beta1K <- tr(pilot_U_cal, W1K) / p
  beta2K <- tr(pilot_V_cal, W2K) / q
  alpha1K <- tr(W1K) / p
  alpha2K <- tr(W2K) / q
  if (min(beta1K, beta2K, alpha1K, alpha2K) <= 0) {
    stop("Pilot calibration produced a non-positive temporal scaling factor.")
  }

  K <- psd_project(sym(Ksig2$K) / (beta1K * beta2K))
  sig2 <- if (noise) Ksig2$sig2 / (alpha1K * alpha2K) else 0

  # Spatial updates use the same frozen pilot weights. Trace normalization
  # removes the unknown signal scaling factors after de-ridging.
  W2U <- get_W(W, 2, "U")
  W3U <- get_W(W, 3, "U")
  alpha2U <- tr(W2U) / q
  alpha3Ubar <- mean(vapply(W3U, tr, numeric(1))) / Nbar
  U_hat <- sym(estimate_U(data, W2U, W3U))
  Gamma_U <- U_hat - diag(sig2 * alpha2U * alpha3Ubar, p)
  U_unprojected <- trace_normalize(Gamma_U, p)
  U <- normalize_psd_trace(Gamma_U, p)

  W1V <- get_W(W, 1, "V")
  W3V <- get_W(W, 3, "V")
  alpha1V <- tr(W1V) / p
  alpha3Vbar <- mean(vapply(W3V, tr, numeric(1))) / Nbar
  V_hat <- sym(estimate_V(data, W1V, W3V))
  Gamma_V <- V_hat - diag(sig2 * alpha1V * alpha3Vbar, q)
  V_unprojected <- trace_normalize(Gamma_V, q)
  V <- normalize_psd_trace(Gamma_V, q)

  list(
    C = create_matrix_covfunc(U, V, K, t_grid, sig2),
    W = W,
    pilot = pilot,
    temporal_fit = Ksig2,
    raw = list(
      Gamma_U = Gamma_U,
      Gamma_V = Gamma_V,
      U_unprojected = U_unprojected,
      V_unprojected = V_unprojected
    )
  )
}

fit_matfcov <- function(
  data,
  wScheme = WEIGHT_SCHEMES,
  noise = FALSE,
  nbasis = "theory",
  ndeg = 3,
  opwIter = 3,
  temporal_method = c("two_step", "joint"),
  nfolds = 10L,
  fold_id = NULL,
  crossfit_diagonal = c("refined", "identity")
) {
  wScheme <- match_weight_scheme(wScheme)
  temporal_method <- match.arg(temporal_method)
  crossfit_diagonal <- match.arg(crossfit_diagonal)
  stopifnot(inherits(data, "matfdata.irreg"))

  if (wScheme == "identity") {
    return(.fit_matfcov_identity(
      data,
      noise = noise,
      nbasis = nbasis,
      ndeg = ndeg,
      temporal_method = temporal_method
    ))
  }

  if (wScheme == "gaussian_refined") {
    if (missing(nfolds)) nfolds <- min(as.integer(nfolds), count_subj(data))
    fold_id <- .make_weight_folds(count_subj(data), nfolds, fold_id)
    if (nfolds > 1L) {
      return(.fit_matfcov_crossfit(
        data,
        fold_id = fold_id,
        wScheme = wScheme,
        noise = noise,
        nbasis = nbasis,
        ndeg = ndeg,
        opwIter = opwIter,
        temporal_method = temporal_method,
        crossfit_diagonal = crossfit_diagonal
      ))
    }
  }

  # Gaussian nfolds = 1 and heuristic weights use the legacy same-sample fit.
  pilot_fit <- .fit_matfcov_identity(
    data,
    noise = noise,
    nbasis = nbasis,
    ndeg = ndeg,
    temporal_method = temporal_method
  )
  .fit_matfcov_from_pilot(
    data,
    pilot = pilot_fit$C,
    pilot_raw = pilot_fit$raw,
    wScheme = wScheme,
    noise = noise,
    nbasis = nbasis,
    ndeg = ndeg,
    opwIter = opwIter,
    temporal_method = temporal_method
  )
}


get_matfpca_scores <- function(
  data,
  eig,
  sig2 = 0,
  K_grid = NULL,
  nugget = 1e-8
) {
  # Best linear prediction of retained product-FPC scores.
  # Under Gaussian scores/errors this is also the conditional expectation.
  n <- count_subj(data)
  I <- ncol(eig$funcs$Phi)
  J <- ncol(eig$funcs$Psi)
  K <- ncol(eig$funcs$Theta)

  lam_phi <- eig$values$lam_phi
  lam_psi <- eig$values$lam_psi
  lam_theta <- eig$values$lam_theta
  stopifnot(length(lam_phi) == I, length(lam_psi) == J, length(lam_theta) == K)
  stopifnot(is.numeric(sig2), length(sig2) == 1L, is.finite(sig2), sig2 >= 0)
  stopifnot(is.numeric(nugget), length(nugget) == 1L, is.finite(nugget), nugget >= 0)

  # The prediction denominator should use the full temporal covariance, not
  # only the retained-rank numerator. This accounts for omitted temporal modes.
  if (is.null(K_grid)) {
    K_grid <- eig$funcs$Theta %*%
      t(sweep(eig$funcs$Theta, 2, lam_theta, "*"))
  }
  stopifnot(
    is.matrix(K_grid),
    nrow(K_grid) == length(data$t_grid),
    ncol(K_grid) == length(data$t_grid)
  )
  K_grid <- sym(K_grid)

  scores <- rTensor::as.tensor(array(NA_real_, dim = c(n, I, J, K)))

  scores_list <- foreach::foreach(
    i = seq_len(n),
    .combine = c,
    .multicombine = TRUE
  ) %do% {
    tid <- data$t_idx_list[[i]]
    Thetai <- eig$funcs$Theta[tid, , drop = FALSE]
    mi <- length(tid)
    Ki <- sym(K_grid[tid, tid, drop = FALSE])

    # Spatial projection gives an I x J x mi tensor of visit vectors.
    projected <- rTensor::ttl(
      get_subj_tnsr(data, i),
      list(t(eig$funcs$Phi), t(eig$funcs$Psi)),
      ms = c(1, 2)
    )@data

    # Diagonalize Ki once; each spatial pair differs only through the ridge
    # sig2 / (lambda_U,j lambda_V,k).
    ee_Ki <- eigen(Ki, symmetric = TRUE)
    Q <- ee_Ki$vectors
    d <- pmax(ee_Ki$values, 0)
    Qt <- t(Q)
    score_numerator <- t(sweep(Thetai, 2, lam_theta, "*"))

    scores_i <- array(0, dim = c(I, J, K))
    for (j in seq_len(I)) {
      for (k in seq_len(J)) {
        spatial_signal <- lam_phi[[j]] * lam_psi[[k]]
        if (!is.finite(spatial_signal) || spatial_signal <= 0) next

        ridge_jk <- max(nugget, sig2 / spatial_signal)
        z_jk <- projected[j, k, ]
        Ainv_z <- Q %*% ((Qt %*% z_jk) / (d + ridge_jk))
        scores_i[j, k, ] <- as.numeric(score_numerator %*% Ainv_z)
      }
    }
    list(rTensor::as.tensor(scores_i))
  }

  for (i in seq_len(n)) {
    scores[i, , , ] <- scores_list[[i]]
  }

  scores
}

recover_matfdata <- function(eigf, scores, mean_func = NULL) {
  # validate arguments
  I <- ncol(eigf$Phi)
  J <- ncol(eigf$Psi)
  K <- ncol(eigf$Theta)
  stopifnot(length(dim(scores)) == 4)
  stopifnot(all(dim(scores)[2:4] == c(I, J, K)))
  if (is.array(scores)) {
    scores <- rTensor::as.tensor(scores)
  }

  X <- with(eigf, ttl(scores, list(Phi, Psi, Theta), 2:4))
  X
}


# Evaluation Metrics -----------------------------------------------------------

quad2d_trap <- function(x, K) {
  stopifnot(is.numeric(x), is.matrix(K), nrow(K) == length(x), ncol(K) == length(x))
  n <- length(x)
  if (is.unsorted(x, strictly = FALSE)) stop("x must be sorted (nondecreasing).")
  if (anyDuplicated(x)) stop("x must have unique points for trapezoid weights.")

  dx <- diff(x)

  # 1D composite trapezoid weights on irregular grid
  w <- numeric(n)
  w[1] <- dx[1] / 2
  w[n] <- dx[n - 1] / 2
  if (n > 2) w[2:(n - 1)] <- (dx[1:(n - 2)] + dx[2:(n - 1)]) / 2

  # 2D integral ≈ w^T K w
  as.numeric(crossprod(w, K %*% w))
}

matFNorm2 <- function(X) sum(X^2)

matFInprod <- function(X, Y) sum(X * Y)

fun2dNorm2 <- function(f, x) quad2d_trap(x, f^2)

fun2dInprod <- function(f, g, x) quad2d_trap(x, f * g)

matfcovNorm2 <- function(cobj) {
  stopifnot(is(cobj, "matfcov"))
  matFNorm2(cobj$U) * matFNorm2(cobj$V) * fun2dNorm2(cobj$K, cobj$t_grid)
}

matfcovInprod <- function(cobj1, cobj2) {
  stopifnot(is(cobj1, "matfcov"))
  stopifnot(is(cobj2, "matfcov"))
  stopifnot(all.equal(cobj1$t_grid, cobj2$t_grid))
  matFInprod(cobj1$U, cobj2$U) * matFInprod(cobj1$V, cobj2$V) * fun2dInprod(cobj1$K, cobj2$K, cobj1$t_grid)
}

matfcovDiffNorm2 <- function(cobj1, cobj2) {
  stopifnot(is(cobj1, "matfcov"))
  stopifnot(is(cobj2, "matfcov"))
  matfcovNorm2(cobj1) + matfcovNorm2(cobj2) - 2 * matfcovInprod(cobj1, cobj2)
}
