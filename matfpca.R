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

  Uinv <- ginv_spd(U, ridge = ridge)
  Vinv <- ginv_spd(V, ridge = ridge)
  LKinv <- lapply(LK, \(Ki) ginv_spd(Ki, ridge = ridge))

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

    W2U <- beta3bar * inv_spd(
      a1 * gamma3bar * V +
        b1 * zeta3bar * diag(1, q) +
        c1 * eta3bar * Vinv,
      ridge = ridge
    )

    W3U <- mapply(
      function(Ki, Ki_inv) {
        beta2 * inv_spd(
          a1 * gamma2 * Ki +
            b1 * zeta2 * diag(1, nrow(Ki)) +
            c1 * eta2 * Ki_inv,
          ridge = ridge
        )
      },
      LK, LKinv,
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

    W1V <- beta3bar * inv_spd(
      a2 * gamma3bar * U +
        b2 * zeta3bar * diag(1, p) +
        c2 * eta3bar * Uinv,
      ridge = ridge
    )

    W3V <- mapply(
      function(Ki, Ki_inv) {
        beta1 * inv_spd(a2 * gamma1 * Ki +
            b2 * zeta1 * diag(1, nrow(Ki)) +
            c2 * eta1 * Ki_inv,
          ridge = ridge
        )
      },
      LK, LKinv,
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

    W1K <- beta2 * inv_spd(a3 * gamma2 * U + b3 * zeta2 * diag(1, p) + c3 * eta2 * Uinv, ridge = ridge)
    if (rescale) W1K <- rescale_trace(W1K, p)

    W2K <- beta1 * inv_spd(a3 * gamma1 * V + b3 * zeta1 * diag(1, q) + c3 * eta1 * Vinv, ridge = ridge)
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
  # Estimate sigma^2 from the spline surface before PSD projection. The
  # projection can inflate the fitted diagonal and bias G(t) - K(t, t)
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

  diagEstFunc <- if (nrow(diag_data) > 1e5) mgcv::bam else mgcv::gam
  formula_G <- stats::as.formula(sprintf(
    "Kval ~ s(t, k = %d, m = %d, bs = 'ps')",
    nbasis[["G"]], ndeg
  ))
  fit_G <- diagEstFunc(
    formula_G,
    data = diag_data,
    weights = diag_data$w,
    method = "REML"
  )
  G_hat <- as.numeric(predict(
    fit_G,
    newdata = data.table::data.table(t = t_grid)
  ))

  sig2.min <- if (is.null(sig2.min)) 0 else as.numeric(sig2.min)[1L]
  sig2_hat <- max(
    sig2.min,
    trapz_average(t_grid, G_hat - diag(K_hat))
  )

  list(
    K = K_hat,
    G = G_hat,
    sig2 = sig2_hat,
    fit = list(K = fit_K, G = fit_G),
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

  formula_K <- stats::as.formula(sprintf(
    paste0(
      "Kval ~ te(t1, t2, k = c(%d, %d), ",
      "m = c(%d, %d), bs = c('ps', 'ps'))%s"
    ),
    nbasis[["K"]], nbasis[["K"]], ndeg, ndeg,
    if (noise) " + diag_ind" else ""
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
  delta_hat <- unname(stats::coef(fit_K)[["diag_ind"]])
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
    U_wk <- get_working_U(get_subj_tnsr(data, i), W2, W3[[i]]) / (q * Nbar)
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
    V_wk <- get_working_V(get_subj_tnsr(data, i), W1, W3[[i]]) / (p * Nbar)
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
  U <- normalize_psd_trace(U_hat - diag(sig2, p), p)
  V_hat <- sym(estimate_V(data, W1, W3))
  V <- normalize_psd_trace(V_hat - diag(sig2, q), q)

  list(
    C = create_matrix_covfunc(U, V, K, t_grid, sig2),
    W = NULL,
    temporal_fit = Ksig2
  )
}

.fit_matfcov_from_pilot <- function(
  data,
  pilot,
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
  beta1K <- tr(pilot$U, W1K) / p
  beta2K <- tr(pilot$V, W2K) / q
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
  U <- normalize_psd_trace(
    U_hat - diag(sig2 * alpha2U * alpha3Ubar, p),
    p
  )

  W1V <- get_W(W, 1, "V")
  W3V <- get_W(W, 3, "V")
  alpha1V <- tr(W1V) / p
  alpha3Vbar <- mean(vapply(W3V, tr, numeric(1))) / Nbar
  V_hat <- sym(estimate_V(data, W1V, W3V))
  V <- normalize_psd_trace(
    V_hat - diag(sig2 * alpha1V * alpha3Vbar, q),
    q
  )

  list(
    C = create_matrix_covfunc(U, V, K, t_grid, sig2),
    W = W,
    pilot = pilot,
    temporal_fit = Ksig2
  )
}

fit_matfcov <- function(
  data,
  wScheme = WEIGHT_SCHEMES,
  noise = FALSE,
  nbasis = "theory",
  ndeg = 3,
  opwIter = 3,
  temporal_method = c("two_step", "joint")
) {
  wScheme <- match_weight_scheme(wScheme)
  temporal_method <- match.arg(temporal_method)
  stopifnot(inherits(data, "matfdata.irreg"))

  pilot_fit <- .fit_matfcov_identity(
    data,
    noise = noise,
    nbasis = nbasis,
    ndeg = ndeg,
    temporal_method = temporal_method
  )
  if (wScheme == "identity") return(pilot_fit)

  .fit_matfcov_from_pilot(
    data,
    pilot = pilot_fit$C,
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
  nugget = 1e-8
) {
  n <- count_subj(data)
  I <- ncol(eig$funcs$Phi)
  J <- ncol(eig$funcs$Psi)
  K <- ncol(eig$funcs$Theta)

  scores <- rTensor::as.tensor(array(NA, dim = c(n, I, J, K)))
  lam_theta <- eig$values$lam_theta

  scores_list <- foreach::foreach(
    i = seq_len(n),
    .combine = c,
    .multicombine = TRUE
  ) %do%
    {
      Thetai <- eig$funcs$Theta[data$t_idx_list[[i]], , drop = F]
      mi <- length(data$t_idx_list[[i]])
      Ki_inv <- solve(
        Thetai %*% t(sweep(Thetai, 2, lam_theta, "*")) + diag(nugget, mi, mi)
      )
      scores_i <- rTensor::ttl(
        get_subj_tnsr(data, i),
        list(
          t(eig$funcs$Phi),
          t(eig$funcs$Psi),
          t(sweep(Thetai, 2, lam_theta, "*")) %*% Ki_inv
        ),
        1:3
      )
      list(scores_i)
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
