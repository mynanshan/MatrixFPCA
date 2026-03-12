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

WEIGHT_SCHEMES <- c("identity", "invcov", "invridcov", "optimal") 

compute_W <- function(
  U, V, LK,
  sig2 = 0,
  wScheme = WEIGHT_SCHEMES,
  nIter.opt = 1,
  rescale.opt = TRUE,
  ridge.opt = 1e-8
) {
  wScheme <- match.arg(wScheme)
  eps <- 10 * .Machine$double.eps

  if (sig2 <= eps && wScheme %in% c("optimal", "invridcov")) {
    wScheme <- "invcov"
  }

  if (wScheme == "optimal") {
    W0 <- compute_W.shared(U, V, LK, sig2, "invridcov", eps)
    return(compute_W.optimal(
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

compute_W.optimal <- function(
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

    # if (rescale) {
    #   W2U <- rescale_trace(W2U, q)
    #   mtrW3U <- mean(vapply(W3U, tr, numeric(1)))
    #   if (is.finite(mtrW3U) && mtrW3U > 1e-14) {
    #     s <- Nbar / mtrW3U
    #     W3U <- lapply(W3U, function(W3i) W3i * s)
    #   }
    # }
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

    # if (rescale) {
    #   W1V <- rescale_trace(W1V, p)
    #   mtrW3V <- mean(vapply(W3V, tr, numeric(1)))
    #   if (is.finite(mtrW3V) && mtrW3V > 1e-14) {
    #     s <- Nbar / mtrW3V
    #     W3V <- lapply(W3V, function(W3i) W3i * s)
    #   }
    # }
  }

  if (rescale) {
    beta1  <- .beta_mode(V, W1V)
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
  class(W) <- c("mode_weights", "optimal")
  W
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

choose.nbasis <- function(
  nbasis = "sparse",
  nbasis.min = 5,
  ndeg = 3,
  nsubj = NULL,
  nobs = NULL,
  nloc = NULL
) {
  # TODO: fix basis choices
  ## Arguments
  # nbasis: Assume smoothness >= ndeg
  #   (1) numeric vector, return nbasis as is
  #   (2) "sparse"
  #   (3) "dense"
  if (is.numeric(nbasis)) {
    return(as.integer(nbasis))
  }
  if (nbasis == "auto") nbasis <- "sparse"
  if (nbasis == "sparse") {
    Nmean <- nobs / nsubj
    nb <- nbasis.min + ceiling((nsubj * Nmean^2)^(1 / (2 * ndeg + 2))) - 1
  }
  if (nbasis == "dense") {
    Nmean <- nobs / nsubj
    scaling <- Nmean / 2 # an intuitive scaling
    nb <- nbasis.min + ceiling(scaling * nsubj^(1 / (2 * ndeg))) - 1
  }
  # restrict the number of marginal basis functions to be
  # smaller than number of unique observation locations
  min(nb, nloc)
}

estimate_K <- function(
  data,
  W1, W2,
  noise = FALSE,
  nbasis = "auto",
  ndeg = 3,
  sig2.min = NULL,
  ...
) {
  UseMethod("estimate_K")
}

estimate_K.matfdata.grid <- function(
  data,
  W1, W2,
  noise = FALSE,
  nbasis = "auto",
  ndeg = 3,
  sig2.min = NULL,
  ...
) {
  if (!is(data, "matfdata.grid")) {
    stop("data needs to be a matfdata.grid")
  }

  n <- dim(data$X)[1]
  nt <- length(data$t_grid)
  p <- data$vdim[1]
  q <- data$vdim[2]
  P <- p * q

  K_hat <- foreach::foreach(
    i = seq_len(n),
    .combine = "+"
  ) %do% {
    # get the working K
    Ki_wk <- get_working_K(get_subj_tnsr(data, i), W1, W2) / P
    Ki_wk / n
  }

  if (!noise) {
    return(list(K = K_hat, sig2 = 0))
  }

  # create meshgrid
  t_cov_mesh <- pracma::meshgrid(data$t_grid, data$t_grid)
  t1 <- as.vector(t_cov_mesh$X)
  t2 <- as.vector(t_cov_mesh$Y)

  # num.basis is almost num.points
  nb <- ifelse(is.null(nbasis) || !is.numeric(nbasis), nt - 2, nbasis)

  # TODO: also switch to joint est of K and sig2

  # smooth off-diagonal values
  offdiag_idx <- t1 != t2
  K_wk_offdiag <- data.table::data.table(
    t1 = t1[offdiag_idx],
    t2 = t2[offdiag_idx],
    Kval = K_hat[offdiag_idx]
  )
  covEstFunc <- if(nrow(K_wk_offdiag) > 1e5) mgcv::bam else mgcv::gam
  fit_K <- covEstFunc(
    Kval ~ te(t1, t2, k = c(nb, nb), m = ndeg, bs = "ps"),
    data = K_wk_offdiag
  )
  K_smooth <- predict(
    fit_K,
    newdata = data.table::data.table(t1 = t1, t2 = t2)
  ) |>
    matrix(nt, nt) |> 
    psd_project()

  # smooth diagonal values
  fit_K_diag <- mgcv::gam(
    Kval ~ s(t, k = nb, m = ndeg, bs = "ps"),
    data = data.table::data.table(
      t = data$t_grid,
      Kval = diag(K_hat)
    )
  )
  K_diag <- predict(
    fit_K_diag,
    newdata = data.table::data.table(t = data$t_grid)
  )

  # noise level estimate
  sig2.min <- ifelse(
    is.null(sig2.min),
    1e-4 * mean(K_smooth^2),
    as.numeric(sig2.min)
  )
  K_diag_smooth <- K_smooth[t1 == t2]
  mid_idx <- seq(floor(nt * 0.1), ceiling(nt * 0.9))
  sig2 <- max(sig2.min, mean(K_diag[mid_idx] - K_diag_smooth[mid_idx]))

  # return K and sig2
  list(K = K_smooth, sig2 = sig2)
}

estimate_K.matfdata.irreg <- function(
  data, W1, W2,
  nbasis = "auto", ndeg = 3,
  sig2.min = 0,
  use_subject_weights = FALSE,
  ...
) {
  stopifnot(inherits(data, "matfdata.irreg"))
  n  <- length(data$X_list)
  nt <- length(data$t_grid)
  P  <- prod(data$vdim)
  Nis <- vapply(data$t_idx_list, length, integer(1))

  nb <- choose.nbasis(
    nbasis = nbasis,
    ndeg = ndeg,
    nsubj = n,
    nobs = sum(Nis),
    nloc = length(unique(unlist(data$t_idx_list)))
  )

  # Build scatter table of all (t1,t2,Kval) points
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
      diag = as.integer(t1 == t2)
    )

    if (use_subject_weights) {
      # prevent big-Ni subjects from dominating (one simple choice)
      dt$w <- 1 / (length(tid)^2)
    } else {
      dt$w <- 1
    }
    dt
  }

  covEstFunc <- if (nrow(K_wk_all) > 1e5) mgcv::bam else mgcv::gam

  # Joint fit: smooth f(t1,t2) + parametric diag term
  fit <- covEstFunc(
    Kval ~ te(t1, t2, k = c(nb, nb), m = c(ndeg, ndeg), bs = c("ps", "ps")) + diag,
    data = K_wk_all,
    weights = K_wk_all$w,
    method = "REML"
  )

  # Predict the smooth covariance surface with diag=0 (so nugget not included)
  grid <- data.table::data.table(
    t1 = rep(data$t_grid, each = nt),
    t2 = rep(data$t_grid, times = nt),
    diag = 0L
  )
  K_hat <- matrix(predict(fit, newdata = grid), nt, nt)
  K_hat <- psd_project(sym(K_hat))

  # diag coefficient estimates delta = sigma^2 * alpha1 * alpha2 on the K_raw scale
  delta_hat <- max(unname(coef(fit)["diag"]), sig2.min)

  list(K = K_hat, sig2 = delta_hat, fit = fit)
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

### Main Estimation Function

fit_matfcov <- function(
  data,
  wScheme = WEIGHT_SCHEMES,
  noise = FALSE,
  nbasis = "auto",
  ndeg = 3,
  maxIter = 2,
  opwIter = 3
) {
  wScheme <- match.arg(wScheme)
  stopifnot(inherits(data, "matfdata.irreg"))

  n <- count_subj(data)
  p <- data$vdim[1]
  q <- data$vdim[2]
  t_grid <- data$t_grid
  nt <- length(t_grid)

  # TODO: some lines are still specific to matfdata.irreg. fix that in the future
  Nis <- sapply(data$t_idx_list, length)
  Nbar <- mean(Nis)

  #==== First iteration: identity weights ====#

  W1 <- diag(1, p); W2 <- diag(1, q)
  W3 <- lapply(Nis, \(Ni) diag(1, Ni))

  # ----- update K -----
  Ksig2 <- estimate_K(data, W1, W2, noise = noise, nbasis = nbasis, ndeg = ndeg)
  K <- sym(Ksig2$K)
  sig2 <- Ksig2$sig2

  LK <- lapply(data$t_idx_list, function(tid) K[tid, tid, drop = FALSE])

  beta3bar <- mean(mapply(function(Ki, W3i) tr(Ki, W3i), LK, W3)) / Nbar

  # ----- update U -----
  U_hat <- sym(estimate_U(data, W2, W3))

  U_tilde <- U_hat - diag(sig2, p)
  U_tilde <- psd_project(U_tilde)

  U <- U_tilde / beta3bar
  U <- p * U / tr(U)

  # ----- update V -----
  V_hat <- sym(estimate_V(data, W1, W3))

  V_tilde <- V_hat - diag(sig2, q)
  V_tilde <- psd_project(V_tilde)

  V <- V_tilde / beta3bar
  V <- q * V / tr(V)

  #==== Rest iterations: wScheme ====#

  for (iter in seq_len(max(maxIter - 1, 0))) {

    W <- compute_W(U, V, LK, sig2, wScheme, nIter.opt = opwIter)
    # FIXME: cheating code
    # W <- compute_W(
    #   cov_true$U, cov_true$V, lapply(data$t_idx_list, function(tid) cov_true$K[tid, tid, drop = FALSE]),
    #   cov_true$sig2, wScheme, nIter.opt = opwIter
    # )

    # ----- update K -----
    W1K <- get_W(W, 1, "K"); W2K <- get_W(W, 2, "K")
    Ksig2 <- estimate_K(data, W1K, W2K, noise = noise, nbasis = nbasis, ndeg = ndeg)
    K_raw <- sym(Ksig2$K)
    sig2_raw <- Ksig2$sig2

    beta1K <- tr(U, W1K) / p
    beta2K <- tr(V, W2K) / q
    alpha1K <- tr(W1K) / p
    alpha2K <- tr(W2K) / q

    K <- psd_project(K_raw / (beta1K * beta2K))
    sig2 <- sig2_raw / (alpha1K * alpha2K)

    # refresh LK after K update
    LK <- lapply(data$t_idx_list, function(tid) K[tid, tid, drop = FALSE])

    # ----- update U -----
    W2U <- get_W(W, 2, "U"); W3U <- get_W(W, 3, "U")
    U_hat <- sym(estimate_U(data, W2U, W3U))

    beta2U <- tr(V, W2U) / q
    beta3Ubar <- mean(mapply(function(Ki, W3i) tr(Ki, W3i), LK, W3U)) / Nbar

    alpha2U <- tr(W2U) / q
    alpha3Ubar <- mean(sapply(W3U, tr)) / Nbar

    sig2_cap <- min(eigen(U_hat, symmetric=TRUE, only.values=TRUE)$values)
    sig2_use <- max(0, min(sig2, 0.99 * sig2_cap / (alpha2U * alpha3Ubar + 1e-12)))
    U_tilde <- U_hat - diag(sig2_use * alpha2U * alpha3Ubar, p)
    U_tilde <- psd_project(U_tilde)

    U <- U_tilde / (beta2U * beta3Ubar)
    U <- p * U / tr(U)

    # ----- update V -----
    W1V <- get_W(W, 1, "V"); W3V <- get_W(W, 3, "V")
    V_hat <- sym(estimate_V(data, W1V, W3V))

    beta1V <- tr(U, W1V) / p
    beta3Vbar <- mean(mapply(function(Ki, W3i) tr(Ki, W3i), LK, W3V)) / Nbar

    alpha1V <- tr(W1V) / p
    alpha3Vbar <- mean(sapply(W3V, tr)) / Nbar

    sig2_cap <- min(eigen(V_hat, symmetric=TRUE, only.values=TRUE)$values)
    sig2_use <- max(0, min(sig2, 0.99 * sig2_cap / (alpha1V * alpha3Vbar + 1e-12)))
    V_tilde <- V_hat - diag(sig2_use * alpha1V * alpha3Vbar, q)
    V_tilde <- psd_project(V_tilde)

    V <- V_tilde / (beta1V * beta3Vbar)
    V <- q * V / tr(V)
  }

  list(
    C = create_matrix_covfunc(U, V, K, t_grid, sig2),
    W = if (wScheme == "identity") NULL else W
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

  return(scores)
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
  return(X)
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
