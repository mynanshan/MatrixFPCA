library(parallel)
library(foreach)
library(doFuture)
library(rTensor)

source("matfdata.R")

# TODO: add support for mean functions

# Cov Estimation ----------------------------------------------------------

### Formulate tilde statistics

get_working_K <- function(X) {
  # \frac{1}{tr(U)tr(V)} unfold_3(X) unfold_3(X)^T
  # Input:
  #   X: one observation tensor, (p, q, mi)
  #   U, V: row/col cov., (p, p) and (q, q)
  # Return:
  #   working K, mi x mi matrix
  res <- unfold(X, row_idx = 3, col_idx = c(1, 2))@data %*%
    t(unfold(X, row_idx = 3, col_idx = c(1, 2))@data)
  return(res)
}

get_working_V <- function(X) {
  # \frac{1}{tr(U)tr(K)} unfold_2(X) unfold_2(X)^T
  # Input:
  #   X: one observation tensor, (p, q, m)
  #   U, K: row/temporal cov, (p, p) and (m, m)
  # Return:
  #   working V, q x q matrix
  res <- unfold(X, row_idx = 2, col_idx = c(1, 3))@data %*%
    t(unfold(X, row_idx = 2, col_idx = c(1, 3))@data)
  return(res)
}

get_working_U <- function(X) {
  # \frac{1}{tr(V)tr(K)} unfold_1(X) unfold_1(X)^T
  # Input:
  #   X: one observation tensor, (p, q, m)
  #   Vinv, Kinv: inverse of col/temporal cov, (q, q) and (m, m)
  # Return:
  #   working V, p x p matrix
  res <- unfold(X, row_idx = 1, col_idx = c(2, 3))@data %*%
    t(unfold(X, row_idx = 1, col_idx = c(2, 3))@data)
  return(res)
}


### Update one of U,V,K with others fixed

estimate_K <- function(
    data, U, V,
    wScheme = c("OBS", "SUBJ"),
    use_parallel = FALSE, ...) {
  UseMethod("estimate_K")
}

estimate_K.matfdata.grid <- function(
    data,
    wScheme = c("OBS", "SUBJ"),
    use_parallel = FALSE) {
  if (!is(data, "matfdata.grid")) stop("data needs to be a matfdata.grid")
  # wScheme is not used in this case

  n <- dim(data$X)[1]
  p <- data$vdim[1]
  q <- data$vdim[2]

  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  K_hat <- foreach(i = seq_len(n), .combine = "+") %op% {
    # get the working K
    Ki_wk <- get_working_K(get_subj_tnsr(data, i))
    # TODO: check weighting
    Ki_wk / (n * p * q)
  }

  return(K_hat)
}

estimate_K.matfdata.irreg <- function(
    data,
    wScheme = c("OBS", "SUBJ"),
    use_parallel = FALSE) {
  if (!is(data, "matfdata.irreg")) stop("data needs to be a matfdata.irreg")
  wScheme <- match.arg(wScheme)

  n <- length(data$X_list)
  nt <- length(data$t_grid)
  p <- data$vdim[1]
  q <- data$vdim[2]

  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  K_wk_all <- foreach::foreach(
    i = seq_len(n),
    .combine = rbind
  ) %op% {
    # get the working K
    Ki_wk <- get_working_K(get_subj_tnsr(data, i))

    # only keep the lower half due to symmetry
    low_tri_id <- lower.tri(Ki_wk, diag = TRUE)
    Ki_wk <- Ki_wk[low_tri_id]

    # create a meshgrid of (t1, t2)
    t_cov_mesh <- pracma::meshgrid(data$t_grid[data$t_idx_list[[i]]])

    # only keep the corresponding lower half (t2 > t1)
    t_cov_grid <- cbind(t_cov_mesh$X[low_tri_id], t_cov_mesh$Y[low_tri_id])

    # points weighting
    if (wScheme == "OBS") {
      w <- rep(1, length(Ki_wk))
    } else {
      w <- rep(1 / length(Ki_wk), length(Ki_wk))
    }

    data.table::data.table(
      t1 = t_cov_grid[, 1],
      t2 = t_cov_grid[, 2],
      Kval = as.vector(Ki_wk) / (p * q), # divide by tr(U)tr(V)
      w = w
    )
  }

  # fit the smooth surface of K
  # TODO: optimize using mgcv::bam
  K_wk_all$w <- K_wk_all$w / mean(K_wk_all$w)
  fit_K <- mgcv::gam(
    Kval ~ ti(t1) + ti(t2) + ti(t1, t2),
    data = K_wk_all,
    weights = K_wk_all$w
  )

  # calculate K(s,t) on a dense mesh
  K_u_hat <- matrix(0, nt, nt)
  low_tri_id <- lower.tri(K_u_hat, diag = TRUE)
  t_cov_mesh <- pracma::meshgrid(data$t_grid)
  t_cov_grid <- cbind(t_cov_mesh$X[low_tri_id], t_cov_mesh$Y[low_tri_id])
  K_u_hat[low_tri_id] <- predict(
    fit_K,
    newdata = data.table::data.table(
      t1 = t_cov_grid[, 1], t2 = t_cov_grid[, 2]
    )
  )
  # copy the lower triangle to upper triangle
  K_hat <- K_u_hat + t(K_u_hat)
  diag(K_hat) <- diag(K_u_hat)

  return(K_hat)
}

estimate_U <- function(
    data,
    wScheme = c("OBS", "SUBJ"),
    use_parallel = FALSE, ...) {
  UseMethod("estimate_U")
}

estimate_U.matfdata.grid <- function(
    data,
    wScheme = c("OBS", "SUBJ"),
    use_parallel = FALSE) {
  if (!is(data, "matfdata.grid")) stop("data needs to be a matfdata.grid")
  # wSceheme is not used in this case

  n <- dim(data$X)[1]
  nt <- length(data$t_grid)
  p <- data$vdim[1]
  q <- data$vdim[2]

  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  U_hat <- foreach(i = seq_len(n), .combine = "+") %op% {
    U_wk <- get_working_U(get_subj_tnsr(data, i))
    U_wk / (n * q * nt)
  }

  return(U_hat)
}

# FIXME: recheck the weighting schemes

estimate_U.matfdata.irreg <- function(
    data,
    wScheme = c("OBS", "SUBJ"),
    use_parallel = FALSE) {
  if (!is(data, "matfdata.irreg")) stop("data needs to be a matfdata.irreg")
  wScheme <- match.arg(wScheme)

  n <- length(data$X_list)
  nt <- length(data$t_grid)
  p <- data$vdim[1]
  q <- data$vdim[2]
  totalN <- sum(sapply(data$t_idx_list, length))

  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  # average all working U
  U_hat <- foreach::foreach(i = seq_len(n), .combine = "+") %op% {
    t_id <- data$t_idx_list[[i]]
    Ni <- length(t_id)
    U_wk <- get_working_U(get_subj_tnsr(data, i))
    if (wScheme == "OBS") {
      return(U_wk / q)
    } else {
      return(U_wk / (Ni * q))
    }
  }

  if (wScheme == "OBS") {
    return(U_hat / totalN)
  } else {
    return(U_hat / n)
  }

  return(U_hat)
}

estimate_V <- function(
    data,
    wScheme = c("OBS", "SUBJ"),
    use_parallel = FALSE, ...) {
  UseMethod("estimate_V")
}

estimate_V.matfdata.grid <- function(
    data,
    wScheme = c("OBS", "SUBJ"),
    use_parallel = FALSE) {
  if (!is(data, "matfdata.grid")) stop("data needs to be a matfdata.grid")
  wScheme <- match.arg(wScheme)

  n <- dim(data$X)[1]
  nt <- length(data$t_grid)
  p <- data$vdim[1]
  q <- data$vdim[2]

  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  V_hat <- foreach(i = seq_len(n), .combine = "+") %op% {
    V_wk <- get_working_V(get_subj_tnsr(data, i))
    V_wk / (n * p * nt)
  }

  return(V_hat)
}

estimate_V.matfdata.irreg <- function(
    data, U, K,
    wScheme = c("OBS", "SUBJ"),
    use_parallel = FALSE) {
  if (!is(data, "matfdata.irreg")) stop("data needs to be a matfdata.irreg")
  wScheme <- match.arg(wScheme)

  n <- length(data$X_list)
  nt <- length(data$t_grid)
  p <- data$vdim[1]
  q <- data$vdim[2]
  totalN <- sum(sapply(data$t_idx_list, length))

  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  # average all working V
  V_hat <- foreach::foreach(i = seq_len(n), .combine = "+") %op% {
    t_id <- data$t_idx_list[[i]]
    Ni <- length(t_id)
    V_wk <- get_working_V(get_subj_tnsr(data, i))
    if (wScheme == "OBS") {
      return(V_wk / p)
    } else {
      return(V_wk / (Ni * p))
    }
  }

  if (wScheme == "OBS") {
    return(V_hat / totalN)
  } else {
    return(V_hat / n)
  }

  return(V_hat)
}

### Main Estimation Function

fit_matfcov <- function(
    data,
    wScheme = c("OBS", "SUBJ"),
    use_parallel = FALSE) {
  wScheme <- match.arg(wScheme)

  p <- data$vdim[1]
  q <- data$vdim[2]
  t_grid <- data$t_grid

  # Iterations
  U <- estimate_U(data, wScheme, use_parallel)
  V <- estimate_V(data, wScheme, use_parallel)
  K <- estimate_K(data, wScheme, use_parallel)
  U <- p * U / sum(diag(U))
  V <- q * V / sum(diag(V))

  return(create_matrix_covfunc(U, V, K, t_grid))
}


get_matfpca_scores <- function(
    data, eig, nugget = 1e-8, use_parallel = FALSE) {
  n <- count_subj(data)
  I <- ncol(eig$funcs$Phi)
  J <- ncol(eig$funcs$Psi)
  K <- ncol(eig$funcs$Rho)

  scores <- rTensor::as.tensor(array(NA, dim = c(n, I, J, K)))
  lam_rho <- eig$values$lam_rho

  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  # for (i in seq_len(n)) {
  scores_list <- foreach::foreach(
    i = seq_len(n), .combine = c, .multicombine = TRUE
  ) %op% {
    Rhoi <- eig$funcs$Rho[data$t_idx_list[[i]], , drop = F]
    mi <- length(data$t_idx_list[[i]])
    Ki_inv <- solve(Rhoi %*% t(sweep(Rhoi, 2, lam_rho, "*")) + diag(nugget, mi, mi))
    scores_i <- rTensor::ttl(
      get_subj_tnsr(data, i),
      list(
        t(eig$funcs$Phi),
        t(eig$funcs$Psi),
        t(sweep(Rhoi, 2, lam_rho, "*")) %*% Ki_inv
      ), 1:3
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
  K <- ncol(eigf$Rho)
  stopifnot(length(dim(scores)) == 4)
  stopifnot(all(dim(scores)[2:4] == c(I, J, K)))
  if (is.array(scores)) scores <- rTensor::as.tensor(scores)

  X <- with(eigf, ttl(scores, list(Phi, Psi, Rho), 2:4))
  return(X)
}


# Other Helpers -----------------------------------------------------------

scale_mat_trace <- function(X) {
  stopifnot(is.matrix(X))
  stopifnot(nrow(X) == ncol(X))
  d <- nrow(X)
  trace_X <- sum(diag(X))
  scale_X <- trace_X / d
  X <- X / scale_X
  attr(X, "scale") <- scale_X
  return(X)
}
