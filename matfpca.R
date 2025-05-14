library(parallel)
library(foreach)
library(doFuture)

source("matfdata.R")

# TODO: add support for mean functions

# Cov Estimation ----------------------------------------------------------

### Formulate tilde statistics

get_working_K.mle = function(X, Uinv, Vinv) {
  # \frac{1}{pq} unfold_3(X) (V^{-1} \otimes U^{-1}) unfold_3(X)^T
  # Input:
  #   X: one observation tensor, (p, q, mi)
  #   Uinv, Vinv: inverse of row/col cov, (p, p) and (q, q)
  # Return:
  #   working K, mi x mi matrix
  p = dim(X)[1]; q = dim(X)[2]
  X_decor_rc = rTensor::ttl(X, list(Uinv, Vinv), c(1,2))
  return(
    unfold(X_decor_rc, row_idx = 3, col_idx = c(1,2))@data %*%
    t(unfold(X, row_idx = 3, col_idx = c(1,2))@data) / (p * q)
  )
}

get_working_K.basic = function(X, U, V) {
  # \frac{1}{tr(U)tr(V)} unfold_3(X) unfold_3(X)^T
  # Input:
  #   X: one observation tensor, (p, q, mi)
  #   U, V: row/col cov., (p, p) and (q, q)
  # Return:
  #   working K, mi x mi matrix
  p = dim(X)[1]; q = dim(X)[2]
  scaling = sum(diag(U)) * sum(diag(V))
  return(
    unfold(X, row_idx = 3, col_idx = c(1,2))@data %*%
    t(unfold(X, row_idx = 3, col_idx = c(1,2))@data) / scaling
  )
}

get_working_V.mle = function(X, Uinv, Kinv) {
  # \frac{1}{pm} unfold_2(X) (K^{-1} \otimes U^{-1}) unfold_2(X)^T
  # Input:
  #   X: one observation tensor, (p, q, m)
  #   Uinv, Kinv: inverse of row/temporal cov, (p, p) and (m, m)
  # Return:
  #   working V, q x q matrix
  p = dim(X)[1]; m = dim(X)[3]
  X_decor_rt = rTensor::ttl(X, list(Uinv, Kinv), c(1,3))
  return(
    unfold(X_decor_rt, row_idx = 2, col_idx = c(1,3))@data %*%
    t(unfold(X, row_idx = 2, col_idx = c(1,3))@data) / (p * m)
  )
}

get_working_V.basic = function(X, U, K) {
  # \frac{1}{tr(U)tr(K)} unfold_2(X) unfold_2(X)^T
  # Input:
  #   X: one observation tensor, (p, q, m)
  #   U, K: row/temporal cov, (p, p) and (m, m)
  # Return:
  #   working V, q x q matrix
  p = dim(X)[1]; m = dim(X)[3]
  scaling = sum(diag(U)) * sum(diag(K))
  return(
    unfold(X, row_idx = 2, col_idx = c(1,3))@data %*%
    t(unfold(X, row_idx = 2, col_idx = c(1,3))@data) / scaling
  )
}

get_working_U.mle = function(X, Vinv, Kinv) {
  # \frac{1}{qm} unfold_1(X) (K^{-1} \otimes V^{-1}) unfold_1(X)^T
  # Input:
  #   X: one observation tensor, (p, q, m)
  #   Vinv, Kinv: inverse of col/temporal cov, (q, q) and (m, m)
  # Return:
  #   working V, p x p matrix
  q = dim(X)[2]; m = dim(X)[3]
  X_decor_ct = rTensor::ttl(X, list(Vinv, Kinv), c(2,3))
  return(
    unfold(X_decor_ct, row_idx = 1, col_idx = c(2,3))@data %*%
    t(unfold(X, row_idx = 1, col_idx = c(2,3))@data) / (q * m)
  )
}

get_working_U.basic = function(X, V, K) {
  # \frac{1}{tr(V)tr(K)} unfold_1(X) unfold_1(X)^T
  # Input:
  #   X: one observation tensor, (p, q, m)
  #   Vinv, Kinv: inverse of col/temporal cov, (q, q) and (m, m)
  # Return:
  #   working V, p x p matrix
  q = dim(X)[2]; m = dim(X)[3]
  scaling = sum(diag(V)) * sum(diag(K))
  return(
    unfold(X, row_idx = 1, col_idx = c(2,3))@data %*%
    t(unfold(X, row_idx = 1, col_idx = c(2,3))@data) / scaling
  )
}


### Update one of U,V,K with others fixed

estimate_K = function(
  data, U, V, nugget = 0, moment_type = c("mle", "basic"),
  use_parallel = FALSE, ...) UseMethod("estimate_K")

estimate_K.matfdata.grid = function(
  data, U, V, nugget = 0, moment_type = c("mle", "basic"),
  use_parallel = FALSE
) {
  if (!is(data, "matfdata.grid")) stop("data needs to be a matfdata.grid")
  moment_type = match.arg(moment_type)
  
  n = dim(data$X)[1]
  nt = length(data$t_grid)
  p = data$vdim[1]
  q = data$vdim[2]
  
  if (moment_type == "mle") {
    Uinv = solve(U + diag(nugget, p, p))
    Vinv = solve(V + diag(nugget, q, q))
  }

  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  K_hat = foreach(i = seq_len(n), .combine = "+") %op% {
    # get the working K
    if (moment_type == "mle") {
      Ki_wk = get_working_K.mle(get_subj_tnsr(data, i), Uinv, Vinv)
    } else if (moment_type == "basic") {
      Ki_wk = get_working_K.basic(get_subj_tnsr(data, i), U, V)
    } else {stop("Unknown moment_type.")}
    Ki_wk / n
  }
  
  return(K_hat)
}

estimate_K.matfdata.irreg = function(
  data, U, V, nugget = 0, moment_type = c("mle", "basic"),
  use_parallel = FALSE
) {
  if (!is(data, "matfdata.irreg")) stop("data needs to be a matfdata.irreg")
  moment_type = match.arg(moment_type)
  
  n = length(data$X_list)
  nt = length(data$t_grid)
  p = data$vdim[1]
  q = data$vdim[2]
  
  if (moment_type == "mle") {
    Uinv = solve(U + diag(nugget, p, p))
    Vinv = solve(V + diag(nugget, q, q))
  }
  
  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  K_wk_all <- foreach::foreach(
    i = seq_len(n),
    .combine = rbind
  ) %op% {
    # get the working K
    if (moment_type == "mle") {
      Ki_wk = get_working_K.mle(get_subj_tnsr(data, i), Uinv, Vinv)
    }
    if (moment_type == "basic") {
      Ki_wk = get_working_K.basic(get_subj_tnsr(data, i), U, V)
    }
    
    # only keep the lower half due to symmetry
    low_tri_id = lower.tri(Ki_wk, diag = TRUE)
    Ki_wk = Ki_wk[low_tri_id]
    
    # create a meshgrid of (t1, t2)
    t_cov_mesh = pracma::meshgrid(data$t_grid[data$t_idx_list[[i]]])
    
    # only keep the corresponding lower half (t2 > t1)
    t_cov_grid = cbind(t_cov_mesh$X[low_tri_id], t_cov_mesh$Y[low_tri_id])
    
    data.table::data.table(
      t1 = t_cov_grid[,1],
      t2 = t_cov_grid[,2],
      Kval = as.vector(Ki_wk))
  }
  # TODO: should we standardize K_wk first and directly use mgcv::gam
  # or keep the unstandardized K_wk but assign weights in mgcv::gam?
  # browser()
  
  # fit the smooth surface of K
  fit_K = mgcv::gam(Kval ~ ti(t1) + ti(t2) + ti(t1, t2), data = K_wk_all)
  
  # calculate K(s,t) on a dense mesh
  K_u_hat = matrix(0, nt, nt)
  low_tri_id = lower.tri(K_u_hat, diag = TRUE)
  t_cov_mesh = pracma::meshgrid(data$t_grid)
  t_cov_grid = cbind(t_cov_mesh$X[low_tri_id], t_cov_mesh$Y[low_tri_id])
  K_u_hat[low_tri_id] = predict(
    fit_K, newdata = data.table::data.table(
    t1 = t_cov_grid[,1], t2 = t_cov_grid[,2]))
  # copy the lower triangle to upper triangle
  K_hat = K_u_hat + t(K_u_hat)
  diag(K_hat) = diag(K_u_hat)
  
  return(K_hat)
}

estimate_U = function(
  data, V, K, nugget = 1e-8, moment_type = c("mle", "basic"),
  use_parallel = FALSE, ...) UseMethod("estimate_U")

estimate_U.matfdata.grid = function(
  data, V, K, nugget = 1e-8, moment_type = c("mle", "basic"),
  use_parallel = FALSE
) {
  if (!is(data, "matfdata.grid")) stop("data needs to be a matfdata.grid")
  moment_type = match.arg(moment_type)
  
  n = dim(data$X)[1]
  nt = length(data$t_grid)
  m = nt
  p = data$vdim[1]
  q = data$vdim[2]
  
  if (moment_type == "mle") {
    Vinv = solve(V + diag(nugget, q, q))
    Kinv = solve(K + diag(nugget, nt, nt))
  }
  if (moment_type == "basic") {
    trace_V = sum(diag(V))
    trace_K = sum(diag(K))
  }
  
  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  U_hat = foreach(i = seq_len(n), .combine = "+") %op% {
    if (moment_type == "mle") {
      U_wk = get_working_U.mle(get_subj_tnsr(data, i), Vinv, Kinv)
    } else if (moment_type == "basic") {
      U_wk = get_working_U.basic(get_subj_tnsr(data, i), V, K)
    }
    U_wk / n
  }
  
  return(U_hat)
}

estimate_U.matfdata.irreg = function(
  data, V, K, nugget = 1e-8, moment_type = c("mle", "basic"),
  use_parallel = FALSE
) {
  if (!is(data, "matfdata.irreg")) stop("data needs to be a matfdata.irreg")
  moment_type = match.arg(moment_type)
  
  n = length(data$X_list)
  nt = length(data$t_grid)
  p = data$vdim[1]
  q = data$vdim[2]
  
  if (moment_type == "mle") {
    Vinv = solve(V + diag(nugget, q, q))
  } else if (moment_type == "basic") {
    trace_V = sum(diag(V))
  }
  
  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  # weights for working U
  wv = switch(moment_type, "mle" = q, "basic" = trace_V)
  wk = sapply(seq_len(n), \(i) {
    t_id = data$t_idx_list[[i]]
    mi = length(t_id)
    switch(moment_type, "mle" = mi, "basic" = sum(diag(K)[t_id]))
  })
  w = wv * wk

  # average all working U
  U_hat = foreach::foreach(i = seq_len(n), .combine = "+") %op% {
    t_id = data$t_idx_list[[i]]
    mi = length(t_id)
    if (moment_type == "mle") {
      Ki_inv = solve(K[t_id, t_id] + diag(nugget, mi, mi))
      U_wk = get_working_U.mle(get_subj_tnsr(data, i), Vinv, Ki_inv)
    } else if (moment_type == "basic") {
      Ki = K[t_id, t_id]
      U_wk = get_working_U.basic(get_subj_tnsr(data, i), V, Ki)
    }
    w[i] * U_wk / n
  }
  U_hat = U_hat / mean(w)
  
  return(U_hat)
}

estimate_V = function(
  data, U, K, nugget = 1e-8, moment_type = c("mle", "basic"),
  use_parallel = FALSE, ...) UseMethod("estimate_V")

estimate_V.matfdata.grid = function(
  data, U, K, nugget = 1e-8, moment_type = c("mle", "basic"),
  use_parallel = FALSE
) {
  if (!is(data, "matfdata.grid")) stop("data needs to be a matfdata.grid")
  moment_type = match.arg(moment_type)
  
  n = dim(data$X)[1]
  nt = length(data$t_grid)
  p = data$vdim[1]
  q = data$vdim[2]
  
  if (moment_type == "mle") {
    Uinv = solve(U + diag(nugget, p, p))
    Kinv = solve(K + diag(nugget, nt, nt))
  } else if (moment_type == "basic") {
    trace_U = sum(diag(U))
    trace_K = sum(diag(K))
  }
  
  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  V_hat = foreach(i = seq_len(n), .combine = "+") %op% {
    if (moment_type == "mle") {
      V_wk = get_working_V.mle(get_subj_tnsr(data, i), Uinv, Kinv)
    } else if (moment_type == "basic") {
      V_wk = get_working_V.basic(get_subj_tnsr(data, i), U, K)
    }
    V_wk / n
  }
  
  return(V_hat)
}

estimate_V.matfdata.irreg = function(
  data, U, K, nugget = 1e-8, moment_type = c("mle", "basic"),
  use_parallel = FALSE
) {
  if (!is(data, "matfdata.irreg")) stop("data needs to be a matfdata.irreg")
  moment_type = match.arg(moment_type)
  
  n = length(data$X_list)
  nt = length(data$t_grid)
  p = data$vdim[1]
  q = data$vdim[2]
  
  if (moment_type == "mle") {
    Uinv = solve(U + diag(nugget, p, p))
  } else if (moment_type == "basic") {
    trace_U = sum(diag(U))
  }
  
  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`

  # weights for working V
  wu = switch(moment_type, "mle" = p, "basic" = trace_U)
  wk = sapply(seq_len(n), \(i) {
    t_id = data$t_idx_list[[i]]
    mi = length(t_id)
    switch(moment_type, "mle" = mi, "basic" = sum(diag(K)[t_id]))
  })
  w = wu * wk

  # average all working V
  V_hat = foreach::foreach(i = seq_len(n), .combine = "+") %op% {
    t_id = data$t_idx_list[[i]]
    mi = length(t_id)
    if (moment_type == "mle") {
      Ki_inv = solve(K[t_id, t_id] + diag(nugget, mi, mi))
      V_wk = get_working_V.mle(get_subj_tnsr(data, i), Uinv, Ki_inv)
    } else if (moment_type == "basic") {
      Ki = K[t_id, t_id]
      V_wk = get_working_V.basic(get_subj_tnsr(data, i), U, Ki)
    }
    w[i] * V_wk / n
  }
  V_hat = V_hat / mean(w)
  
  return(V_hat)
}

### Main Estimation Function

set_cov_inits = function(p, q, t_grid, U=NULL, V=NULL, K=NULL) {
  if (is.null(U)) {
    U = diag(1, p, p)
  }
  if (is.null(V)) {
    V = diag(1, q, q)
  }
  if (is.null(K)) {
    t_span = diff(range(t_grid))
    K = exp(-2 * (outer(t_grid, t_grid, "-") / t_span)^2)
  }
  cov_inits = create_matrix_covfunc(U, V, K, t_grid)
  stopifnot(cov_inits$vdim[1] == p)
  stopifnot(cov_inits$vdim[2] == q)
  return(cov_inits)
}


fit_matfcov = function(
  data, U = NULL, V = NULL, K = NULL,
  max_iter = 20, tol = 1e-4, nugget = rep(1e-8, 3),
  moment_type = c("mle", "basic"),
  use_parallel = FALSE
) {
  moment_type = match.arg(moment_type)
  
  p = data$vdim[1]
  q = data$vdim[2]
  t_grid = data$t_grid
  
  # Initialization
  cov_inits = set_cov_inits(p, q, t_grid, U, V, K)
  U = cov_inits$U
  V = cov_inits$V
  K = cov_inits$K
  U_old = U
  V_old = V
  K_old = K
  
  # Iterations
  for (it in seq_len(max_iter)) {
    K = estimate_K(data, U, V, nugget[3], moment_type, use_parallel)
    U = estimate_U(data, V, K, nugget[1], moment_type, use_parallel)
    V = estimate_V(data, U, K, nugget[2], moment_type, use_parallel)
    U_factor = sum(diag(U)) / p
    U = U / U_factor
    K = K * U_factor
    V_factor = sum(diag(V)) / q
    V = V / V_factor
    K = K * V_factor
    
    # check convergence
    is_converged = mean((U - U_old)^2) < tol &&
      mean((V - V_old)^2) < tol &&
      mean((K - K_old)^2) < tol
    if (is_converged) break
    U_old = U
    V_old = V
    K_old = K
  
    # cat("Iter:", it)
    # print(mean((scale_mat_trace(U) - scale_mat_trace(cov_true$U))^2))
    # print(mean((scale_mat_trace(V) - scale_mat_trace(cov_true$V))^2))
    # print(mean((scale_mat_trace(K) - scale_mat_trace(cov_true$K))^2))
  }
  
  return(list(
    cov = create_matrix_covfunc(U, V, K, t_grid),
    is_converged = is_converged, iter = it
  ))
}


get_matfpca_scores = function(
  data, eig, nugget=1e-8, use_parallel = FALSE) {
  
  n = count_subj(data)
  I = ncol(eig$funcs$Phi)
  J = ncol(eig$funcs$Psi)
  K = ncol(eig$funcs$Rho)
  
  scores = rTensor::as.tensor(array(NA, dim=c(n, I, J, K)))
  lam_rho = eig$values$lam_rho
  
  `%op%` <- if (use_parallel) `%dofuture%` else `%do%`
  
  # for (i in seq_len(n)) {
  scores_list = foreach::foreach(
    i = seq_len(n), .combine = c, .multicombine = TRUE
  ) %op% {
    Rhoi = eig$funcs$Rho[data$t_idx_list[[i]],,drop=F]
    mi = length(data$t_idx_list[[i]])
    Ki_inv = solve(Rhoi %*% t(sweep(Rhoi, 2, lam_rho, "*")) + diag(nugget, mi, mi))
    scores_i = rTensor::ttl(
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
    scores[i,,,] = scores_list[[i]]
  }
  
  return(scores)
}

recover_matfdata = function(eigf, scores, mean_func = NULL) {
  
  # validate arguments
  I = ncol(eigf$Phi)
  J = ncol(eigf$Psi)
  K = ncol(eigf$Rho)
  stopifnot(length(dim(scores)) == 4)
  stopifnot(all(dim(scores)[2:4] == c(I, J, K)))
  if (is.array(scores)) scores = rTensor::as.tensor(scores)
  
  X = with(eigf, ttl(scores, list(Phi, Psi, Rho), 2:4))
  return(X)
}


# Other Helpers -----------------------------------------------------------

scale_mat_trace = function(X) {
  stopifnot(is.matrix(X))
  stopifnot(nrow(X) == ncol(X))
  d = nrow(X)
  trace_X = sum(diag(X))
  scale_X = trace_X / d
  X = X / scale_X
  attr(X, "scale") = scale_X
  return(X)
}
