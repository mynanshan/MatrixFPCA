library(fda)
library(rTensor)

source("matfdata.R")

gen_data_simple <- function(
  n = 10, p = 10, q = 10,
  npc_row = p, npc_col = q, npc_time = 3,
  nt = 51, t_range = c(0,1),
  sampling_ratio = NULL,
  eigset = 1
) {
  
  t_grid = seq(t_range[1], t_range[2], length.out=nt)
  
  # dimension 1: row
  # dimension 2: column
  # dimension 3: time
  
  # eigenfunction
  # FIXME: when npc < corresponding dimension
  # the covariance matrix is singular
  # then following procedure will be incorrect
  npc = npc_col * npc_row * npc_time
  Phi = pracma::randortho(p)[,1:npc_row,drop=FALSE]
  Psi = pracma::randortho(q)[,1:npc_col,drop=FALSE]
  rho_funcs = function(x, npc) fda::fourier(x, nbasis = npc)
  # FIXME: nt fourier bases evaluated on the nt-grid gives a constant zero basis.
  rho_mat = rho_funcs(t_grid, npc_time)
  
  # eigenvalues
  if (eigset == 1) {
    lam_phi = seq(2, 0.2, length.out=npc_row)
    lam_psi = seq(2, 0.2, length.out=npc_col)
  } else if (eigset == 2) {
    lam_phi = 2 * (1:npc_row)^(-1.2)
    lam_psi = 2 * (1:npc_col)^(-1.3)
  }
  lam_rho = 2 * (1:npc_time)^(-2.0)
  
  eig_true = create_matrix_eigobj(Phi, Psi, rho_mat, lam_phi, lam_psi, lam_rho, t_grid)
  
  U = Phi %*% diag(lam_phi,npc_row,npc_row) %*% t(Phi)
  V = Psi %*% diag(lam_psi,npc_col,npc_col) %*% t(Psi)
  K = rho_mat %*% diag(lam_rho,npc_time,npc_time) %*% t(rho_mat)
  
  # re-scale covariances
  scale_U <- sum(lam_phi) / p
  scale_V <- sum(lam_psi) / q
  U <- U / scale_U
  V <- V / scale_V
  K <- K * scale_U * scale_V
  lam_phi <- lam_phi / scale_U
  lam_psi <- lam_psi / scale_V
  lam_rho <- lam_rho * scale_U * scale_V
  
  cov_true = create_matrix_covfunc(U, V, K, t_grid)
  
  # FPC scores
  Z = rnorm(n * npc) |>
    array(dim=c(n, npc_col, npc_row, npc_time)) |> 
    as.tensor()
  scores = ttl(
    Z, list(
      diag(sqrt(lam_phi), npc_row, npc_row),
      diag(sqrt(lam_psi), npc_col, npc_col), 
      diag(sqrt(lam_rho), npc_time, npc_time)
    ), c(2:4)
  )
  
  # tensor X
  X = ttl(scores, list(Phi, Psi, rho_mat), c(2:4))
  
  if (is.null(sampling_ratio)) {
    # data on a dense grid
    data = create_matrix_funcdata_grid(X, t_grid)
  } else {
    # sparse irregular sampling, sampling_ratio = (Pmin, Pmax)
    m_min = ceiling(nt * sampling_ratio[1])
    m_max = max(floor(nt * sampling_ratio[2]), m_min)
    ms = runif(n, min = m_min, max = m_max) |> round()
    t_idx_sampler = function(n, m, interval_min = 5) {
      sum_interval_min = interval_min * (m - 1)
      if (n < sum_interval_min + m) stop("No enough points.")
      out = sort(sample(1:(n - sum_interval_min), m, replace = FALSE))
      out = as.integer(out + interval_min * (0:(m-1)))
      return(out)
    }
    # interval_min: minimum interval width
    # make sure sparse samples span at least 10% of the domain
    interval_min = ceiling((nt-1)*0.1 / (m_min-1)) - 1
    t_idx_list = mapply(
      t_idx_sampler, ms, n = nt,
      interval_min = interval_min
    )
    
    data = create_matrix_funcdata_irreg(
      X_list = lapply(1:n, function(i) X[i,,,t_idx_list[[i]]]),
      t_grid = t_grid, t_idx_list = t_idx_list
    )
  }
  
  list(
    data = data,
    cov = cov_true,
    eig = eig_true,
    scores = scores
  )
  
}