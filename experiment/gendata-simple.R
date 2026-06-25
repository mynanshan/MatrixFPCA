library(fda)
library(rTensor)

source("matfdata.R")

gen_data_simple <- function(
  n = 10, p = 10, q = 10,
  npc_row = p, npc_col = q, npc_time = 3,
  nt = 51, t_range = c(0,1),
  sampling_ratio = NULL,
  sig = 0,
  eigvecSet = 1,
  eigvalSet = 1,
  scores_distr = c("normal", "uniform", "beta32", "t5"),
  score_dependence = c("independent", "common_amplitude"),
  amplitude_var = 1
) {
  
  t_grid = seq(t_range[1], t_range[2], length.out=nt)
  
  # dimension 1: row
  # dimension 2: column
  # dimension 3: time
  
  # eigenfunction
  # NOTE: when npc < corresponding dimension
  # the covariance matrix is singular
  # then following procedure will be incorrect
  npc = npc_col * npc_row * npc_time
  if (eigvecSet == 1) {
    Phi = pracma::randortho(p)[,1:npc_row,drop=FALSE]
    Psi = pracma::randortho(q)[,1:npc_col,drop=FALSE]
  } else if (eigvecSet == 2) {
    Phi = fda::fourier(
      seq(0,1,length.out=p), nbasis = npc_row+1
    )[,2:(npc_row+1)]
    Psi = fda::fourier(
      seq(0,1,length.out=q), nbasis = npc_col+1
    )[,2:(npc_col+1)]
    Phi = qr.Q(qr(Phi))
    Psi = qr.Q(qr(Psi))
  }
  theta_funcs = function(x, npc) fda::fourier(x, nbasis = npc)
  theta_mat = theta_funcs(t_grid, npc_time)
  
  # eigenvalues
  if (eigvalSet == 1) {
    lam_phi = seq(2, 0.2, length.out=npc_row)
    lam_psi = seq(2, 0.2, length.out=npc_col)
  } else if (eigvalSet == 2) {
    lam_phi = 2 * (1:npc_row)^(-1.2)
    lam_psi = 2 * (1:npc_col)^(-1.3)
  }
  lam_theta = 2 * (1:npc_time)^(-2.0)

  scale_U <- sum(lam_phi) / p
  scale_V <- sum(lam_psi) / q
  lam_phi <- lam_phi / scale_U
  lam_psi <- lam_psi / scale_V
  lam_theta <- lam_theta * scale_U * scale_V
  
  eig_true = create_matrix_eigobj(Phi, Psi, theta_mat, lam_phi, lam_psi, lam_theta, t_grid)
  
  U = Phi %*% diag(lam_phi,npc_row,npc_row) %*% t(Phi)
  V = Psi %*% diag(lam_psi,npc_col,npc_col) %*% t(Psi)
  K = theta_mat %*% diag(lam_theta,npc_time,npc_time) %*% t(theta_mat)
  
  cov_true = create_matrix_covfunc(U, V, K, t_grid, sig^2)
  
  # FPC scores
  scores_distr <- match.arg(scores_distr)
  score_dependence <- match.arg(score_dependence)
  Z_base <- if (scores_distr == "normal") {
    rnorm(n * npc)
  } else if (scores_distr == "uniform") {
    runif(n * npc, -1, 1) / sqrt(1 / 3)
  } else if (scores_distr == "beta32") {
    (rbeta(n * npc, 3, 2) - 0.6) * 5
  } else if (scores_distr == "t5") {
    rt(n * npc, 5) / sqrt(5 / 3)
  }

  amplitude_var <- max(as.numeric(amplitude_var)[1L], 0)
  A <- rep(1, n)
  if (score_dependence == "common_amplitude") {
    if (amplitude_var <= 0) {
      warning("amplitude_var is zero; common_amplitude reduces to independence.")
    } else {
      # A_i^2 has mean one and variance amplitude_var. Multiplication by A_i
      # preserves pairwise uncorrelated scores but makes their squares dependent.
      shape <- 1 / amplitude_var
      A <- sqrt(rgamma(n, shape = shape, scale = 1 / shape))
    }
  } else {
    amplitude_var <- 0
  }

  Z <- array(
    Z_base,
    dim = c(n, npc_row, npc_col, npc_time)
  )
  Z <- sweep(Z, 1L, A, "*") |> as.tensor()
  scores = ttl(
    Z, list(
      diag(sqrt(lam_phi), npc_row, npc_row),
      diag(sqrt(lam_psi), npc_col, npc_col), 
      diag(sqrt(lam_theta), npc_time, npc_time)
    ), c(2:4)
  )
  
  # tensor X
  X = ttl(scores, list(Phi, Psi, theta_mat), c(2:4))
  
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

    X_list <- lapply(1:n, function(i) X[i,,,t_idx_list[[i]]])
    if (sig > 0) {
      for (i in seq_along(X_list)) {
        X_list[[i]] <- X_list[[i]] + rnorm(prod(dim(X_list[[i]])), 0, sig)
      }
    }
    
    data = create_matrix_funcdata_irreg(
      X_list = X_list,
      t_grid = t_grid,
      t_idx_list = t_idx_list
    )
  }
  
  list(
    data = data,
    cov = cov_true,
    eig = eig_true,
    scores = scores,
    score_dependence = score_dependence,
    amplitude_var = amplitude_var,
    subject_amplitude = A
  )
  
}
