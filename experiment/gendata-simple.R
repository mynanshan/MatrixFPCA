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
  spatial_model = c("ar1", "spectral"),
  scores_distr = c("normal", "uniform", "beta32", "t5"),
  score_dependence = c(
    "independent",
    "common_amplitude",
    "common_signed",
    "heavy_amplitude"
  ),
  dependence_param = 1
) {
  spatial_model <- match.arg(spatial_model)
  
  t_grid = seq(t_range[1], t_range[2], length.out=nt)
  
  # dimension 1: row
  # dimension 2: column
  # dimension 3: time
  
  # eigenfunction
  # NOTE: when npc < corresponding dimension
  # the covariance matrix is singular
  # then following procedure will be incorrect
  npc = npc_col * npc_row * npc_time
  if (spatial_model == "ar1") {
    ar1_eigensystem <- function(d, rank, rho = 0.7) {
      covariance <- rho^abs(outer(seq_len(d), seq_len(d), "-"))
      ee <- eigen(covariance, symmetric = TRUE)
      list(
        vectors = ee$vectors[, seq_len(rank), drop = FALSE],
        values = ee$values[seq_len(rank)]
      )
    }
    row_eig <- ar1_eigensystem(p, npc_row)
    col_eig <- ar1_eigensystem(q, npc_col)
    Phi <- row_eig$vectors
    Psi <- col_eig$vectors
    lam_phi <- row_eig$values
    lam_psi <- col_eig$values
  } else {
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
    } else {
      stop("eigvecSet must be 1 or 2.")
    }
  }
  theta_funcs = function(x, npc) fda::fourier(x, nbasis = npc)
  theta_mat = theta_funcs(t_grid, npc_time)
  
  # eigenvalues
  if (spatial_model == "spectral") {
    if (eigvalSet == 1) {
      lam_phi = seq(2, 0.2, length.out=npc_row)
      lam_psi = seq(2, 0.2, length.out=npc_col)
    } else if (eigvalSet == 2) {
      lam_phi = 2 * (1:npc_row)^(-1.2)
      lam_psi = 2 * (1:npc_col)^(-1.3)
    } else {
      stop("eigvalSet must be 1 or 2.")
    }
  }
  if (eigvalSet == 1) {
    lam_theta = 2 * (1:npc_time)^(-2.0)
  } else if (eigvalSet == 2) {
    lam_theta = 2 * (1:npc_time)^(-2.5)
  } else {
    stop("eigvalSet must be 1 or 2.")
  }

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

  draw_standardized <- function(size) {
    if (scores_distr == "normal") {
      rnorm(size)
    } else if (scores_distr == "uniform") {
      runif(size, -1, 1) / sqrt(1 / 3)
    } else if (scores_distr == "beta32") {
      (rbeta(size, 3, 2) - 0.6) * 5
    } else if (scores_distr == "t5") {
      rt(size, 5) / sqrt(5 / 3)
    }
  }

  dependence_param <- as.numeric(dependence_param)[1L]
  if (!is.finite(dependence_param)) {
    stop("dependence_param must be finite.")
  }

  A <- rep(1, n)
  Z_base <- draw_standardized(n * npc)

  if (score_dependence == "common_amplitude") {
    # Z_{i,alpha} = A_i xi_{i,alpha}, E(A_i^2) = 1. Scores remain
    # pairwise uncorrelated, but all squared scores share A_i^2.
    if (dependence_param <= 0) {
      warning("dependence_param is non-positive; common_amplitude reduces to independence.")
      dependence_param <- 0
    } else {
      shape <- 1 / dependence_param
      A <- sqrt(rgamma(n, shape = shape, scale = 1 / shape))
    }
  } else if (score_dependence == "common_signed") {
    # The user's common-variable construction:
    # Z_{i,alpha} = A_i S_{i,alpha}, with independent Rademacher signs.
    # Every score has unit variance and distinct scores are uncorrelated, while
    # Z_{i,alpha}^2 = A_i^2 for every alpha, giving maximal squared dependence.
    A <- abs(draw_standardized(n))
    Z_base <- sample(c(-1, 1), n * npc, replace = TRUE)
    dependence_param <- 0
  } else if (score_dependence == "heavy_amplitude") {
    # A_i^2 is normalized Pareto with tail index dependence_param. For
    # 1 < dependence_param <= 2, E(A_i^2) exists but E(A_i^4) is infinite,
    # deliberately violating the fourth-moment assumptions.
    tail_index <- dependence_param
    if (tail_index <= 1) {
      stop("heavy_amplitude requires dependence_param > 1.")
    }
    pareto <- runif(n)^(-1 / tail_index)
    A <- sqrt(pareto * (tail_index - 1) / tail_index)
  } else {
    dependence_param <- 0
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
    visit_counts = rep(nt, n)
  } else {
    # Sparse irregular sampling. The ratios specify the approximate fraction
    # of the common nt-point visit grid retained for each subject; every pixel
    # in an image is observed at each retained visit.
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
    visit_counts = lengths(t_idx_list)
  }
  
  list(
    data = data,
    cov = cov_true,
    eig = eig_true,
    scores = scores,
    score_dependence = score_dependence,
    dependence_param = dependence_param,
    spatial_model = spatial_model,
    sampling_ratio = sampling_ratio,
    visit_counts = visit_counts,
    subject_amplitude = A
  )
  
}
