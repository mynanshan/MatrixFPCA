library(rTensor)

source("matfpca.R")

t_grid <- c(0, 1)
observed <- as.tensor(array(c(1, 1), dim = c(1, 1, 2)))
data <- create_matrix_funcdata_irreg(
  X_list = list(observed),
  t_grid = t_grid,
  t_idx_list = list(1:2)
)
eig <- create_matrix_eigobj(
  Phi = matrix(1),
  Psi = matrix(1),
  Theta = diag(2),
  lam_phi = 1,
  lam_psi = 1,
  lam_theta = c(1, 1e-8),
  t_grid = t_grid
)
K_grid <- diag(c(1, 1e-8))

score_default <- get_matfpca_scores(data, eig, K_grid = K_grid)
score_auto <- get_matfpca_scores(
  data, eig, K_grid = K_grid,
  inverse_strategy = "auto"
)
score_ginv <- get_matfpca_scores(
  data, eig, K_grid = K_grid,
  inverse_strategy = "ginv",
  ginv_tol = 1e-3
)
score_ridge <- get_matfpca_scores(
  data, eig, K_grid = K_grid,
  inverse_strategy = "ridge"
)
score_less_truncated <- get_matfpca_scores(
  data, eig, K_grid = K_grid,
  inverse_strategy = "ginv",
  ginv_tol = 1e-9
)

manual_ginv <- t(sweep(eig$funcs$Theta, 2, eig$values$lam_theta, "*")) %*%
  ginv_spd(K_grid, tol = 1e-3, ridge = 1e-8) %*%
  c(1, 1)

stopifnot(
  isTRUE(all.equal(score_default@data, score_auto@data)),
  isTRUE(all.equal(score_auto@data, score_ginv@data)),
  isTRUE(all.equal(as.numeric(score_ginv@data), as.numeric(manual_ginv))),
  score_ginv@data[1, 1, 1, 2] == 0,
  score_ridge@data[1, 1, 1, 2] > 0,
  score_less_truncated@data[1, 1, 1, 2] > 0
)

# Positive noise keeps the existing ridge BLUP under the automatic strategy.
score_auto_positive <- get_matfpca_scores(
  data, eig, sig2 = 0.1, K_grid = K_grid,
  inverse_strategy = "auto"
)
score_ridge_positive <- get_matfpca_scores(
  data, eig, sig2 = 0.1, K_grid = K_grid,
  inverse_strategy = "ridge"
)
stopifnot(isTRUE(all.equal(
  score_auto_positive@data,
  score_ridge_positive@data
)))

stopifnot(
  inherits(
    try(get_matfpca_scores(
      data, eig, K_grid = K_grid,
      inverse_strategy = "unknown"
    ), silent = TRUE),
    "try-error"
  ),
  inherits(
    try(get_matfpca_scores(
      data, eig, K_grid = K_grid,
      ginv_tol = 1
    ), silent = TRUE),
    "try-error"
  )
)

cat("Score inverse-strategy checks passed.\n")
