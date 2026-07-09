source("matfpca.R")
source(file.path("experiment", "gendata-simple.R"))

set.seed(20260625)

generated <- gen_data_simple(
  n = 20,
  p = 6,
  q = 5,
  npc_row = 6,
  npc_col = 5,
  npc_time = 3,
  nt = 21,
  sampling_ratio = c(0.20, 0.35),
  sig = 0.5,
  eigvecSet = 2,
  scores_distr = "normal",
  score_dependence = "common_amplitude",
  dependence_param = 1
)

stopifnot(
  identical(dim(generated$scores), c(20L, 6L, 5L, 3L)),
  abs(mean(generated$subject_amplitude^2) - 1) < 0.75,
  !"maxIter" %in% names(formals(fit_matfcov)),
  "temporal_method" %in% names(formals(fit_matfcov))
)

fit_identity <- fit_matfcov(
  generated$data,
  wScheme = "identity",
  noise = TRUE,
  nbasis = "theory"
)
fit_refined <- fit_matfcov(
  generated$data,
  wScheme = "gaussian_refined",
  noise = TRUE,
  nbasis = "theory",
  opwIter = 2
)
fit_joint_identity <- fit_matfcov(
  generated$data,
  wScheme = "identity",
  noise = TRUE,
  nbasis = "theory",
  temporal_method = "joint"
)
fit_joint_refined <- fit_matfcov(
  generated$data,
  wScheme = "gaussian_refined",
  noise = TRUE,
  nbasis = "theory",
  opwIter = 2,
  temporal_method = "joint"
)

for (fit in list(fit_identity, fit_refined, fit_joint_identity, fit_joint_refined)) {
  covariance <- fit$C
  stopifnot(
    abs(tr(covariance$U) - 6) < 1e-8,
    abs(tr(covariance$V) - 5) < 1e-8,
    min(eigen(covariance$U, symmetric = TRUE, only.values = TRUE)$values) > -1e-8,
    min(eigen(covariance$V, symmetric = TRUE, only.values = TRUE)$values) > -1e-8,
    min(eigen(covariance$K, symmetric = TRUE, only.values = TRUE)$values) > -1e-8,
    covariance$sig2 >= 0,
    identical(names(fit$temporal_fit$fit), c("K", "G")),
    identical(names(fit$temporal_fit$nbasis), c("K", "G")),
    fit$temporal_fit$temporal_method %in% c("two_step", "joint")
  )
}

stopifnot(
  is.null(fit_identity$W),
  inherits(fit_refined$W, "mode_weights"),
  inherits(fit_refined$W, "gaussian_refined"),
  is(fit_refined$pilot, "matfcov")
)

# At exactly zero noise, the Gaussian-refined formula must reduce to the
# plug-in inverse-covariance scheme.
U_test <- diag(c(2, 1))
V_test <- diag(c(3, 1.5, 0.75))
LK_test <- list(diag(c(2, 1)), diag(c(1.5, 0.5, 0.25)))
weights_refined_zero <- compute_W(
  U_test, V_test, LK_test,
  sig2 = 0,
  wScheme = "gaussian_refined"
)
weights_inverse <- compute_W(
  U_test, V_test, LK_test,
  sig2 = 0,
  wScheme = "invcov"
)
stopifnot(
  inherits(weights_refined_zero, "shared"),
  isTRUE(all.equal(weights_refined_zero, weights_inverse))
)

cat("Temporal estimator checks passed.\n")
