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
  !"maxIter" %in% names(formals(fit_matfcov))
)

fit_identity <- fit_matfcov(
  generated$data,
  wScheme = "identity",
  noise = TRUE,
  nbasis = "theory"
)
fit_optimal <- fit_matfcov(
  generated$data,
  wScheme = "optimal",
  noise = TRUE,
  nbasis = "theory",
  opwIter = 2
)

for (fit in list(fit_identity, fit_optimal)) {
  covariance <- fit$C
  stopifnot(
    abs(tr(covariance$U) - 6) < 1e-8,
    abs(tr(covariance$V) - 5) < 1e-8,
    min(eigen(covariance$U, symmetric = TRUE, only.values = TRUE)$values) > -1e-8,
    min(eigen(covariance$V, symmetric = TRUE, only.values = TRUE)$values) > -1e-8,
    min(eigen(covariance$K, symmetric = TRUE, only.values = TRUE)$values) > -1e-8,
    covariance$sig2 >= 0,
    identical(names(fit$temporal_fit$fit), c("K", "G")),
    identical(names(fit$temporal_fit$nbasis), c("K", "G"))
  )
}

stopifnot(
  is.null(fit_identity$W),
  inherits(fit_optimal$W, "mode_weights"),
  is(fit_optimal$pilot, "matfcov")
)

cat("Two-step estimator checks passed.\n")
