library(rTensor)

source("matfdata.R")
source("matfpca.R")
source(file.path("experiment", "competitors", "core.R"))

tgrid <- seq(0, 1, length.out = 21)
theta <- sqrt(2) * sin(pi * tgrid)
K <- tcrossprod(theta)
U <- diag(c(3, 1))
V <- diag(c(3, 1))
truth <- create_matrix_covfunc(U, V, K, tgrid, sig2 = 0)

# A complete product eigendecomposition must reconstruct the covariance under
# the time-quadrature Hilbert metric.
truth_fit <- product_covariance_eigenfit(truth, ncomp = 4L)
truth_metrics <- evaluate_full_eigenfit(truth_fit, truth, eval_ncomp = 2L)
stopifnot(
  truth_metrics[["rel_hs_cov"]] < 1e-10,
  truth_metrics[["mean_sin2"]] < 1e-10,
  truth_metrics[["rel_eigenvalue_sse"]] < 1e-10
)

# A fixed-rank comparison must truncate a product-covariance estimator too;
# otherwise MatFPCA would be scored at full rank while competitors are scored
# at the rank supplied through --cov_ncomp.
rank_two_product_metrics <- evaluate_product_fit(
  list(C = truth), truth, eval_ncomp = 2L, covariance_ncomp = 2L
)

stopifnot(
  rank_two_product_metrics[["estimated_rank"]] == 2L,
  rank_two_product_metrics[["oracle_rank_rel_hs"]] > 0,
  abs(rank_two_product_metrics[["rel_hs_cov"]] -
      rank_two_product_metrics[["oracle_rank_rel_hs"]]) < 1e-10
)

# The requested rank two cuts the repeated eigenvalue 3 * 1 = 1 * 3, so the
# identifiable evaluation target must include the complete three-dimensional
# leading eigenspace.
stopifnot(
  .identifiable_evaluation_rank(truth, requested = 2L) == 3L,
  truth_metrics[["evaluated_components"]] == 3L
)

# Signs and labels of individual FPCs cannot change either covariance or
# eigenspace accuracy.
permutation <- c(3L, 1L, 2L, 4L)
equivalent_fit <- truth_fit
equivalent_fit$psi <- equivalent_fit$psi[, , , permutation, drop = FALSE]
equivalent_fit$psi[, , , 1L] <- -equivalent_fit$psi[, , , 1L]
equivalent_fit$nu <- equivalent_fit$nu[permutation]
equivalent_metrics <- evaluate_full_eigenfit(equivalent_fit, truth, eval_ncomp = 2L)
stopifnot(
  equivalent_metrics[["rel_hs_cov"]] < 1e-10,
  equivalent_metrics[["mean_sin2"]] < 1e-10
)

# Rescaling the three product factors without changing their Kronecker product
# must have exactly zero covariance error.
rescaled <- create_matrix_covfunc(2 * U, V / 3, 1.5 * K, tgrid, sig2 = 0)
stopifnot(matfcovDiffNorm2(rescaled, truth) < 1e-10)

cat("competitor method metric tests passed\n")
