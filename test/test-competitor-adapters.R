library(rTensor)

source(file.path("R", "matfdata.R"))
source(file.path("R", "matfpca.R"))
source(file.path("simulation", "gendata-simple.R"))
source(file.path("simulation", "competitors", "core.R"))

set.seed(19)
generated <- gen_data_simple(
  n = 12L, p = 3L, q = 3L, npc_row = 3L, npc_col = 3L, npc_time = 5L,
  nt = 31L, t_range = c(0, 1), sampling_ratio = c(0.15, 0.30), sig = 0.1,
  eigvecSet = 2L, eigvalSet = 1L, spatial_model = "ar1",
  scores_distr = "normal", score_dependence = "independent",
  dependence_param = 1
)

fits <- list(
  shi_mdfpca = fit_shi_mdfpca(
    generated$data, ncomp = 2L, nbasis_time = 5L,
    nbasis_row = 3L, nbasis_col = 3L, maxit = 20L,
    score_ridge = "auto"
  ),
  spaco_minus = fit_spaco_minus(
    generated$data, ncomp = 2L, nbasis = 5L, maxit = 10L
  )
)
if (requireNamespace("face", quietly = TRUE)) {
  fits$happ_face <- fit_happ_face(
    generated$data, npc = 2L, ncomp = 2L, knots = 3L,
    search_length = 3L, two_step = FALSE
  )
}
if (requireNamespace("mfaces", quietly = TRUE)) {
  fits$mfaces <- fit_mfaces(
    generated$data, ncomp = 2L, knots = 3L, search_length = 3L
  )
}
if (requireNamespace("mGSFPCA", quietly = TRUE)) {
  fits$mgsfpca <- fit_mgsfpca(
    generated$data, ncomp = 2L, rank_candidates = 2L,
    basis_candidates = 5L, mu_nbasis = 5L, maxit = 30L
  )
}
if (requireNamespace("fdarep", quietly = TRUE)) {
  fits$sparse_product_fpca <- fit_sparse_product_fpca(
    generated$data, ncomp = 2L, spatial_npc = 2L, temporal_npc = 2L
  )
}

for (method in names(fits)) {
  fit <- fits[[method]]
  stopifnot(
    identical(dim(fit$psi)[1L:3L], c(31L, 3L, 3L)),
    length(fit$nu) == dim(fit$psi)[4L],
    all(is.finite(fit$psi)), all(is.finite(fit$nu)), all(fit$nu > 0)
  )
  weights <- rep(.trapz_weights(fit$tgrid), 9L)
  functions <- matrix(fit$psi, nrow = length(weights))
  gram <- crossprod(sqrt(weights) * functions)
  stopifnot(max(abs(gram - diag(ncol(gram)))) < 1e-7)
  metrics <- evaluate_full_eigenfit(fit, generated$cov, eval_ncomp = 2L)
  stopifnot(all(is.finite(metrics[c("rel_hs_cov", "mean_sin2")])))
}
stopifnot(
  is.finite(fits$shi_mdfpca$diagnostics$selected_score_ridge),
  fits$shi_mdfpca$diagnostics$selected_score_ridge >= 0
)
auto_spaco <- fit_spaco_minus(
  generated$data, ncomp = 2L, nbasis = 5L, maxit = 10L,
  score_ridge = "auto"
)
stopifnot(
  is.finite(auto_spaco$diagnostics$selected_score_ridge),
  auto_spaco$diagnostics$selected_score_ridge >= 0,
  length(auto_spaco$diagnostics$score_ridge_cv_error) == 18L,
  length(auto_spaco$diagnostics$score_ridge_cv_se) == 18L
)

cat("competitor adapter smoke tests passed\n")
