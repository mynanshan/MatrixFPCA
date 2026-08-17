# Zhang et al.'s multivariate generalized sparse FPCA (mGSFPCA), used here
# with Gaussian outcomes. Each image pixel is one functional coordinate.

fit_mgsfpca <- function(data, ncomp = NULL, rank_candidates = 3L,
                        basis_candidates = 7L, mu_nbasis = 7L,
                        maxit = 300L, optim_tol = 1e-5,
                        bin_size = NULL) {
  .require_competitor_package("mGSFPCA", "mgsfpca")
  longitudinal <- .as_longitudinal_arrays(data)
  p <- dim(longitudinal$Ly[[1L]])[2L]
  q <- dim(longitudinal$Ly[[1L]])[3L]
  noutcome <- p * q
  ranks <- .integer_candidates(rank_candidates, "rank_candidates")
  bases <- .integer_candidates(basis_candidates, "basis_candidates")
  outcomes <- lapply(.pixel_long_data(data), function(outcome) {
    data.frame(ID = outcome$subj, time = outcome$argvals, value = outcome$y)
  })
  if (is.null(bin_size)) bin_size <- length(data$t_grid)

  fit <- mGSFPCA::spMultFPCA(
    dataCell = outcomes,
    r = rep(list(ranks), noutcome),
    G = rep(list(bases), noutcome),
    basis_type = "bspline",
    nRegGrid = length(data$t_grid),
    grid_range = rep(list(range(data$t_grid)), noutcome),
    mu_nbasis = max(4L, as.integer(mu_nbasis)),
    M_npc = NULL,
    maxit = max(1L, as.integer(maxit)),
    optim_tol = optim_tol,
    optim_trace = 0L,
    bin_size = max(5L, as.integer(bin_size))
  )
  if (is.null(fit$fullEigFun) || is.null(fit$values)) {
    stop("mGSFPCA returned no multivariate eigendecomposition.")
  }
  noise <- vapply(fit$Ckk, function(univariate_fit) {
    candidates <- c(univariate_fit$sigma2, univariate_fit$sig2,
                    univariate_fit$var_error, univariate_fit$error_var)
    candidates <- as.numeric(candidates)
    candidates <- candidates[is.finite(candidates) & candidates >= 0]
    if (length(candidates)) candidates[1L] else NA_real_
  }, numeric(1L))
  selected_rank <- vapply(fit$Ckk, function(x) {
    if (!is.null(x$pars$p)) as.integer(x$pars$p) else NA_integer_
  }, integer(1L))
  selected_basis <- vapply(fit$Ckk, function(x) {
    if (!is.null(x$pars$k)) as.integer(x$pars$k) else NA_integer_
  }, integer(1L))
  answer <- list(
    psi = .reshape_outcome_eigenfunctions(
      fit$fullEigFun, length(data$t_grid), p, q
    ),
    nu = as.numeric(fit$values),
    sig2 = if (any(is.finite(noise))) mean(noise, na.rm = TRUE) else NA_real_,
    tgrid = data$t_grid,
    diagnostics = list(
      package_version = as.character(utils::packageVersion("mGSFPCA")),
      rank_candidates = ranks,
      basis_candidates = bases,
      selected_rank = selected_rank,
      selected_basis = selected_basis,
      raw_fit = fit
    )
  )
  .truncate_eigenfit(answer, ncomp)
}
