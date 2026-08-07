# Li, Xiao and Luo (2020), fast covariance estimation for multivariate sparse
# functional data. Every image pixel is treated as one functional outcome.

fit_mfaces <- function(data, ncomp = NULL, knots = 5L, pve = 0.999,
                       search_length = 10L, lambda = NULL,
                       lambda_mean = NULL, lambda_cross = NULL,
                       max_outcomes = 25L, allow_large = FALSE) {
  .require_competitor_package("mfaces", "mfaces")
  longitudinal <- .as_longitudinal_arrays(data)
  p <- dim(longitudinal$Ly[[1L]])[2L]
  q <- dim(longitudinal$Ly[[1L]])[3L]
  noutcome <- p * q
  max_outcomes <- max(1L, as.integer(max_outcomes))
  if (!allow_large && noutcome > max_outcomes) {
    stop(
      "mfaces was not run for ", noutcome, " pixel outcomes. It smooths ",
      noutcome * (noutcome + 1L) / 2L, " covariance surfaces; the current ",
      "safety limit is ", max_outcomes, " outcomes. Reduce p*q or pass ",
      "--allow_slow_competitors 1.", call. = FALSE
    )
  }
  if (!is.finite(pve) || pve <= 0 || pve > 1) stop("pve must be in (0, 1].")

  outcomes <- .pixel_long_data(data)
  fit <- mfaces::mface.sparse(
    data = outcomes,
    argvals.new = data$t_grid,
    knots = max(2L, as.integer(knots)),
    lambda = lambda,
    lambda_mean = lambda_mean,
    lambda_bps = lambda_cross,
    search.length = max(2L, as.integer(search_length)),
    calculate.scores = FALSE,
    pve = pve
  )
  if (is.null(fit$eigenfunctions) || is.null(fit$eigenvalues)) {
    stop("mfaces returned no multivariate eigendecomposition.")
  }
  sig2 <- mean(as.numeric(fit$var.error.new), na.rm = TRUE)
  if (!is.finite(sig2)) sig2 <- NA_real_
  answer <- list(
    psi = .reshape_outcome_eigenfunctions(
      fit$eigenfunctions, length(data$t_grid), p, q
    ),
    nu = as.numeric(fit$eigenvalues),
    sig2 = sig2,
    tgrid = data$t_grid,
    diagnostics = list(
      package_version = as.character(utils::packageVersion("mfaces")),
      outcomes = noutcome,
      covariance_surfaces = noutcome * (noutcome + 1L) / 2L,
      knots = knots,
      search_length = search_length,
      raw_fit = fit
    )
  )
  .truncate_eigenfit(answer, ncomp)
}
