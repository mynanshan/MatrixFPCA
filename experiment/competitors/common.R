# Shared adapters for methods that expect vector-valued longitudinal data.

.require_competitor_package <- function(package, method) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      method, " requires package '", package, "'. Run ",
      "Rscript experiment/competitors/install-dependencies.R --r-only.",
      call. = FALSE
    )
  }
}

.pixel_long_data <- function(data) {
  longitudinal <- .as_longitudinal_arrays(data)
  Ly <- longitudinal$Ly
  Lt <- longitudinal$Lt
  n <- length(Ly)
  p <- dim(Ly[[1L]])[2L]
  q <- dim(Ly[[1L]])[3L]

  outcomes <- vector("list", p * q)
  names(outcomes) <- sprintf("pixel_%03d", seq_along(outcomes))
  outcome <- 0L
  # This order matches matrix(Y, nrow = n_time): row index varies first.
  for (k in seq_len(q)) {
    for (j in seq_len(p)) {
      outcome <- outcome + 1L
      outcomes[[outcome]] <- do.call(rbind, lapply(seq_len(n), function(i) {
        data.frame(
          subj = i,
          argvals = Lt[[i]],
          y = as.numeric(Ly[[i]][, j, k])
        )
      }))
    }
  }
  outcomes
}

.reshape_outcome_eigenfunctions <- function(eigenfunctions, nt, p, q) {
  eigenfunctions <- as.matrix(eigenfunctions)
  if (nrow(eigenfunctions) != nt * p * q) {
    stop(
      "Expected ", nt * p * q, " eigenfunction rows but received ",
      nrow(eigenfunctions), "."
    )
  }
  array(eigenfunctions, c(nt, p, q, ncol(eigenfunctions)))
}

.truncate_eigenfit <- function(fit, ncomp = NULL) {
  fit <- canonicalize_eigenfit(fit)
  if (is.null(ncomp) || length(fit$nu) <= ncomp) return(fit)
  keep <- seq_len(max(1L, as.integer(ncomp)))
  fit$psi <- fit$psi[, , , keep, drop = FALSE]
  fit$nu <- fit$nu[keep]
  fit
}

.interpolate_columns <- function(from, to, values) {
  values <- as.matrix(values)
  if (length(from) == length(to) && isTRUE(all.equal(as.numeric(from), as.numeric(to)))) {
    return(values)
  }
  vapply(seq_len(ncol(values)), function(component) {
    stats::approx(from, values[, component], xout = to, rule = 2L)$y
  }, numeric(length(to)))
}

.integer_candidates <- function(x, name) {
  x <- unique(as.integer(x))
  if (!length(x) || any(!is.finite(x)) || any(x < 1L)) {
    stop(name, " must contain positive integers.")
  }
  x
}
