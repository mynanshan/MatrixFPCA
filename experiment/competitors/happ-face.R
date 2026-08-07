# Happ--Greven MFPCA with FACE rather than PACE for each pixel's univariate
# expansion. This isolates whether HAPP's error is driven by local PACE fits.

fit_happ_face <- function(data, npc = 3L, ncomp = NULL, knots = 5L,
                          pve = 0.999, search_length = 10L,
                          two_step = TRUE, ncores = 1L,
                          backend = c("multisession", "multicore")) {
  .require_competitor_package("face", "happ_face")
  backend <- match.arg(backend)
  longitudinal <- .as_longitudinal_arrays(data)
  Ly <- longitudinal$Ly
  n <- length(Ly)
  p <- dim(Ly[[1L]])[2L]
  q <- dim(Ly[[1L]])[3L]
  nt <- length(data$t_grid)
  npc <- max(1L, as.integer(npc))
  weights <- .trapz_weights(data$t_grid)
  outcomes <- .pixel_long_data(data)

  pixel_task <- function(outcome_index) {
    k <- (outcome_index - 1L) %/% p + 1L
    j <- (outcome_index - 1L) %% p + 1L
    zero_fit <- function(status) {
      .zero_pixel_fit(j, k, nt, n, npc, data$t_grid, weights, status)
    }
    outcome <- outcomes[[outcome_index]]
    if (max(abs(outcome$y - mean(outcome$y))) <=
        1e-12 + 1e-8 * max(1, abs(mean(outcome$y)))) {
      return(zero_fit("degenerate"))
    }
    tryCatch({
      local <- suppressWarnings(face::face.sparse(
        data = outcome,
        argvals.new = data$t_grid,
        knots = max(2L, as.integer(knots)),
        search.length = max(2L, as.integer(search_length)),
        calculate.scores = TRUE,
        pve = pve,
        two_step = two_step
      ))
      phi_raw <- as.matrix(local$eigenfunctions)
      lambda_raw <- pmax(as.numeric(local$eigenvalues), 0)
      scores_raw <- as.matrix(local$rand_eff$scores)
      available <- min(ncol(phi_raw), length(lambda_raw), ncol(scores_raw), npc)
      if (available < 1L || !any(lambda_raw[seq_len(available)] > 0)) {
        return(zero_fit("zero_eigen"))
      }
      phi_raw <- phi_raw[, seq_len(available), drop = FALSE]
      qr_fit <- qr(sqrt(weights) * phi_raw)
      phi <- qr.Q(qr_fit) / sqrt(weights)
      rotation <- qr.R(qr_fit)
      scores <- scores_raw[, seq_len(available), drop = FALSE] %*%
        t(rotation[seq_len(available), seq_len(available), drop = FALSE])
      if (available < npc) {
        phi <- cbind(phi, .weighted_completion(
          phi, data$t_grid, weights, npc - available
        ))
        scores <- cbind(scores, matrix(0, n, npc - available))
      }
      sig2 <- mean(as.numeric(local$var.error.new), na.rm = TRUE)
      if (!is.finite(sig2)) sig2 <- NA_real_
      list(
        j = j, k = k, phi = phi[, seq_len(npc), drop = FALSE],
        lambda = c(lambda_raw[seq_len(available)], rep(0, npc - available)),
        sig2 = sig2, scores = scores[, seq_len(npc), drop = FALSE],
        status = "ok"
      )
    }, error = function(error) {
      answer <- zero_fit("face_error")
      answer$error <- conditionMessage(error)
      answer
    })
  }

  indices <- seq_along(outcomes)
  pixel_fits <- if (ncores <= 1L) {
    lapply(indices, pixel_task)
  } else {
    .require_competitor_package("future", "parallel happ_face")
    .require_competitor_package("future.apply", "parallel happ_face")
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    strategy <- if (backend == "multicore") future::multicore else future::multisession
    future::plan(strategy, workers = as.integer(ncores))
    future.apply::future_lapply(
      indices, pixel_task, future.seed = TRUE, future.packages = "face"
    )
  }

  local_phi <- array(0, c(nt, p, q, npc))
  score_matrix <- matrix(0, n, p * q * npc)
  for (local in pixel_fits) {
    local_phi[, local$j, local$k, ] <- local$phi
    block <- ((local$j - 1L) * q + local$k - 1L) * npc + seq_len(npc)
    score_matrix[, block] <- local$scores
  }
  score_matrix <- scale(score_matrix, center = TRUE, scale = FALSE)
  score_svd <- svd(score_matrix / sqrt(n - 1L), nu = 0L)
  tolerance_svd <- max(dim(score_matrix)) * .Machine$double.eps * max(score_svd$d, 0)
  rank_effective <- sum(score_svd$d > tolerance_svd)
  if (rank_effective < 1L) stop("The HAPP-FACE score covariance has zero rank.")
  if (is.null(ncomp)) ncomp <- rank_effective
  ncomp <- min(max(1L, as.integer(ncomp)), rank_effective)
  coefficients <- score_svd$v[, seq_len(ncomp), drop = FALSE]
  eigenfunctions <- array(0, c(nt, p, q, ncomp))
  for (j in seq_len(p)) for (k in seq_len(q)) {
    block <- ((j - 1L) * q + k - 1L) * npc + seq_len(npc)
    eigenfunctions[, j, k, ] <- local_phi[, j, k, ] %*%
      coefficients[block, , drop = FALSE]
  }
  status <- vapply(pixel_fits, `[[`, character(1L), "status")
  noise <- vapply(pixel_fits, `[[`, numeric(1L), "sig2")
  .truncate_eigenfit(list(
    psi = eigenfunctions,
    nu = score_svd$d[seq_len(ncomp)]^2,
    sig2 = if (any(is.finite(noise))) mean(noise, na.rm = TRUE) else NA_real_,
    tgrid = data$t_grid,
    diagnostics = list(
      package_version = as.character(utils::packageVersion("face")),
      pixel_status = table(status), knots = knots, pve = pve,
      two_step = two_step
    )
  ), ncomp)
}
