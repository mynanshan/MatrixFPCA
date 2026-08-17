# fdarep::SparseProductFPCA is a two-stage sparse product FPCA. It treats the
# vectorized image as a one-dimensional spatial function; this is deliberately
# recorded because the artificial pixel ordering is not invariant.

fit_sparse_product_fpca <- function(data, ncomp = NULL, spatial_npc = 5L,
                                    temporal_npc = 3L) {
  .require_competitor_package("fdarep", "sparse_product_fpca")
  longitudinal <- .as_longitudinal_arrays(data)
  Ly <- longitudinal$Ly
  p <- dim(Ly[[1L]])[2L]
  q <- dim(Ly[[1L]])[3L]
  npix <- p * q
  spatial_npc <- min(npix, max(1L, as.integer(spatial_npc)))
  temporal_npc <- max(1L, as.integer(temporal_npc))
  spatial_grid <- seq(0, 1, length.out = npix)
  image_list <- lapply(Ly, function(Y) {
    # columns are visits, rows use R's column-major image vectorization.
    t(matrix(Y, nrow = nrow(Y), ncol = npix))
  })

  fit <- fdarep::SparseProductFPCA(
    sSup = spatial_grid,
    Lt = longitudinal$Lt,
    Ly = image_list,
    fpca.op1 = list(
      maxK = spatial_npc,
      FVEthreshold = 0.999,
      usergrid = TRUE,
      methodSelectK = "FVE"
    ),
    fpca.op2 = list(
      maxK = temporal_npc,
      FVEthreshold = 0.999,
      usergrid = FALSE,
      methodSelectK = "FVE"
    ),
    pc.j = spatial_npc,
    pc.k = temporal_npc
  )
  spatial_phi <- as.matrix(fit$psi)
  temporal_phi <- as.matrix(fit$phi)
  temporal_phi <- .interpolate_columns(fit$year.grid, data$t_grid, temporal_phi)
  if (nrow(spatial_phi) != npix) {
    spatial_phi <- .interpolate_columns(fit$s.grid, spatial_grid, spatial_phi)
  }

  components <- list()
  eigenvalues <- numeric()
  index <- 0L
  for (j in seq_len(min(ncol(spatial_phi), length(fit$VarOrdered)))) {
    lambda_j <- as.numeric(fit$VarOrdered[[j]])
    nk <- min(ncol(temporal_phi), length(lambda_j), temporal_npc)
    for (k in seq_len(nk)) {
      if (!is.finite(lambda_j[k]) || lambda_j[k] <= 0) next
      index <- index + 1L
      components[[index]] <- array(
        temporal_phi[, k] %o% spatial_phi[, j],
        c(length(data$t_grid), p, q)
      )
      eigenvalues[index] <- lambda_j[k]
    }
  }
  if (!length(components)) stop("SparseProductFPCA returned no positive component.")
  psi <- array(unlist(components, use.names = FALSE), c(
    length(data$t_grid), p, q, length(components)
  ))
  answer <- list(
    psi = psi,
    nu = eigenvalues,
    sig2 = NA_real_,
    tgrid = data$t_grid,
    diagnostics = list(
      package_version = as.character(utils::packageVersion("fdarep")),
      spatial_npc = spatial_npc,
      temporal_npc = temporal_npc,
      spatial_order = "column-major vectorization of the image",
      raw_fit = fit
    )
  )
  .truncate_eigenfit(answer, ncomp)
}
