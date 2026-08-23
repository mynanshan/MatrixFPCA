#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]))
} else {
  normalizePath(file.path("test", "test-adni-latent-utils.R"))
}
project_root <- normalizePath(file.path(dirname(script_path), ".."))
setwd(project_root)

suppressPackageStartupMessages({
  library(data.table)
  library(rTensor)
})
source(file.path("R", "matfpca.R"))
source(file.path("application", "adni-latent", "adni_latent_utils.R"))

# The R reader must preserve NumPy C-order rows without reticulate/Python.
write_test_npy <- function(path, matrix_value) {
  base_header <- sprintf(
    "{'descr': '<f4', 'fortran_order': False, 'shape': (%d, %d), }",
    nrow(matrix_value), ncol(matrix_value)
  )
  padding <- (64L - ((10L + nchar(base_header) + 1L) %% 64L)) %% 64L
  header <- paste0(base_header, strrep(" ", padding), "\n")
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(as.raw(c(0x93, utf8ToInt("NUMPY"))), connection)
  writeBin(as.raw(c(1L, 0L)), connection)
  writeBin(nchar(header), connection, size = 2L, endian = "little")
  writeBin(charToRaw(header), connection)
  writeBin(as.numeric(t(matrix_value)), connection, size = 4L, endian = "little")
}

expected <- matrix(as.numeric(seq_len(12L)), nrow = 3L, byrow = TRUE)
npy_path <- tempfile(fileext = ".npy")
on.exit(unlink(npy_path), add = TRUE)
write_test_npy(npy_path, expected)
stopifnot(identical(read_npy_matrix(npy_path), expected))

# Feature PCA gives every subject equal aggregate weight even when visit
# counts differ. Subject A has two visits and subject B has one.
projection_paths <- vapply(seq_len(3L), function(index) {
  path <- tempfile(fileext = ".npy")
  write_test_npy(path, list(
    matrix(c(0, 0), nrow = 1L),
    matrix(c(2, 0), nrow = 1L),
    matrix(c(0, 4), nrow = 1L)
  )[[index]])
  path
}, character(1L))
on.exit(unlink(projection_paths), add = TRUE)
projection_metadata <- data.table(
  subject_id = c("A", "A", "B"),
  latent_path = projection_paths
)
projection <- fit_feature_projection(projection_metadata, pve = 1, rank_cap = 2L)
stopifnot(
  isTRUE(projection$subject_weighted),
  all.equal(projection$center, c(0.5, 2), tolerance = 1e-12),
  projection$training_token_rows == 3L,
  projection$effective_weighted_token_rows == 2,
  projection$retained_pve >= projection$requested_pve - 1e-12,
  positive_mode_rank(c(3, 1, 1e-15)) == 2L
)

# MONAI flattening has width varying fastest, then height, then depth.
cube <- token_vector_to_cube(seq_len(216L))
stopifnot(
  identical(dim(cube), c(6L, 6L, 6L)),
  cube[1L, 1L, 1L] == 1L,
  cube[1L, 1L, 2L] == 2L,
  cube[1L, 2L, 1L] == 7L,
  cube[2L, 1L, 1L] == 37L
)

# Temporal eigenfunctions must be orthonormal under irregular-grid trapezoid
# quadrature, not under a uniform-grid scaling.
grid <- c(0, 0.25, 1)
K <- exp(-outer(grid, grid, function(left, right) (left - right)^2))
covariance <- create_matrix_covfunc(diag(c(2, 1)), diag(c(3, 2, 1)), K, grid)
eig <- eigen_matfcov_quadrature(covariance)
weights <- trapezoid_weights(grid)
inner_products <- crossprod(eig$funcs$Theta, weights * eig$funcs$Theta)
stopifnot(
  max(abs(inner_products - diag(nrow(inner_products)))) < 1e-10,
  abs(sum(weights) - diff(range(grid))) < 1e-12,
  all(diff(eig$values$lam_theta) <= 1e-12)
)

cat("ADNI latent utility regression tests passed.\n")
