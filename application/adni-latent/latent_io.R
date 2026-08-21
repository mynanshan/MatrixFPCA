read_npy_matrix <- function(path, require_finite = TRUE) {
  if (!file.exists(path)) stop("Latent matrix does not exist: ", path)
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)

  magic <- readBin(connection, what = "raw", n = 6L)
  expected_magic <- as.raw(c(0x93, utf8ToInt("NUMPY")))
  if (!identical(magic, expected_magic)) stop("Not a NumPy .npy file: ", path)

  version <- readBin(connection, what = "integer", n = 2L, size = 1L,
                     signed = FALSE, endian = "little")
  if (version[[1L]] == 1L) {
    header_length <- readBin(
      connection, what = "integer", n = 1L, size = 2L,
      signed = FALSE, endian = "little"
    )
  } else if (version[[1L]] %in% c(2L, 3L)) {
    header_length <- readBin(
      connection, what = "integer", n = 1L, size = 4L,
      signed = TRUE, endian = "little"
    )
  } else {
    stop("Unsupported .npy version ", paste(version, collapse = "."), ": ", path)
  }
  header <- rawToChar(readBin(connection, what = "raw", n = header_length))

  capture <- function(pattern) {
    match <- regexec(pattern, header, perl = TRUE)
    value <- regmatches(header, match)[[1L]]
    if (length(value) < 2L) stop("Malformed .npy header in ", path)
    value[[2L]]
  }
  descriptor <- capture("['\"]descr['\"]\\s*:\\s*['\"]([^'\"]+)['\"]")
  fortran_order <- capture("['\"]fortran_order['\"]\\s*:\\s*(True|False)")
  shape_text <- capture("['\"]shape['\"]\\s*:\\s*\\(([^)]*)\\)")
  if (fortran_order != "False") stop("Fortran-order .npy arrays are unsupported: ", path)

  shape <- as.integer(trimws(strsplit(shape_text, ",", fixed = TRUE)[[1L]]))
  shape <- shape[!is.na(shape)]
  if (length(shape) != 2L || any(shape < 1L)) {
    stop("Expected a two-dimensional nonempty .npy array: ", path)
  }
  if (!descriptor %in% c("<f4", "|f4", "=f4", "<f8", "|f8", "=f8", ">f4", ">f8")) {
    stop("Expected float32 or float64 .npy data, found ", descriptor, ": ", path)
  }
  size <- if (endsWith(descriptor, "f4")) 4L else 8L
  endian <- if (startsWith(descriptor, ">")) "big" else "little"
  values <- readBin(
    connection, what = "numeric", n = prod(shape), size = size,
    signed = TRUE, endian = endian
  )
  if (length(values) != prod(shape)) stop("Truncated .npy payload: ", path)
  result <- matrix(values, nrow = shape[[1L]], ncol = shape[[2L]], byrow = TRUE)
  if (require_finite && any(!is.finite(result))) {
    stop("Non-finite values in latent matrix: ", path)
  }
  result
}

resolve_latent_path <- function(path, metadata_path) {
  if (grepl("^/", path)) return(path)
  normalizePath(file.path(dirname(metadata_path), path), mustWork = FALSE)
}

token_vector_to_cube <- function(values, token_grid = c(6L, 6L, 6L)) {
  token_grid <- as.integer(token_grid)
  if (length(values) != prod(token_grid)) {
    stop("Token vector length does not match the requested three-dimensional grid.")
  }
  # MONAI flattens a D x H x W convolutional grid with W varying fastest.
  aperm(array(values, dim = rev(token_grid)), c(3L, 2L, 1L))
}
