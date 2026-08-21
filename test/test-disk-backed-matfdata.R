#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rTensor))
source(file.path("R", "matfdata.R"))

temporary_directory <- tempfile("matrixfpca-disk-backed-")
dir.create(temporary_directory)
on.exit(unlink(temporary_directory, recursive = TRUE), add = TRUE)

paths <- file.path(temporary_directory, paste0("subject-", 1:2, ".rds"))
saveRDS(as.tensor(array(seq_len(24), dim = c(3, 4, 2))), paths[[1L]])
saveRDS(as.tensor(array(seq_len(36), dim = c(3, 4, 3))), paths[[2L]])

data <- create_matrix_funcdata_irreg(
  X_list = paths,
  t_grid = c(0, 1, 2),
  t_idx_list = list(c(1L, 3L), c(1L, 2L, 3L))
)

stopifnot(
  identical(data$vdim, c(3L, 4L)),
  identical(data$in_mem, FALSE),
  identical(dim(get_subj_tnsr(data, 1L)), c(3L, 4L, 2L)),
  identical(dim(get_subj_tnsr(data, 2L)), c(3L, 4L, 3L))
)

message("Disk-backed matfdata regression test passed.")
