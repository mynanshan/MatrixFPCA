#!/usr/bin/env Rscript
# EchoNet RDS conversion (step 2): one rTensor (p x q x n_i) per subject.
#
#   Rscript application/echonet/echonet_to_rds.R <export_dir>
#
# Reads frames.f32 / frames.csv / subjects.csv / mask.csv written by
# export_echonet.py and writes <export_dir>/rds/<id>.rds. Each tensor holds
# all retained systolic frames of one video, ordered by phase, already
# gain-normalised and masked (not centred). The frame table gains an
# `rds_file` column and a 1-based `t_idx` into the 51-point phase grid, which
# is what create_matrix_funcdata_irreg() expects.

suppressPackageStartupMessages({
  library(data.table)
  library(rTensor)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: echonet_to_rds.R <export_dir>")
export_dir <- normalizePath(args[[1L]])
read_export_config <- function(path) {
  # tiny reader for the flat export.json (avoids a jsonlite dependency)
  txt <- paste(readLines(path), collapse = "")
  get <- function(key) as.numeric(sub(
    sprintf('.*"%s": *([0-9.]+).*', key), "\\1", txt))
  list(p = get("p"), q = get("q"), n_frames = get("n_frames"))
}
cfg <- read_export_config(file.path(export_dir, "export.json"))
p <- as.integer(cfg$p); q <- as.integer(cfg$q)

frames <- fread(file.path(export_dir, "frames.csv"))
subjects <- fread(file.path(export_dir, "subjects.csv"))
stopifnot(nrow(frames) == cfg$n_frames)

con <- file(file.path(export_dir, "frames.f32"), "rb")
values <- readBin(con, what = "numeric", size = 4L, n = p * q * nrow(frames))
close(con)
stopifnot(length(values) == p * q * nrow(frames))
dim(values) <- c(p, q, nrow(frames))

rds_dir <- file.path(export_dir, "rds")
dir.create(rds_dir, showWarnings = FALSE)
frames[, t_idx := grid_idx + 1L]
setorder(frames, id, t_idx)
frames[, rds_file := file.path(rds_dir, paste0(id, ".rds"))]
for (sid in subjects$id) {
  rows <- frames[id == sid]
  x <- values[, , rows$offset + 1L, drop = FALSE]
  saveRDS(rTensor::as.tensor(x), rows$rds_file[[1L]], compress = FALSE)
}
fwrite(frames, file.path(export_dir, "frames_rds.csv"))
cat(sprintf("wrote %d subject tensors (%d x %d, %d frames) to %s\n",
            nrow(subjects), p, q, nrow(frames), rds_dir))
