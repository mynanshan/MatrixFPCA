#!/usr/bin/env Rscript
# How many TRAIN subjects does MatFPCA need? Compares fits on nested
# stratified samples (and one disjoint sample) with a reference fit.
#   Rscript application/echonet/compare_sample_sizes.R <ref_dir> <out_csv> <label>=<dir> ...
# For each fit: run time, temporal rank, dynamic share, agreement with the
# reference (temporal eigenfunctions, row/column eigenspaces, mean video,
# TEST score subspace) and the downstream results (EF TEST R^2 with k chosen
# on VAL; held-out frame RMSE). All fits must share VAL/TEST subjects.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
ref_dir <- a[1]; out_csv <- a[2]
runs <- do.call(rbind, strsplit(a[-(1:2)], "=", fixed = TRUE))
ref <- readRDS(file.path(ref_dir, "model.rds"))
h <- diff(ref$t_grid)[1]
tcos <- function(f, g) abs(sum(f * g) * h) / sqrt(sum(f^2) * h * sum(g^2) * h)
pcos <- function(A, B, r) mean(svd(crossprod(A[, 1:r], B[, 1:r]))$d)
ref_sc <- fread(file.path(ref_dir, "scores.csv"))[design == "dense" & split == "TEST"]
fpc <- sprintf("FPC%03d", 1:50)
rows <- lapply(seq_len(nrow(runs)), function(i) {
  d <- runs[i, 2]; f <- readRDS(file.path(d, "model.rds"))
  td <- fread(file.path(d, "temporal_diagnostics.csv"))
  sel <- fread(file.path(d, "prediction_selected_k.csv"))[target == "EF"]
  rec <- fread(file.path(d, "reconstruction_summary.csv"))
  sc <- fread(file.path(d, "scores.csv"))[design == "dense" & split == "TEST"]
  sc <- sc[match(ref_sc$id, id)]
  X <- as.matrix(sc[, ..fpc])
  r2 <- vapply(fpc, function(v) summary(lm(ref_sc[[v]] ~ X))$r.squared, 1)
  Th <- f$eig$funcs$Theta; Tr <- ref$eig$funcs$Theta
  data.table(
    run = runs[i, 1],
    n_fit = if (!is.null(f$args$n_fit) && f$args$n_fit > 0) f$args$n_fit else
      uniqueN(fread(file.path(d, "scores.csv"), select = c("id", "split"))[split == "TRAIN", id]),
    n_regression = uniqueN(fread(file.path(d, "scores.csv"), select = c("id", "split"))[split == "TRAIN", id]),
    fit_seconds = f$elapsed,
    temporal_rank = ncol(Th),
    dynamic_share = td[dynamic_screen == TRUE, sum(marginal_variance_fraction)],
    theta2_cos = tcos(Th[, 2], Tr[, 2]),
    theta3_cos = tcos(Th[, 3], Tr[, 3]),
    row_space10_cos = pcos(f$eig$funcs$Phi, ref$eig$funcs$Phi, 10),
    col_space10_cos = pcos(f$eig$funcs$Psi, ref$eig$funcs$Psi, 10),
    mean_rel_diff = sqrt(sum((f$mu - ref$mu)^2) / sum(ref$mu^2)),
    test_score_R2_mean = mean(r2), test_score_R2_min = min(r2),
    EF_R2_matfpca_sparse = sel[method == "MatFPCA predicted video -> PCA (sparse frames)", R2],
    EF_k_matfpca_sparse = sel[method == "MatFPCA predicted video -> PCA (sparse frames)", k],
    EF_R2_interp_sparse = sel[method == "PCA: Stacked sparse video (linear interpolation, 11 phases)", R2],
    EF_R2_dense_video = sel[method == "PCA: Stacked dense video (11 phases)", R2],
    heldout_rmse_matfpca = rec[frames == "held-out" & method == "MatFPCA", rmse])
})
out <- rbindlist(rows)
fwrite(out, out_csv)
print(out[, lapply(.SD, function(v) if (is.numeric(v)) signif(v, 3) else v)])
