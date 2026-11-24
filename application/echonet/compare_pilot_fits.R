#!/usr/bin/env Rscript
# Sparse-vs-dense comparison of two MatFPCA fits (PLAN.md Section 2.3).
#   Rscript application/echonet/compare_pilot_fits.R <sparse_dir> <dense_dir> <out_csv>
# Reports, for the fits trained on sparse (3-8 frames) and on all frames:
#   * |cosine| between matched temporal eigenfunctions (L2 on the phase grid);
#   * principal-angle cosines between the leading row / column eigenspaces;
#   * relative L2 difference of the mean videos;
#   * R^2 of each dense-fit TEST score (FPC1-50) on the sparse-fit FPC1-50 scores.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
fs <- readRDS(file.path(a[1], "model.rds")); fd <- readRDS(file.path(a[2], "model.rds"))
tg <- fs$t_grid; h <- diff(tg)[1]
rows <- list()
nK <- min(4L, ncol(fs$eig$funcs$Theta), ncol(fd$eig$funcs$Theta))
for (l in seq_len(nK)) {
  cs <- abs(sum(fs$eig$funcs$Theta[, l] * fd$eig$funcs$Theta[, l]) * h) /
    sqrt(sum(fs$eig$funcs$Theta[, l]^2) * h * sum(fd$eig$funcs$Theta[, l]^2) * h)
  rows[[length(rows) + 1]] <- data.table(quantity = sprintf("temporal eigenfunction %d, |cos|", l), value = cs)
}
for (d in c("Phi", "Psi")) for (r in c(5L, 10L, 20L)) {
  cs <- svd(crossprod(fs$eig$funcs[[d]][, 1:r], fd$eig$funcs[[d]][, 1:r]))$d
  rows[[length(rows) + 1]] <- data.table(
    quantity = sprintf("%s leading %d-dim subspace, mean / min principal cosine",
                       if (d == "Phi") "row" else "column", r),
    value = mean(cs), value2 = min(cs))
}
rows[[length(rows) + 1]] <- data.table(quantity = "mean video, relative L2 difference",
                                       value = sqrt(sum((fs$mu - fd$mu)^2) / sum(fd$mu^2)))
td <- function(dir) fread(file.path(dir, "temporal_diagnostics.csv"))
rows[[length(rows) + 1]] <- data.table(quantity = "temporal variance in dynamic modes (sparse fit / dense fit)",
  value = td(a[1])[dynamic_screen == TRUE, sum(marginal_variance_fraction)],
  value2 = td(a[2])[dynamic_screen == TRUE, sum(marginal_variance_fraction)])
ss <- fread(file.path(a[1], "scores.csv"))[design == "dense" & split == "TEST"]
sd_ <- fread(file.path(a[2], "scores.csv"))[design == "dense" & split == "TEST"]
sd_ <- sd_[match(ss$id, id)]
# FPCs with near-tied eigenvalues can swap order between fits, so compare
# score subspaces: R^2 of each dense-fit score on all top-50 sparse-fit scores
Xs <- as.matrix(ss[, sprintf("FPC%03d", 1:50), with = FALSE])
r2 <- vapply(1:50, function(f) summary(lm(sd_[[sprintf("FPC%03d", f)]] ~ Xs))$r.squared, 1)
rows[[length(rows) + 1]] <- data.table(
  quantity = "TEST scores: R2 of dense-fit FPC1-50 on sparse-fit FPC1-50 (mean / min)",
  value = mean(r2), value2 = min(r2))
out <- rbindlist(rows, fill = TRUE)
fwrite(out, a[3]); print(out)
