#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
model_path <- if (length(args)) args[[1]] else "application/olives/output/retfound/matfpca/model.rds"
output <- if (length(args) > 1) args[[2]] else dirname(model_path)
library(data.table)
fit <- readRDS(model_path)
meta <- as.data.table(fit$metadata)
ids <- fit$patient_ids
A <- fit$scores@data
stopifnot(length(dim(A)) == 4L, dim(A)[1] == length(ids))
dimnames(A) <- NULL
v <- apply(A, c(2, 3, 4), var)
keep <- order(as.vector(v), decreasing = TRUE)[seq_len(min(20L, length(v)))]
idx <- arrayInd(keep, .dim = dim(v))
score_names <- sprintf("FPC_r%03d_c%04d_t%02d", idx[,1], idx[,2], idx[,3])
base <- meta[order(visit), .SD[1L], by = patient_id]
end <- meta[visit == 13L, .SD[1L], by = patient_id]
clinical <- merge(base[, .(patient_id, arm, base_bcva = as.numeric(bcva), base_cst = as.numeric(cst))],
                   end[, .(patient_id, end_bcva = as.numeric(bcva), end_cst = as.numeric(cst))],
                   by = "patient_id", all.x = TRUE)
clinical[, `:=`(bcva_change = end_bcva - base_bcva, cst_change = end_cst - base_cst)]
clinical <- clinical[match(ids, patient_id)]
adjusted <- function(x, y, arm) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5L) return(c(r = NA_real_, p = NA_real_, n = sum(ok)))
  z <- model.matrix(~ arm, data.frame(arm = factor(arm[ok])))
  rx <- residuals(lm(x[ok] ~ z - 1)); ry <- residuals(lm(y[ok] ~ z - 1))
  ct <- suppressWarnings(cor.test(rx, ry, method = "pearson"))
  c(r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}
outcomes <- list(baseline_bcva = clinical$base_bcva, baseline_cst = clinical$base_cst,
                 v13_bcva_change = clinical$bcva_change, v13_cst_change = clinical$cst_change)
rows <- list(); pos <- 1L
for (s in seq_len(nrow(idx))) {
  score <- A[, idx[s,1], idx[s,2], idx[s,3]]
  for (outcome in names(outcomes)) {
    y <- outcomes[[outcome]]; ok <- is.finite(score) & is.finite(y)
    rho <- if (sum(ok) >= 5L) suppressWarnings(cor.test(score[ok], y[ok], method = "spearman")) else NULL
    adj <- adjusted(score, y, clinical$arm)
    rows[[pos]] <- data.table(score = score_names[s], row_patch = idx[s,1], feature = idx[s,2],
                              temporal_component = idx[s,3], outcome = outcome, n = adj[["n"]],
                              spearman_rho = if (is.null(rho)) NA_real_ else unname(rho$estimate),
                              spearman_p = if (is.null(rho)) NA_real_ else rho$p.value,
                              arm_adjusted_r = adj[["r"]], arm_adjusted_p = adj[["p"]],
                              product_variance = v[idx[s,1], idx[s,2], idx[s,3]])
    pos <- pos + 1L
  }
}
result <- rbindlist(rows)
result[, arm_adjusted_q := p.adjust(arm_adjusted_p, method = "BH")]
result[, spearman_q := p.adjust(spearman_p, method = "BH")]
fwrite(result, file.path(output, "score_clinical_associations.csv"))
fwrite(data.table(score = score_names, row_patch = idx[,1], feature = idx[,2],
                  temporal_component = idx[,3], product_variance = v[keep]),
       file.path(output, "top_product_scores.csv"))
