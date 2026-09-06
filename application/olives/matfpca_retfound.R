#!/usr/bin/env Rscript
# Fit MatFPCA to frozen RETFound patch-token matrices.
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[[1]] else "application/olives/output/retfound/features"
output <- if (length(args) > 1) args[[2]] else file.path(input, "matfpca")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

library(data.table)
library(rTensor)
source(file.path("R", "matfdata.R"))
source(file.path("R", "matfpca.R"))

manifest <- fread(file.path(dirname(input), "manifest.csv"))
stopifnot(nrow(manifest) > 0L)
meta <- manifest[order(patient_id, visit)]
n <- nrow(meta); P <- 196L; D <- 1024L
con <- file(file.path(input, "patch_tokens.bin"), "rb")
raw <- readBin(con, "numeric", n = n * P * D, size = 4L, endian = "little")
close(con)
stopifnot(length(raw) == n * P * D, all(is.finite(raw)))
# Binary order is observation, patch, feature; R's first dimension varies
# fastest, so this is feature x patch x observation before transposition.
Y <- t(matrix(raw, nrow = P * D))
rm(raw)

meta[, arm := factor(arm)]
max_visit <- max(meta$visit)
meta[, time := (visit - 1) / max(1, max_visit - 1)]
basis <- splines::ns(meta$time, knots = c(1/3, 2/3), Boundary.knots = c(0, 1))
design <- model.matrix(~ arm * basis, meta)
patient_tab <- table(meta$patient_id)
weight <- 1 / as.numeric(patient_tab[as.character(meta$patient_id)])
B <- solve(crossprod(design, design * weight) + diag(1e-8, ncol(design)),
           crossprod(design, Y * weight))
residual <- Y - design %*% B
patient_ids <- sort(unique(meta$patient_id))
X_list <- lapply(patient_ids, function(id) {
  ix <- which(meta$patient_id == id)
  rTensor::as.tensor(array(t(residual[ix, , drop = FALSE]), dim = c(P, D, length(ix))))
})
t_grid <- seq(0, 1, length.out = max_visit)
t_idx <- lapply(patient_ids, function(id) as.integer(meta$visit[meta$patient_id == id]))
mat_data <- create_matrix_funcdata_irreg(X_list = X_list, t_grid = t_grid, t_idx_list = t_idx)

message("Fitting RETFound patch x feature MatFPCA: ", n, " visits, ", length(patient_ids), " patients")
timing <- system.time(covariance_fit <- fit_matfcov(
  mat_data, wScheme = "identity", noise = TRUE, nbasis = "theory",
  opwIter = 1L, temporal_method = "two_step", crossfit_diagonal = "identity"
))
eig_full <- eigen.matfcov(covariance_fit$C)
eig <- trunc.matfeigen(eig_full, pve = 0.95, pve_rule = "temporal", verbose = TRUE)
eig <- trunc.matfeigen(eig_full, I = P, J = D, K = max(2L, ncol(eig$funcs$Theta)), verbose = FALSE)
scores <- get_matfpca_scores(mat_data, eig, sig2 = covariance_fit$C$sig2,
                             K_grid = covariance_fit$C$K)
saveRDS(list(covariance_fit = covariance_fit, eig = eig, scores = scores,
             metadata = meta, patient_ids = patient_ids,
             dimensions = c(patch = P, feature = D, time = max_visit),
             elapsed_seconds = unname(timing[["elapsed"]])),
        file.path(output, "model.rds"), compress = FALSE)
score_matrix <- matrix(scores@data, nrow = length(patient_ids))
score_names <- as.data.table(as.data.frame(score_matrix))
setnames(score_names, paste0("score_", seq_len(ncol(score_names))))
score_names[, patient_id := patient_ids]
fwrite(score_names, file.path(output, "scores.csv"))
fwrite(data.table(component = seq_along(eig$values$lam_theta),
                  temporal_eigenvalue = eig$values$lam_theta),
       file.path(output, "temporal_spectrum.csv"))
fwrite(data.table(n_visits = n, n_patients = length(patient_ids), patch = P,
                  feature = D, temporal_rank = ncol(eig$funcs$Theta),
                  noise_variance = covariance_fit$C$sig2,
                  elapsed_seconds = unname(timing[["elapsed"]])),
       file.path(output, "run_summary.csv"))
