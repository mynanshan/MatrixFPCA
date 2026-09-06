#!/usr/bin/env Rscript
# Joint-score tests and an ADNI-style principal curve through RETFound scores.
args <- commandArgs(trailingOnly = TRUE)
model_path <- if (length(args)) args[[1]] else "application/olives/output/retfound/matfpca/model.rds"
output <- if (length(args) > 1) args[[2]] else dirname(model_path)
demographics_path <- if (length(args) > 2) args[[3]] else file.path(dirname(dirname(model_path)), "demographics.csv")
biomarkers_path <- if (length(args) > 3) args[[4]] else file.path(dirname(dirname(model_path)), "visit_biomarkers.csv")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

library(data.table)
source(file.path("application", "adni", "adni_continuum_utils.R"))
fit <- readRDS(model_path)
meta <- as.data.table(fit$metadata)
ids <- as.character(fit$patient_ids)
A <- fit$scores@data
stopifnot(length(dim(A)) == 4L, dim(A)[1] == length(ids))
n <- length(ids)

# The Gram construction gives an exact Euclidean coordinate for the complete
# retained MatFPCA score vector, up to an orthogonal rotation. It avoids making
# a 600,000-dimensional principal-curve smoother directly.
S <- matrix(A, nrow = n)
S <- scale(S, center = TRUE, scale = FALSE)
G <- tcrossprod(S)
eg <- eigen((G + t(G)) / 2, symmetric = TRUE)
keep <- which(eg$values > max(eg$values) * 1e-10)
joint_pc <- sweep(eg$vectors[, keep, drop = FALSE], 2L,
                  sqrt(eg$values[keep]), "*")
joint_var <- eg$values[keep]
joint_pve <- joint_var / sum(joint_var)
joint_rank <- min(10L, ncol(joint_pc))

demog <- fread(demographics_path)
demog[, patient_id := as.character(patient_id)]
demog[, age_at_enrollment := as.numeric(age_at_enrollment)]
demographic_names <- setdiff(names(demog), c("patient_id", "study_id", "age_at_enrollment"))
biomarkers <- fread(biomarkers_path)
biomarkers[, patient_id := as.character(patient_id)]
biomarker_names <- setdiff(names(biomarkers),
                           c("patient_id", "eye_id", "arm", "study_id", "visit", "bcva", "cst"))
biomarker_base <- biomarkers[visit == 1L, c("patient_id", biomarker_names), with = FALSE]
biomarker_base <- biomarker_base[!duplicated(patient_id)]
base <- meta[order(visit), .SD[1L], by = patient_id]
end <- meta[visit == 13L, .SD[1L], by = patient_id]
clinical <- merge(base[, .(patient_id = as.character(patient_id), arm,
                            baseline_bcva = as.numeric(bcva), baseline_cst = as.numeric(cst))],
                   end[, .(patient_id = as.character(patient_id), end_bcva = as.numeric(bcva),
                           end_cst = as.numeric(cst))], by = "patient_id", all.x = TRUE)
clinical[, `:=`(bcva_change = end_bcva - baseline_bcva,
                cst_change = end_cst - baseline_cst)]
clinical <- merge(data.table(patient_id = ids), clinical, by = "patient_id", all.x = TRUE)
clinical <- merge(clinical, demog[, .(patient_id, age_at_enrollment)], by = "patient_id", all.x = TRUE)
clinical <- merge(clinical, demog[, c("patient_id", demographic_names), with = FALSE],
                  by = "patient_id", all.x = TRUE)
clinical <- merge(clinical, biomarker_base, by = "patient_id", all.x = TRUE)
clinical <- clinical[match(ids, patient_id)]

joint_test <- function(y, arm, X, permutations = 4999L, seed = 20260905L) {
  ok <- is.finite(y) & apply(X, 1L, function(z) all(is.finite(z)))
  y <- y[ok]; arm <- factor(arm[ok]); X <- X[ok, , drop = FALSE]
  if (length(y) < 10L || nlevels(arm) < 2L) {
    return(c(n = length(y), rank = ncol(X), r2 = NA, p = NA))
  }
  reduced <- lm(y ~ arm)
  full <- lm(y ~ arm + X)
  rss0 <- sum(residuals(reduced)^2); rss1 <- sum(residuals(full)^2)
  observed <- max(0, 1 - rss1 / rss0)
  fitted0 <- fitted(reduced); residual0 <- residuals(reduced)
  set.seed(seed)
  permuted <- numeric(permutations)
  groups <- split(seq_along(y), arm)
  for (b in seq_len(permutations)) {
    rp <- residual0
    for (g in groups) rp[g] <- sample(residual0[g], length(g), replace = FALSE)
    stat <- lm(fitted0 + rp ~ arm + X)
    permuted[b] <- max(0, 1 - sum(residuals(stat)^2) / rss0)
  }
  c(n = length(y), rank = ncol(X), r2 = observed,
    p = (1 + sum(permuted >= observed)) / (permutations + 1))
}

# Distance/kernel association test using every nonzero Gram coordinate. This
# avoids treating the first ten coordinates as the only representation of the
# score vector in the primary joint screen.
full_geometry_test <- function(y, arm, X, permutations = 4999L, seed = 20260907L) {
  ok <- is.finite(y)
  y <- y[ok]; arm <- factor(arm[ok]); X <- X[ok, , drop = FALSE]
  if (length(y) < 10L || nlevels(arm) < 2L) return(c(stat = NA, p = NA, n = length(y)))
  residual <- residuals(lm(y ~ arm))
  H <- diag(length(y)) - 1 / length(y)
  K <- -0.5 * H %*% (as.matrix(dist(X))^2) %*% H
  statistic <- sum(K * tcrossprod(residual)) /
    sqrt(sum(K^2) * sum(residual^2)^2)
  set.seed(seed)
  permuted <- numeric(permutations)
  groups <- split(seq_along(y), arm)
  for (b in seq_len(permutations)) {
    rp <- residual
    for (g in groups) rp[g] <- sample(residual[g], length(g), replace = FALSE)
    permuted[b] <- sum(K * tcrossprod(rp)) /
      sqrt(sum(K^2) * sum(rp^2)^2)
  }
  c(stat = statistic, p = (1 + sum(permuted >= statistic)) / (permutations + 1), n = length(y))
}

outcomes <- list(age = clinical$age_at_enrollment,
                 baseline_bcva = clinical$baseline_bcva,
                 baseline_cst = clinical$baseline_cst,
                 v13_bcva_change = clinical$bcva_change,
                 v13_cst_change = clinical$cst_change)
outcomes <- c(outcomes, setNames(lapply(biomarker_names, function(name) clinical[[name]]),
                                 paste0("baseline biomarker: ", biomarker_names)))
demographic_numeric <- list()
for (name in demographic_names) {
  value <- suppressWarnings(as.numeric(clinical[[name]]))
  if (sum(is.finite(value)) >= 10L && length(unique(value[is.finite(value)])) > 1L) {
    demographic_numeric[[paste0("demographic: ", name)]] <- value
  }
}
outcomes <- c(outcomes, demographic_numeric)
joint_rows <- rbindlist(lapply(names(outcomes), function(name) {
  z <- joint_test(outcomes[[name]], clinical$arm, joint_pc[, seq_len(joint_rank), drop = FALSE])
  g <- full_geometry_test(outcomes[[name]], clinical$arm, joint_pc)
  data.table(outcome = name, n = z[["n"]], rank = z[["rank"]], joint_r2 = z[["r2"]],
             permutation_p = z[["p"]], full_geometry_stat = g[["stat"]],
             full_geometry_p = g[["p"]])
}))
joint_rows[, permutation_q := p.adjust(permutation_p, method = "BH")]
joint_rows[, full_geometry_q := p.adjust(full_geometry_p, method = "BH")]
fwrite(joint_rows, file.path(output, "joint_score_omnibus_tests.csv"))

curve <- fit_score_principal_curve(joint_pc, spline_df = 5L, max_iterations = 20L)
curve <- orient_principal_curve(curve, clinical$age_at_enrollment)
coordinate <- curve$coordinate
curve_table <- data.table(patient_id = ids, morphology_coordinate = coordinate,
                          curve_residual = curve$residual_distance,
                          age_at_enrollment = clinical$age_at_enrollment,
                          arm = clinical$arm, baseline_bcva = clinical$baseline_bcva,
                          baseline_cst = clinical$baseline_cst,
                          v13_bcva_change = clinical$bcva_change,
                          v13_cst_change = clinical$cst_change)
for (name in biomarker_names) curve_table[[name]] <- clinical[[name]]
for (name in demographic_names) curve_table[[name]] <- clinical[[name]]
fwrite(curve_table, file.path(output, "morphology_coordinates_clinical.csv"))
curve_points <- as.data.table(curve$curve_scores)
setnames(curve_points, paste0("joint_pc_", seq_len(ncol(curve_points))))
curve_points[, curve_coordinate := curve$curve_coordinate]
fwrite(curve_points, file.path(output, "principal_curve_scores.csv"))

curve_associations <- rbindlist(lapply(names(outcomes), function(name) {
  y <- outcomes[[name]]; ok <- is.finite(y) & is.finite(coordinate)
  if (sum(ok) < 5L || length(unique(y[ok])) < 2L) {
    return(data.table(outcome = name, n = sum(ok), spearman_rho = NA_real_,
                      spearman_p = NA_real_, arm_adjusted_r = NA_real_,
                      arm_adjusted_p = NA_real_))
  }
  rho <- suppressWarnings(cor.test(coordinate[ok], y[ok], method = "spearman"))
  z <- joint_test(y, clinical$arm, matrix(coordinate, ncol = 1L), permutations = 4999L, seed = 20260906L)
  data.table(outcome = name, n = sum(ok), spearman_rho = unname(rho$estimate),
             spearman_p = rho$p.value, arm_adjusted_r = sqrt(z[["r2"]]) * sign(cor(coordinate[ok], y[ok])),
             arm_adjusted_p = z[["p"]])
}))
curve_associations[, spearman_q := p.adjust(spearman_p, method = "BH")]
curve_associations[, arm_adjusted_q := p.adjust(arm_adjusted_p, method = "BH")]
fwrite(curve_associations, file.path(output, "coordinate_associations.csv"))
fwrite(data.table(n_patients = n, score_vector_rank = ncol(joint_pc), joint_test_rank = joint_rank,
                  score_variance_captured = sum(joint_pve), curve_variance_explained = curve$variance_explained,
                  linear_variance_explained = curve$linear_variance_explained,
                  curve_iterations = curve$iterations,
                  age_orientation_spearman = curve$orientation_reference_spearman),
       file.path(output, "curve_diagnostics.csv"))
