#!/usr/bin/env Rscript
# PCA modes of the complete MatFPCA score vector versus clinical quantities.
args <- commandArgs(trailingOnly = TRUE)
model_path <- if (length(args)) args[[1]] else "application/olives/output/retfound/matfpca/model.rds"
output <- if (length(args) > 1) args[[2]] else dirname(model_path)
demographics_path <- if (length(args) > 2) args[[3]] else file.path(dirname(dirname(model_path)), "demographics.csv")
biomarkers_path <- if (length(args) > 3) args[[4]] else file.path(dirname(dirname(model_path)), "visit_biomarkers.csv")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
library(data.table)

fit <- readRDS(model_path)
meta <- as.data.table(fit$metadata)
ids <- as.character(fit$patient_ids)
A <- fit$scores@data
S <- scale(matrix(A, nrow = length(ids)), center = TRUE, scale = FALSE)
G <- tcrossprod(S)
eg <- eigen((G + t(G)) / 2, symmetric = TRUE)
keep <- which(eg$values > max(eg$values) * 1e-10)
pc <- sweep(eg$vectors[, keep, drop = FALSE], 2L, sqrt(eg$values[keep]), "*")
pc_var <- eg$values[keep]
pc_pve <- pc_var / sum(pc_var)

demog <- fread(demographics_path)
demog[, patient_id := as.character(patient_id)]
demog[, age_at_enrollment := as.numeric(age_at_enrollment)]
demographic_names <- setdiff(names(demog), c("patient_id", "study_id", "age_at_enrollment"))
bio <- fread(biomarkers_path)
bio[, patient_id := as.character(patient_id)]
biomarker_names <- setdiff(names(bio), c("patient_id", "eye_id", "arm", "study_id", "visit", "bcva", "cst"))
bio_base <- bio[visit == 1L, c("patient_id", biomarker_names), with = FALSE]
bio_base <- bio_base[!duplicated(patient_id)]
base <- meta[order(visit), .SD[1L], by = patient_id]
end <- meta[visit == 13L, .SD[1L], by = patient_id]
clinical <- merge(base[, .(patient_id = as.character(patient_id), arm,
                           baseline_bcva = as.numeric(bcva), baseline_cst = as.numeric(cst))],
                  end[, .(patient_id = as.character(patient_id), end_bcva = as.numeric(bcva),
                          end_cst = as.numeric(cst))], by = "patient_id", all.x = TRUE)
clinical[, c("bcva_change", "cst_change") := list(end_bcva - baseline_bcva,
                                                   end_cst - baseline_cst)]
clinical <- merge(data.table(patient_id = ids), clinical, by = "patient_id", all.x = TRUE)
clinical <- merge(clinical, demog[, c("patient_id", "age_at_enrollment", demographic_names), with = FALSE],
                  by = "patient_id", all.x = TRUE)
clinical <- merge(clinical, bio_base, by = "patient_id", all.x = TRUE)
clinical <- clinical[match(ids, patient_id)]

outcomes <- list(age = clinical$age_at_enrollment,
                 baseline_bcva = clinical$baseline_bcva,
                 baseline_cst = clinical$baseline_cst,
                 v13_bcva_change = clinical$bcva_change,
                 v13_cst_change = clinical$cst_change)
outcomes <- c(outcomes, setNames(lapply(biomarker_names, function(name) clinical[[name]]),
                                 paste0("baseline biomarker: ", biomarker_names)))
for (name in demographic_names) {
  value <- suppressWarnings(as.numeric(clinical[[name]]))
  if (sum(is.finite(value)) >= 10L && length(unique(value[is.finite(value)])) > 1L) {
    outcomes[[paste0("demographic: ", name)]] <- value
  }
}

# A one-positive binary label can produce a small residual-correlation p-value
# through leverage alone. Keep such outcomes in the audit table, but flag them
# so that the clinically interpretable screen is not confused with a sparse
# label artifact. Binary/low-cardinality outcomes require at least five
# observations in each level; continuous outcomes only require variation.
outcome_meta <- rbindlist(lapply(names(outcomes), function(name) {
  y <- outcomes[[name]]
  vals <- y[is.finite(y)]
  counts <- table(vals)
  low_cardinality <- length(counts) <= 5L
  min_level_n <- if (length(counts)) min(counts) else 0L
  data.table(outcome = name,
             outcome_n = length(vals),
             outcome_unique = length(counts),
             outcome_min_level_n = min_level_n,
             outcome_arm_levels = length(unique(clinical$arm[is.finite(y)])),
             outcome_screenable = length(vals) >= 10L &&
               length(unique(clinical$arm[is.finite(y)])) >= 2L &&
               length(counts) > 1L &&
               (!low_cardinality || min_level_n >= 5L))
}))

adjusted <- function(x, y, arm) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 8L || length(unique(arm[ok])) < 2L) {
    return(c(n = sum(ok), r = NA, p = NA, rho = NA, rho_p = NA))
  }
  z <- model.matrix(~ arm, data.frame(arm = factor(arm[ok])))
  rx <- residuals(lm(x[ok] ~ z - 1)); ry <- residuals(lm(y[ok] ~ z - 1))
  pearson <- suppressWarnings(cor.test(rx, ry, method = "pearson"))
  spearman <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  c(n = sum(ok), r = unname(pearson$estimate), p = pearson$p.value,
    rho = unname(spearman$estimate), rho_p = spearman$p.value)
}
rows <- list(); pos <- 1L
for (j in seq_len(ncol(pc))) {
  for (name in names(outcomes)) {
    z <- adjusted(pc[, j], outcomes[[name]], clinical$arm)
    rows[[pos]] <- data.table(mode = j, outcome = name, n = z[["n"]],
                               pve = pc_pve[j], arm_adjusted_r = z[["r"]],
                               arm_adjusted_p = z[["p"]], spearman_rho = z[["rho"]],
                               spearman_p = z[["rho_p"]])
    pos <- pos + 1L
  }
}
result <- rbindlist(rows)
result <- merge(result, outcome_meta, by = "outcome", all.x = TRUE, sort = FALSE)
result[, arm_adjusted_q := p.adjust(arm_adjusted_p, method = "BH")]
result[, spearman_q := p.adjust(spearman_p, method = "BH")]
result[, arm_adjusted_q_within_outcome := ave(arm_adjusted_p, outcome,
                                              FUN = function(p) p.adjust(p, "BH"))]
# This family excludes constant or extremely sparse labels. It is a useful
# sensitivity correction, while the all-tests q-value remains the primary
# conservative correction reported in the markdown report.
screenable_p <- ifelse(result$outcome_screenable, result$arm_adjusted_p, NA_real_)
result[, arm_adjusted_q_screenable_global := p.adjust(screenable_p, method = "BH")]
result[, arm_adjusted_q_within_screenable_outcome := ifelse(
  outcome_screenable,
  ave(arm_adjusted_p, outcome, FUN = function(p) p.adjust(p, "BH")), NA_real_)]
result[, rank_by_variance := mode]
fwrite(result, file.path(output, "pca_mode_clinical_associations.csv"))
fwrite(data.table(mode = seq_len(ncol(pc)), pve = pc_pve,
                  cumulative_pve = cumsum(pc_pve)),
       file.path(output, "score_pca_spectrum.csv"))
