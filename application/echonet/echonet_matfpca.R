#!/usr/bin/env Rscript
# EchoNet-Dynamic MatFPCA fit and evaluation (PLAN.md Sections 2-4).
#
#   Rscript application/echonet/echonet_matfpca.R --input=<export_dir> \
#     --output=<dir> [--design=sparse] [--m-min=3] [--m-max=8] ...
#
# Input: an export directory after echonet_to_rds.R (per-subject rTensor RDS
# files of gain-normalised, masked systolic frames on the 51-point phase grid).
#
# Steps
#   1. sampling design: TRAIN subjects keep all frames (dense) or a random
#      m_i ~ U{m_min..m_max} subset (sparse);
#   2. mean mu(t): natural spline in phase, subject-balanced, TRAIN only;
#   3. MatFPCA: fit_matfcov -> eigen.matfcov -> trunc.matfeigen (same settings
#      as application/adni/adni_matfpca.R);
#   4. temporal diagnostics (ADNI constant-vs-dynamic screen);
#   5. scores (BLUP) for TRAIN/VAL/TEST from sparse and from dense frames;
#   6. EF prediction: ridge on the leading product-FPC scores, penalty tuned
#      on VAL, evaluated on TEST; PCA baselines on the same subjects;
#   7. reconstruction of every TEST frame from the sparse frames;
#   8. figures.

suppressPackageStartupMessages({
  library(data.table)
  library(rTensor)
})

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]))
} else normalizePath(file.path("application", "echonet", "echonet_matfpca.R"))
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
setwd(project_root)
source(file.path("R", "matfpca.R"))

# ---------------------------------------------------------------- options ----
defaults <- list(
  input = "application/echonet/output/export_pilot56",
  output = "application/echonet/output/pilot56_sparse",
  design = "sparse",
  m_min = 3L, m_max = 8L,
  mean_df = 6L,
  weight_scheme = "gaussian_refined",
  weight_folds = 3L,
  opw_iter = 1L,
  nbasis = "theory",
  mode_pve = 0.95,
  n_fpcs = "10,20,50,100,200,500",
  seed = 20261124L,
  cache = TRUE
)
parse_options <- function(d) {
  for (arg in commandArgs(trailingOnly = TRUE)) {
    if (!grepl("^--[^=]+=", arg)) stop("Use --name=value: ", arg)
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", kv[[1L]], fixed = TRUE)
    if (!key %in% names(d)) stop("Unknown option --", kv[[1L]])
    val <- paste(kv[-1L], collapse = "=")
    d[[key]] <- if (is.logical(d[[key]])) tolower(val) %in% c("true", "1", "yes")
      else if (is.integer(d[[key]])) as.integer(val)
      else if (is.numeric(d[[key]])) as.numeric(val) else val
  }
  d
}
args <- parse_options(defaults)
stopifnot(args$design %in% c("sparse", "dense"))
n_fpcs <- as.integer(strsplit(args$n_fpcs, ",")[[1L]])
dir.create(args$output, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(args$output, "figures")
dir.create(fig_dir, showWarnings = FALSE)
log_message <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = ""); flush.console()
}
writeLines(capture.output(str(args)), file.path(args$output, "options.txt"))

# ------------------------------------------------------------------- data ----
frames <- fread(file.path(args$input, "frames_rds.csv"))
subjects <- fread(file.path(args$input, "subjects.csv"))
mask <- as.matrix(fread(file.path(args$input, "mask.csv"), header = FALSE)) == 1
p <- nrow(mask); q <- ncol(mask)
t_grid <- seq(0, 1, by = 0.02)
nt <- length(t_grid)
setorder(subjects, split, id)
subjects[, split := factor(split, c("TRAIN", "VAL", "TEST"))]
setorder(subjects, split, id)
log_message("Subjects: ", paste(names(table(subjects$split)), table(subjects$split),
                                 collapse = ", "), "; image ", p, " x ", q,
            "; mask ", sum(mask), " px")

log_message("Loading subject tensors...")
raw <- lapply(subjects$id, function(s) readRDS(file.path(args$input, "rds", paste0(s, ".rds")))@data)
names(raw) <- subjects$id
t_idx_all <- split(frames$t_idx, factor(frames$id, levels = subjects$id))
stopifnot(all(vapply(raw, function(x) dim(x)[3], 1L) == lengths(t_idx_all)))

# sparse subset of frames per subject (TRAIN, VAL and TEST all get one, so the
# sparse-design evaluation is identical across splits)
set.seed(args$seed)
sparse_pos <- lapply(t_idx_all, function(ti) {
  n <- length(ti)
  m <- min(n, sample(args$m_min:args$m_max, 1L))
  sort(sample.int(n, m))
})

# ---------------------------------------------------------- mean function ----
fit_mean_function <- function(times, values, subject_ids, df) {
  basis <- splines::ns(times, df = df - 1L, intercept = FALSE, Boundary.knots = c(0, 1))
  design <- cbind(1, basis)
  w <- 1 / as.numeric(table(subject_ids)[subject_ids])
  nm <- crossprod(design, design * w)
  coef <- solve(nm + diag(1e-8 * mean(diag(nm)), ncol(design)), crossprod(design, values * w))
  list(coef = coef, knots = attr(basis, "knots"))
}
mean_on_grid <- function(model, times) {
  cbind(1, splines::ns(times, knots = model$knots, intercept = FALSE,
                       Boundary.knots = c(0, 1))) %*% model$coef
}
train_ids <- subjects[split == "TRAIN", id]
train_pos <- if (args$design == "sparse") sparse_pos[train_ids] else
  lapply(t_idx_all[train_ids], seq_along)
mean_times <- unlist(Map(function(ti, pos) t_grid[ti[pos]], t_idx_all[train_ids], train_pos))
mean_vals <- do.call(rbind, Map(function(x, pos) t(matrix(x[, , pos, drop = FALSE], p * q)),
                                raw[train_ids], train_pos))
mean_subj <- rep(train_ids, lengths(train_pos))
mean_model <- fit_mean_function(mean_times, mean_vals, mean_subj, args$mean_df)
mu <- mean_on_grid(mean_model, t_grid)            # nt x (p*q)
mu[, !as.vector(mask)] <- 0
mu_arr <- array(t(mu), c(p, q, nt))
rm(mean_vals)
log_message("Mean fitted on ", length(mean_times), " TRAIN frames (",
            args$design, " design)")

centre <- function(x, ti) x - mu_arr[, , ti, drop = FALSE]
make_data <- function(ids, which = c("sparse", "dense", "train")) {
  which <- match.arg(which)
  pos <- switch(which,
    sparse = sparse_pos[ids],
    dense = lapply(t_idx_all[ids], seq_along),
    train = train_pos[ids])
  X <- Map(function(x, ti, ps) rTensor::as.tensor(centre(x[, , ps, drop = FALSE], ti[ps])),
           raw[ids], t_idx_all[ids], pos)
  create_matrix_funcdata_irreg(X_list = unname(X), t_grid = t_grid,
                               t_idx_list = unname(Map(function(ti, ps) as.integer(ti[ps]),
                                                       t_idx_all[ids], pos)))
}

# ------------------------------------------------------------------- fit -----
fit_file <- file.path(args$output, "model.rds")
if (args$cache && file.exists(fit_file)) {
  fit <- readRDS(fit_file)
  log_message("Reusing fit: ", fit_file)
} else {
  train_data <- make_data(train_ids, "train")
  set.seed(args$seed + 2L)
  fold_id <- sample(rep(seq_len(args$weight_folds), length.out = length(train_ids)))
  nb <- if (grepl("^[0-9]+$", args$nbasis)) as.integer(args$nbasis) else args$nbasis
  log_message("Fitting MatFPCA on ", length(train_ids), " subjects, ",
              sum(lengths(train_pos)), " frames...")
  timing <- system.time({
    cov_fit <- fit_matfcov(train_data, wScheme = args$weight_scheme, noise = TRUE,
                           nbasis = nb, opwIter = args$opw_iter,
                           temporal_method = "two_step", nfolds = args$weight_folds,
                           fold_id = fold_id, crossfit_diagonal = "refined")
  })
  eig_full <- eigen.matfcov(cov_fit$C)
  eig <- trunc.matfeigen(eig_full, pve = args$mode_pve, pve_rule = "temporal",
                         verbose = FALSE)
  if (ncol(eig$funcs$Theta) < 2L) {
    eig <- trunc.matfeigen(eig_full, I = ncol(eig_full$funcs$Phi),
                           J = ncol(eig_full$funcs$Psi), K = 2L, verbose = FALSE)
  }
  fit <- list(C = cov_fit$C, eig = eig, eig_full = eig_full, mean_model = mean_model,
              mu = mu, mask = mask, t_grid = t_grid, args = args,
              elapsed = unname(timing[["elapsed"]]))
  saveRDS(fit, fit_file)
  log_message("Fit done in ", round(fit$elapsed), " s")
}
eig <- fit$eig
C <- fit$C
ranks <- c(I = ncol(eig$funcs$Phi), J = ncol(eig$funcs$Psi), K = ncol(eig$funcs$Theta))
log_message("Ranks ", paste(ranks, collapse = " x "), "; sigma^2 = ", signif(C$sig2, 3))

# --------------------------------------------------- temporal diagnostics ----
# Same screen as ADNI: eigen-decomposition of the temporal covariance with
# trapezoid weights; a mode is 'dynamic' if it carries >= 0.1% of the temporal
# variance and >= 50% of its energy is outside its best constant.
tw <- c(diff(t_grid)[1] / 2, (t_grid[-(1:2)] - t_grid[1:(nt - 2)]) / 2, diff(t_grid)[nt - 1] / 2)
sw <- sqrt(tw)
ek <- eigen(sym(sweep(sweep(C$K, 1, sw, "*"), 2, sw, "*")), symmetric = TRUE)
tvals <- pmax(ek$values, 0)
tfun <- sweep(ek$vectors, 1, sw, "/")
n_tmodes <- min(10L, sum(tvals > max(tvals) * 1e-10))
temporal_diag <- rbindlist(lapply(seq_len(n_tmodes), function(k) {
  f <- tfun[, k]
  wm <- sum(tw * f) / sum(tw)
  data.table(time_component = k, eigenvalue = tvals[k],
             marginal_variance_fraction = tvals[k] / sum(tvals),
             nonconstant_fraction = sum(tw * (f - wm)^2) / sum(tw * f^2))
}))
temporal_diag[, cumulative := cumsum(marginal_variance_fraction)]
temporal_diag[, dynamic_screen := marginal_variance_fraction >= 0.001 & nonconstant_fraction >= 0.5]
fwrite(temporal_diag, file.path(args$output, "temporal_diagnostics.csv"))
dyn_share <- temporal_diag[dynamic_screen == TRUE, sum(marginal_variance_fraction)]
log_message("Temporal variance in dynamic modes: ", sprintf("%.1f%%", 100 * dyn_share))

# ---------------------------------------------------------------- scores -----
# Vectorised version of get_matfpca_scores() (ridge inverse, sig2 > 0); it is
# checked against the package function on a few subjects below.
lam_spatial <- outer(eig$values$lam_phi, eig$values$lam_psi)      # I x J
K_grid <- sym(C$K)
score_one <- function(x, ti) {
  proj <- rTensor::ttl(rTensor::as.tensor(x), list(t(eig$funcs$Phi), t(eig$funcs$Psi)),
                       ms = c(1, 2))@data
  m <- length(ti)
  Z <- matrix(proj, ranks[["I"]] * ranks[["J"]], m)                 # IJ x m
  e <- eigen(K_grid[ti, ti, drop = FALSE], symmetric = TRUE)
  d <- pmax(e$values, 0)
  ridge <- pmax(1e-8, C$sig2 / as.vector(lam_spatial))
  ridge[!(is.finite(lam_spatial) & lam_spatial > 0)] <- Inf
  W <- (Z %*% e$vectors) / outer(ridge, d, "+")                     # IJ x m
  num <- sweep(eig$funcs$Theta[ti, , drop = FALSE], 2, eig$values$lam_theta, "*")
  S <- W %*% t(e$vectors) %*% num                                   # IJ x K
  array(S, c(ranks[["I"]], ranks[["J"]], ranks[["K"]]))
}
score_set <- function(ids, which) {
  pos <- if (which == "sparse") sparse_pos[ids] else lapply(t_idx_all[ids], seq_along)
  out <- array(NA_real_, c(length(ids), ranks))
  for (i in seq_along(ids)) {
    ps <- pos[[i]]; ti <- t_idx_all[[ids[i]]][ps]
    out[i, , , ] <- score_one(centre(raw[[ids[i]]][, , ps, drop = FALSE], ti), ti)
  }
  out
}
# check against the package implementation
chk_ids <- subjects$id[1:3]
ref <- get_matfpca_scores(make_data(chk_ids, "sparse"), eig, sig2 = C$sig2, K_grid = C$K)@data
mine <- score_set(chk_ids, "sparse")
score_check <- max(abs(ref - mine)) / max(abs(ref))
log_message("Score implementation check: max relative difference ", signif(score_check, 3))
stopifnot(score_check < 1e-6)

# product FPCs ordered by eigenvalue
pv <- eig$values$lam_phi %o% eig$values$lam_psi %o% eig$values$lam_theta
prod_idx <- as.data.table(arrayInd(seq_along(pv), dim(pv)))
setnames(prod_idx, c("j", "k", "l"))
prod_idx[, eigenvalue := as.vector(pv)]
prod_idx <- prod_idx[eigenvalue > 0][order(-eigenvalue)]
prod_idx[, `:=`(fpc = .I, pve = eigenvalue / sum(eigenvalue))]
prod_idx[, cum_pve := cumsum(pve)]
prod_col <- with(prod_idx, j + (k - 1L) * ranks[["I"]] + (l - 1L) * ranks[["I"]] * ranks[["J"]])
fwrite(prod_idx[1:200], file.path(args$output, "product_fpcs.csv"))

log_message("Scoring subjects (sparse and dense frames)...")
ids_by <- split(subjects$id, subjects$split)
scores <- list(); scores_full <- list()
for (w in c("sparse", "dense")) for (s in names(ids_by)) {
  a <- score_set(ids_by[[s]], w)
  # all product-FPC scores, ordered by eigenvalue (= coordinates of the
  # predicted 51-phase video in an orthonormal basis)
  scores_full[[w]][[s]] <- matrix(a, length(ids_by[[s]]))[, prod_col, drop = FALSE]
  scores[[w]][[s]] <- scores_full[[w]][[s]][, 1:max(n_fpcs), drop = FALSE]
}
sc_out <- rbindlist(lapply(c("sparse", "dense"), function(w) rbindlist(lapply(names(ids_by), function(s) {
  data.table(design = w, id = ids_by[[s]], split = s,
             as.data.table(`colnames<-`(scores[[w]][[s]], sprintf("FPC%03d", seq_len(max(n_fpcs))))))
}))))
sc_out <- merge(sc_out, subjects[, .(id, EF, ESV, EDV)], by = "id")
fwrite(sc_out, file.path(args$output, "scores.csv"))

# ------------------------------------------------------- EF prediction -------
auc <- function(y, s) {
  r <- rank(c(s[y], s[!y]))
  (sum(r[seq_len(sum(y))]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(!y))
}
ridge_eval <- function(Ztr, ytr, Zva, yva, Zte, yte) {
  m <- colMeans(Ztr); s <- apply(Ztr, 2, sd); s[s == 0] <- 1
  f <- function(Z) cbind(1, sweep(sweep(Z, 2, m), 2, s, "/"))
  A <- f(Ztr); best <- NULL
  for (lam in 10^seq(-2, 4, by = 0.5)) {
    b <- solve(crossprod(A) + diag(c(0, rep(lam, ncol(Ztr)))), crossprod(A, ytr))
    rm <- sqrt(mean((f(Zva) %*% b - yva)^2))
    if (is.null(best) || rm < best$rm) best <- list(rm = rm, lam = lam, b = b)
  }
  pr <- as.vector(f(Zte) %*% best$b)
  list(R2 = 1 - sum((yte - pr)^2) / sum((yte - mean(yte))^2), MAE = mean(abs(yte - pr)),
       AUC_EF40 = auc(yte < 40, -pr), lambda = best$lam, pred = pr)
}
pca_scores <- function(Ftr, Flist, k) {
  m <- colMeans(Ftr)
  k <- min(k, nrow(Ftr) - 1L)
  sv <- svd(sweep(Ftr, 2, m), nu = 0, nv = k)
  lapply(Flist, function(F) sweep(F, 2, m) %*% sv$v)
}
y <- lapply(ids_by, function(ids) subjects[match(ids, id), EF])

# image features for baselines (in-mask pixels, gain-normalised, not centred)
mvec <- as.vector(mask)
frame_mat <- function(x) t(matrix(x, p * q))[, mvec, drop = FALSE]   # m x P
interp11 <- function(x, ti) {
  ph <- t_grid[ti]; X <- frame_mat(x)
  apply(X, 2, function(col) approx(ph, col, xout = seq(0, 1, 0.1), rule = 2)$y)
}
feat <- function(ids, kind) {
  t(vapply(ids, function(s) {
    x <- raw[[s]]; ti <- t_idx_all[[s]]
    switch(kind,
      ed = frame_mat(x[, , 1, drop = FALSE])[1, ],
      diff = { X <- frame_mat(x); X[nrow(X), ] - X[1, ] },
      avg_dense = colMeans(frame_mat(x)),
      avg_sparse = colMeans(frame_mat(x[, , sparse_pos[[s]], drop = FALSE])),
      stacked = as.vector(interp11(x, ti)),
      stacked_sparse = { ps <- sparse_pos[[s]]; as.vector(interp11(x[, , ps, drop = FALSE], ti[ps])) })
  }, numeric(if (grepl("stacked", kind)) 11 * sum(mvec) else sum(mvec))))
}
baseline_kinds <- c(
  "ED frame only (static)" = "ed",
  "Average of sparse frames (time discarded)" = "avg_sparse",
  "Average of all frames (time discarded)" = "avg_dense",
  "ES - ED difference image" = "diff",
  "Stacked dense video (11 phases)" = "stacked",
  "Stacked sparse video (linear interpolation, 11 phases)" = "stacked_sparse")
log_message("EF prediction: MatFPCA scores and PCA baselines...")
ef_rows <- list(); preds <- list()
for (w in c("sparse", "dense")) for (k in n_fpcs) {
  r <- ridge_eval(scores[[w]]$TRAIN[, 1:k], y$TRAIN, scores[[w]]$VAL[, 1:k], y$VAL,
                  scores[[w]]$TEST[, 1:k], y$TEST)
  ef_rows[[length(ef_rows) + 1L]] <- data.table(
    method = sprintf("MatFPCA scores (%s frames)", w), frames = w, k = k,
    R2 = r$R2, MAE = r$MAE, AUC_EF40 = r$AUC_EF40, lambda = r$lambda)
  if (k == max(n_fpcs)) preds[[paste("MatFPCA", w)]] <- r$pred
}
# MatFPCA as a sparse-to-full video predictor, followed by ordinary PCA of
# the predicted videos (PCA in the full score space)
pca_k <- unique(pmin(n_fpcs, length(ids_by$TRAIN) - 1L))
k_show <- if (50L %in% pca_k) 50L else max(pca_k)   # k used in the scatter plot
for (w in c("sparse", "dense")) for (k in pca_k) {
  Z <- pca_scores(scores_full[[w]]$TRAIN, scores_full[[w]], k)
  r <- ridge_eval(Z$TRAIN, y$TRAIN, Z$VAL, y$VAL, Z$TEST, y$TEST)
  ef_rows[[length(ef_rows) + 1L]] <- data.table(
    method = sprintf("MatFPCA predicted video -> PCA (%s frames)", w), frames = w, k = k,
    R2 = r$R2, MAE = r$MAE, AUC_EF40 = r$AUC_EF40, lambda = r$lambda)
  if (k == k_show) preds[[paste("MatFPCA video PCA", w)]] <- r$pred
}
for (nm in names(baseline_kinds)) {
  F <- lapply(ids_by, feat, kind = baseline_kinds[[nm]])
  for (k in pca_k) {
    Z <- pca_scores(F$TRAIN, F, k)
    r <- ridge_eval(Z$TRAIN, y$TRAIN, Z$VAL, y$VAL, Z$TEST, y$TEST)
    ef_rows[[length(ef_rows) + 1L]] <- data.table(
      method = paste("PCA:", nm),
      frames = if (baseline_kinds[[nm]] == "avg_sparse") "sparse" else
        if (baseline_kinds[[nm]] %in% c("ed", "diff")) "ED/ES landmarks" else "dense",
      k = k, R2 = r$R2, MAE = r$MAE, AUC_EF40 = r$AUC_EF40, lambda = r$lambda)
    if (k == k_show) preds[[nm]] <- r$pred
  }
}
ef_table <- rbindlist(ef_rows)
fwrite(ef_table, file.path(args$output, "ef_prediction.csv"))
print(ef_table[, .(method, k, R2 = round(R2, 3), MAE = round(MAE, 2), AUC = round(AUC_EF40, 3))])

# Spearman correlation of leading scores with EF, ESV, EDV (TRAIN, sparse)
cor_tab <- rbindlist(lapply(1:10, function(f) {
  z <- scores$sparse$TRAIN[, f]; d <- subjects[match(ids_by$TRAIN, id)]
  data.table(fpc = f, row = prod_idx$j[f], col = prod_idx$k[f], time = prod_idx$l[f],
             pve = prod_idx$pve[f],
             rho_EF = cor(z, d$EF, method = "spearman"),
             rho_ESV = cor(z, d$ESV, method = "spearman"),
             rho_EDV = cor(z, d$EDV, method = "spearman"))
}))
fwrite(cor_tab, file.path(args$output, "score_correlations.csv"))

# --------------------------------------------------------- reconstruction ----
# For each TEST subject: scores from the sparse frames, then predict every
# retained frame. Compare with the mean only and with mean + static PCA
# (time-averaged training images; the subject's sparse average is projected).
log_message("Reconstructing TEST videos from sparse frames...")
Phi <- eig$funcs$Phi; Psi <- eig$funcs$Psi; Theta <- eig$funcs$Theta
static_tr <- feat(ids_by$TRAIN, "avg_dense")
mu_bar <- colMeans(mu[, mvec, drop = FALSE])
static_sv <- svd(sweep(static_tr, 2, colMeans(static_tr)), nu = 0, nv = 20)$v
static_mean <- colMeans(static_tr)
rec_rows <- list(); rec_keep <- list()
for (s in ids_by$TEST) {
  x <- raw[[s]]; ti <- t_idx_all[[s]]; ps <- sparse_pos[[s]]
  sc <- score_one(centre(x[, , ps, drop = FALSE], ti[ps]), ti[ps])
  rec_c <- rTensor::ttl(rTensor::as.tensor(sc), list(Phi, Psi, Theta[ti, , drop = FALSE]), ms = 1:3)@data
  obs <- frame_mat(x)
  mean_only <- mu[ti, mvec, drop = FALSE]
  mfp <- mean_only + frame_mat(rec_c)
  zs <- (colMeans(obs[ps, , drop = FALSE]) - static_mean) %*% static_sv
  stat <- sweep(mean_only, 2, as.vector(static_sv %*% t(zs)) + static_mean - mu_bar, "+")
  held <- !(seq_along(ti) %in% ps)
  for (nm in c("Mean only", "Mean + static PCA (20)", "MatFPCA")) {
    pr <- switch(nm, "Mean only" = mean_only, "Mean + static PCA (20)" = stat, MatFPCA = mfp)
    for (sub in c("observed", "held-out")) {
      rows <- if (sub == "observed") !held else held
      if (!any(rows)) next
      e <- pr[rows, , drop = FALSE] - obs[rows, , drop = FALSE]
      dyn <- obs[rows, , drop = FALSE] - matrix(colMeans(obs), sum(rows), ncol(obs), byrow = TRUE)
      pd <- pr[rows, , drop = FALSE] - matrix(colMeans(pr), sum(rows), ncol(pr), byrow = TRUE)
      rec_rows[[length(rec_rows) + 1L]] <- data.table(
        id = s, method = nm, frames = sub, n = sum(rows),
        rmse = sqrt(mean(e^2)),
        corr = cor(as.vector(pr[rows, ]), as.vector(obs[rows, ])),
        # how well the within-video motion (frame minus video average) is captured
        motion_corr = if (sum(rows) > 1) cor(as.vector(pd), as.vector(dyn)) else NA_real_)
    }
  }
  if (length(rec_keep) < 3L && s %in% ids_by$TEST[c(1, 400, 800)]) rec_keep[[s]] <- list(obs = obs, rec = mfp, ti = ti, ps = ps)
}
rec <- rbindlist(rec_rows)
fwrite(rec, file.path(args$output, "reconstruction_by_subject.csv"))
rec_sum <- rec[, .(subjects = .N, rmse = mean(rmse), corr = mean(corr),
                   motion_corr = mean(motion_corr, na.rm = TRUE)), by = .(frames, method)]
fwrite(rec_sum, file.path(args$output, "reconstruction_summary.csv"))
print(rec_sum)

# ---------------------------------------------------------------- figures ----
to_img <- function(v) { z <- matrix(NA_real_, p, q); z[mvec] <- v; z }
show <- function(z, zlim = range(z, na.rm = TRUE), col = gray.colors(128, 0, 1), main = "") {
  image(t(z)[, p:1], col = col, zlim = zlim, axes = FALSE, main = main, asp = 1, useRaster = TRUE)
}
div_col <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(129)

# 1. mean dynamics
png(file.path(fig_dir, "mean_dynamics.png"), 1500, 330, res = 110)
par(mfrow = c(1, 6), mar = c(0.5, 0.5, 2, 0.5))
gi <- vapply(c(0, 0.25, 0.5, 0.75, 1), function(v) which.min(abs(t_grid - v)), 1L)
zl <- range(mu[gi, mvec])
for (g in gi) show(to_img(mu[g, mvec]), zl, main = sprintf("mu(t = %.2f)", t_grid[g]))
d <- mu[nt, mvec] - mu[1, mvec]
show(to_img(d), c(-1, 1) * max(abs(d)), div_col, "mu(1) - mu(0)")
dev.off()

# 2. temporal eigenfunctions
png(file.path(fig_dir, "temporal_modes.png"), 1300, 420, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
nshow <- min(4L, n_tmodes)
cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
matplot(t_grid, tfun[, 1:nshow], type = "l", lty = 1, lwd = 2, col = cols[1:nshow],
        xlab = "cardiac phase t (0 = ED, 1 = ES)", ylab = "theta_l(t)",
        main = "Temporal eigenfunctions")
abline(h = 0, col = "grey70")
legend("topleft", sprintf("l = %d (%.1f%%)", 1:nshow, 100 * temporal_diag$marginal_variance_fraction[1:nshow]),
       col = cols[1:nshow], lwd = 2, bty = "n", cex = 0.85)
barplot(100 * temporal_diag$marginal_variance_fraction, names.arg = temporal_diag$time_component,
        col = ifelse(temporal_diag$dynamic_screen, "#d1495b", "#c8ccd2"),
        ylab = "% of temporal variance", xlab = "temporal mode",
        main = sprintf("Dynamic modes carry %.1f%% (ADNI: 0.98%%)", 100 * dyn_share))
legend("topright", c("dynamic (non-constant)", "near-constant"), fill = c("#d1495b", "#c8ccd2"), bty = "n")
dev.off()

# 3. leading product FPCs: eigenimage phi_j psi_k^T and theta_l
nprod <- 8L
png(file.path(fig_dir, "leading_product_fpcs.png"), 1700, 470, res = 110)
par(mfrow = c(2, nprod), mar = c(0.5, 0.5, 2.2, 0.5))
for (f in 1:nprod) {
  e <- Phi[, prod_idx$j[f]] %o% Psi[, prod_idx$k[f]]
  show(e, c(-1, 1) * max(abs(e[mask])), div_col,
       sprintf("FPC%d (%.1f%%)\nrow %d, col %d, time %d", f, 100 * prod_idx$pve[f],
               prod_idx$j[f], prod_idx$k[f], prod_idx$l[f]))
}
par(mar = c(2.5, 2.5, 0.5, 0.5))
for (f in 1:nprod) {
  plot(t_grid, Theta[, prod_idx$l[f]], type = "l", lwd = 2, col = cols[min(prod_idx$l[f], 4)],
       xlab = "", ylab = "", ylim = range(Theta[, 1:min(4, ncol(Theta))]))
  abline(h = 0, col = "grey70")
}
dev.off()

# 4. leading "eigen-videos": PCs of the MatFPCA-predicted videos (sparse
# frames), each a combination of many product FPCs. Shown as the deviation
# from the mean video at five phases (+1 sd), with its correlation with EF.
png(file.path(fig_dir, "video_pcs.png"), 1500, 1050, res = 110)
Ftr <- scores_full$sparse$TRAIN
vpc <- svd(sweep(Ftr, 2, colMeans(Ftr)), nu = 0, nv = 4)
vz <- sweep(Ftr, 2, colMeans(Ftr)) %*% vpc$v
n_full <- prod(ranks)
par(mfrow = c(4, 6), mar = c(0.3, 0.3, 2.2, 0.3))
for (pc in 1:4) {
  full <- numeric(n_full); full[prod_col] <- vpc$v[, pc] * sd(vz[, pc])
  dev_v <- rTensor::ttl(rTensor::as.tensor(array(full, ranks)),
                        list(Phi, Psi, Theta[gi, , drop = FALSE]), ms = 1:3)@data
  lim <- c(-1, 1) * max(abs(dev_v[rep(mask, length(gi))]))
  rho <- cor(vz[, pc], y$TRAIN, method = "spearman")
  for (g in seq_along(gi)) {
    z <- dev_v[, , g]; z[!mask] <- NA
    show(z, lim, div_col, sprintf("video PC%d, t = %.2f", pc, t_grid[gi[g]]))
  }
  dev_all <- rTensor::ttl(rTensor::as.tensor(array(full, ranks)),
                          list(Phi, Psi, Theta), ms = 1:3)@data      # p x q x nt
  amp <- apply(dev_all, 3, function(z) mean(abs(z[mask])))
  share <- var(vz[, pc]) / sum(apply(Ftr, 2, var))
  plot(t_grid, amp, type = "l", lwd = 2, col = "#1f3f8f", xlab = "", ylab = "",
       axes = FALSE, ylim = c(0, max(amp)),
       main = sprintf("rho(EF) %+.2f, %.1f%% var", rho, 100 * share))
  mtext("mean |deviation| vs phase", side = 1, line = 0.2, cex = 0.7)
  axis(1, at = c(0, 0.5, 1), labels = c("ED", "", "ES")); box()
}
dev.off()

# 5. EF predicted vs true (TEST)
png(file.path(fig_dir, "ef_test_scatter.png"), 1800, 500, res = 110)
show_p <- intersect(c("ED frame only (static)", "Stacked sparse video (linear interpolation, 11 phases)",
                      "MatFPCA video PCA sparse", "Stacked dense video (11 phases)"), names(preds))
par(mfrow = c(1, length(show_p)), mar = c(4, 4, 2.5, 1))
for (nm in show_p) {
  r2 <- 1 - sum((y$TEST - preds[[nm]])^2) / sum((y$TEST - mean(y$TEST))^2)
  plot(y$TEST, preds[[nm]], pch = 16, cex = 0.5, col = adjustcolor("#1f3f8f", 0.4),
       xlab = "true EF (%)", ylab = "predicted EF (%)", xlim = c(5, 95), ylim = c(5, 95),
       main = sprintf("%s\nk = %d, R2 = %.2f", sub(" \\(linear.*", " (interp.)", nm), k_show, r2), cex.main = 0.9)
  abline(0, 1, col = "grey50"); abline(v = 40, h = 40, lty = 3, col = "grey60")
}
dev.off()

# 6. reconstruction strips
png(file.path(fig_dir, "reconstruction_examples.png"), 1500, 270 * 2 * length(rec_keep), res = 110)
par(mfrow = c(2 * length(rec_keep), 6), mar = c(0.3, 0.3, 1.8, 0.3))
for (s in names(rec_keep)) {
  r <- rec_keep[[s]]
  sel <- unique(round(seq(1, length(r$ti), length.out = 6)))
  zl2 <- quantile(r$obs, c(0.01, 0.99))
  r$obs <- pmin(pmax(r$obs, zl2[1]), zl2[2]); r$rec <- pmin(pmax(r$rec, zl2[1]), zl2[2])
  for (j in sel) show(to_img(r$obs[j, ]), zl2,
                      main = sprintf("%s t=%.2f%s", substr(s, 1, 8), t_grid[r$ti[j]],
                                     if (j %in% r$ps) " (seen)" else ""))
  for (j in sel) show(to_img(r$rec[j, ]), zl2, main = "MatFPCA prediction")
}
dev.off()

# ---------------------------------------------------------------- summary ----
best <- function(pattern, frames = NULL) {
  d <- ef_table[grepl(pattern, method, fixed = TRUE)]
  d[which.max(R2)]
}
summary_lines <- c(
  sprintf("design: %s | train subjects: %d | train frames: %d | image %dx%d | mask %d px",
          args$design, length(train_ids), sum(lengths(train_pos)), p, q, sum(mask)),
  sprintf("fit time: %.0f s | ranks I x J x K = %s | sigma2 = %.4g",
          fit$elapsed, paste(ranks, collapse = " x "), C$sig2),
  sprintf("temporal variance in dynamic modes: %.1f%%  (go if >= 10%%)", 100 * dyn_share),
  sprintf("EF TEST R2, MatFPCA sparse frames (best k): %.3f", best("MatFPCA scores (sparse")$R2),
  sprintf("EF TEST R2, MatFPCA dense frames (best k):  %.3f", best("MatFPCA scores (dense")$R2),
  sprintf("EF TEST R2, static ED-frame PCA (best k):   %.3f", best("ED frame only")$R2),
  sprintf("score implementation check (relative diff): %.2g", score_check))
writeLines(summary_lines, file.path(args$output, "run_summary.txt"))
cat(summary_lines, sep = "\n")
log_message("Done.")
