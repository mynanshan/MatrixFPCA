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
#   6. EF / ESV / EDV prediction: ridge on MatFPCA representations and on PCA
#      baselines; penalty and number of components k chosen on VAL, TEST
#      reported;
#   7. reconstruction of every TEST frame from the sparse frames;
#   8. interpretation: spatial and temporal modes, principal "eigen-videos"
#      (all and dynamic part only) with their correlates, EF-group means;
#   9. figures.

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
  n_train = 0L,        # use TRAIN subjects with train_rank in (train_skip, train_skip + n_train]
  train_skip = 0L,     # (0 = all exported TRAIN subjects)
  n_fit = 0L,          # fit mean + covariance on TRAIN subjects with train_rank in
  fit_skip = 0L,       # (fit_skip, fit_skip + n_fit]; 0 = all TRAIN. The EF regression
                       # and baselines always use all selected TRAIN subjects.
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
if (args$n_train > 0L) {
  subjects <- subjects[split != "TRAIN" |
                         (train_rank > args$train_skip & train_rank <= args$train_skip + args$n_train)]
  stopifnot(sum(subjects$split == "TRAIN") == args$n_train)
}
frames <- frames[id %in% subjects$id]
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

# sparse subset of frames per subject (TRAIN, VAL and TEST all get one). The
# draw is seeded per subject, so a subject keeps the same frames whichever
# other subjects are in the run (needed to compare fits on different samples).
sparse_pos <- Map(function(ti, sid) {
  set.seed(args$seed + strtoi(substr(sid, 3L, 9L), 16L))
  n <- length(ti)
  m <- min(n, sample(args$m_min:args$m_max, 1L))
  sort(sample.int(n, m))
}, t_idx_all, names(t_idx_all))

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
train_ids <- if (args$n_fit > 0L) {
  subjects[split == "TRAIN" & train_rank > args$fit_skip &
             train_rank <= args$fit_skip + args$n_fit, id]
} else subjects[split == "TRAIN", id]
if (args$n_fit > 0L) stopifnot(length(train_ids) == args$n_fit)
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
# Every representation is reduced to k features (PCA on TRAIN, or the top-k
# product-FPC scores), then ridge regression; the penalty and k are both
# chosen on VAL (RMSE) and TEST is reported. Targets: EF, ESV and EDV.
auc <- function(y, s) {
  r <- rank(c(s[y], s[!y]))
  (sum(r[seq_len(sum(y))]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(!y))
}
ridge_eval <- function(Ztr, ytr, Zva, yva, Zte, yte) {
  m <- colMeans(Ztr); s <- apply(Ztr, 2, sd); s[s == 0] <- 1
  f <- function(Z) cbind(1, sweep(sweep(Z, 2, m), 2, s, "/"))
  A <- f(Ztr); AtA <- crossprod(A); Aty <- crossprod(A, ytr); Ava <- f(Zva); best <- NULL
  for (lam in 10^seq(-2, 4, by = 0.5)) {
    b <- solve(AtA + diag(c(0, rep(lam, ncol(Ztr)))), Aty)
    rm <- sqrt(mean((Ava %*% b - yva)^2))
    if (is.null(best) || rm < best$rm) best <- list(rm = rm, lam = lam, b = b)
  }
  pr <- as.vector(f(Zte) %*% best$b)
  list(R2 = 1 - sum((yte - pr)^2) / sum((yte - mean(yte))^2), MAE = mean(abs(yte - pr)),
       lambda = best$lam, val_rmse = best$rm, pred = pr)
}
pca_all <- function(Flist, kmax) {     # PCA scores (TRAIN-fitted) for every split
  m <- colMeans(Flist$TRAIN)
  kmax <- min(kmax, nrow(Flist$TRAIN) - 1L, ncol(Flist$TRAIN))
  sv <- svd(sweep(Flist$TRAIN, 2, m), nu = 0, nv = kmax)
  lapply(Flist, function(F) sweep(F, 2, m) %*% sv$v)
}
targets <- c("EF", "ESV", "EDV")
yy <- setNames(lapply(targets, function(tg) lapply(ids_by, function(ids) subjects[match(ids, id)][[tg]])),
               targets)
y <- yy$EF

# image features for baselines (in-mask pixels, gain-normalised, not centred)
mvec <- as.vector(mask)
frame_mat <- function(x) t(matrix(x, p * q))[, mvec, drop = FALSE]   # m x P
interp_weights <- function(ph, xout) {   # linear interpolation, constant beyond the ends
  m <- length(ph); W <- matrix(0, length(xout), m)
  if (m == 1L) { W[, 1L] <- 1; return(W) }
  x <- pmin(pmax(xout, ph[1L]), ph[m])
  j <- pmin(findInterval(x, ph), m - 1L)
  w <- (x - ph[j]) / (ph[j + 1L] - ph[j])
  W[cbind(seq_along(x), j)] <- 1 - w
  W[cbind(seq_along(x), j + 1L)] <- W[cbind(seq_along(x), j + 1L)] + w
  W
}
interp11 <- function(x, ti) interp_weights(t_grid[ti], seq(0, 1, 0.1)) %*% frame_mat(x)
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
baseline_frames <- c(ed = "ED/ES landmarks", diff = "ED/ES landmarks", avg_sparse = "sparse",
                     avg_dense = "dense", stacked = "dense", stacked_sparse = "sparse")

ef_rows <- list(); preds <- list()
evaluate <- function(method, frames_lab, Z) {
  for (tg in targets) for (k in unique(pmin(n_fpcs, ncol(Z$TRAIN)))) {
    yt <- yy[[tg]]
    r <- ridge_eval(Z$TRAIN[, 1:k, drop = FALSE], yt$TRAIN, Z$VAL[, 1:k, drop = FALSE], yt$VAL,
                    Z$TEST[, 1:k, drop = FALSE], yt$TEST)
    ef_rows[[length(ef_rows) + 1L]] <<- data.table(
      target = tg, method = method, frames = frames_lab, k = k, R2 = r$R2, MAE = r$MAE,
      AUC_EF40 = if (tg == "EF") auc(yt$TEST < 40, -r$pred) else NA_real_,
      lambda = r$lambda, val_rmse = r$val_rmse)
    if (tg == "EF") preds[[method]][[as.character(k)]] <<- r$pred
  }
}
log_message("Prediction of EF / ESV / EDV: MatFPCA representations and PCA baselines...")
static_cols <- prod_idx$l == 1L
for (w in c("sparse", "dense")) {
  evaluate(sprintf("MatFPCA product-FPC scores, top k (%s frames)", w), w, scores_full[[w]])
  evaluate(sprintf("MatFPCA predicted video -> PCA (%s frames)", w), w,
           pca_all(scores_full[[w]], max(n_fpcs)))
}
# the separable model splits each predicted video into a part carried by the
# (near-constant) first temporal mode and a part carried by the dynamic modes
evaluate("MatFPCA static part only (time mode 1) -> PCA (sparse frames)", "sparse",
         pca_all(lapply(scores_full$sparse, function(S) S[, static_cols, drop = FALSE]), max(n_fpcs)))
evaluate("MatFPCA dynamic part only (time modes 2+) -> PCA (sparse frames)", "sparse",
         pca_all(lapply(scores_full$sparse, function(S) S[, !static_cols, drop = FALSE]), max(n_fpcs)))
for (nm in names(baseline_kinds)) {
  F <- lapply(ids_by, feat, kind = baseline_kinds[[nm]])
  evaluate(paste("PCA:", nm), baseline_frames[[baseline_kinds[[nm]]]], pca_all(F, max(n_fpcs)))
  rm(F); gc(verbose = FALSE)
}
ef_all <- rbindlist(ef_rows)
fwrite(ef_all, file.path(args$output, "prediction_all_k.csv"))
ef_table <- ef_all[target == "EF"]
fwrite(ef_table, file.path(args$output, "ef_prediction.csv"))
# k chosen on VAL
ef_sel <- ef_all[ef_all[, .I[which.min(val_rmse)], by = .(target, method)]$V1]
fwrite(ef_sel, file.path(args$output, "prediction_selected_k.csv"))
print(ef_sel[target == "EF", .(method, k, R2 = round(R2, 3), MAE = round(MAE, 2), AUC = round(AUC_EF40, 3))])
sel_pred <- function(m) preds[[m]][[as.character(ef_sel[target == "EF" & method == m, k])]]

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
  z <- pmin(pmax(z, zlim[1]), zlim[2])            # clamp, so out-of-range pixels are not blank
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

# 5. EF predicted vs true (TEST), k chosen on VAL
png(file.path(fig_dir, "ef_test_scatter.png"), 1800, 500, res = 110)
show_p <- c("PCA: ED frame only (static)",
            "PCA: Stacked sparse video (linear interpolation, 11 phases)",
            "MatFPCA predicted video -> PCA (sparse frames)",
            "PCA: Stacked dense video (11 phases)")
show_lab <- c("ED frame (static)", "sparse frames, linear interpolation",
              "MatFPCA video, sparse frames", "all frames (dense video)")
par(mfrow = c(1, length(show_p)), mar = c(4, 4, 2.5, 1))
for (i in seq_along(show_p)) {
  pr <- sel_pred(show_p[i]); r <- ef_sel[target == "EF" & method == show_p[i]]
  plot(y$TEST, pr, pch = 16, cex = 0.5, col = adjustcolor("#1f3f8f", 0.4),
       xlab = "true EF (%)", ylab = "predicted EF (%)", xlim = c(5, 95), ylim = c(5, 95),
       main = sprintf("%s\nk = %d (chosen on VAL), R2 = %.2f", show_lab[i], r$k, r$R2),
       cex.main = 0.9)
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

# --------------------------------------------------------- interpretation ----
# What does each mode show?
#  (a) spatial modes: leading row (top-bottom) and column (left-right)
#      eigenvectors next to the mean ED image;
#  (b) principal eigen-videos: PCA of the MatFPCA-predicted videos (all
#      parts, and the dynamic part only), each shown as the mean video
#      -/+ 2 SD at ED and ES, with Spearman correlations against clinical
#      and acquisition variables (all splits, sparse frames);
#  (c) EF groups: mean ED / ES images and systolic change for reduced
#      (EF < 40%) and normal (EF >= 55%) TEST hearts, observed from all
#      frames vs MatFPCA-predicted from the sparse frames.
log_message("Interpretation outputs...")
gray_col <- gray.colors(128, 0, 1)
all_ids <- unlist(ids_by, use.names = FALSE)
covars <- c(EF = "EF", EDV = "EDV", ESV = "ESV", gain = "norm_scale", LV_length = "axis_len",
            LV_tilt = "axis_tilt", fps = "FPS", systole_s = "gap_s")
cov_dt <- subjects[match(all_ids, id), ..covars]
setnames(cov_dt, names(covars))
score_video <- function(v, tt = seq_len(nt)) {     # product-score vector -> p x q x |tt|
  full <- numeric(prod(ranks)); full[prod_col] <- v
  rTensor::ttl(rTensor::as.tensor(array(full, ranks)),
               list(Phi, Psi, Theta[tt, , drop = FALSE]), ms = 1:3)@data
}
ed_es <- c(1L, nt)

# (a) spatial modes
png(file.path(fig_dir, "spatial_modes.png"), 1300, 1000, res = 110)
layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE), widths = c(1, 0.8), heights = c(1, 0.8))
par(mar = c(1, 1, 2.5, 1))
show(to_img(mu[1, mvec]), main = "mean video at ED (t = 0)")
pve_phi <- eig$values$lam_phi / sum(eig$values$lam_phi)
pve_psi <- eig$values$lam_psi / sum(eig$values$lam_psi)
par(mar = c(4, 4, 2.5, 1))
matplot(Phi[, 1:4], seq_len(p), type = "l", lty = 1, lwd = 2, col = cols, ylim = c(p, 1),
        xlab = "loading", ylab = "image row (top = apex, bottom = atria)",
        main = "row eigenvectors phi_j (vertical profiles)")
abline(v = 0, col = "grey70")
legend("bottomright", sprintf("j = %d (%.1f%%)", 1:4, 100 * pve_phi[1:4]), col = cols, lwd = 2,
       bty = "n", cex = 0.8)
matplot(seq_len(q), Psi[, 1:4], type = "l", lty = 1, lwd = 2, col = cols,
        xlab = "image column (left = septum / right ventricle, right = lateral wall)",
        ylab = "loading", main = "column eigenvectors psi_k (horizontal profiles)")
abline(h = 0, col = "grey70")
legend("topright", sprintf("k = %d (%.1f%%)", 1:4, 100 * pve_psi[1:4]), col = cols, lwd = 2,
       bty = "n", cex = 0.8)
par(mar = c(4, 4, 2.5, 1))
plot(cumsum(pve_phi)[1:min(30, length(pve_phi))], type = "b", pch = 16, cex = 0.6, col = "#1f3f8f",
     ylim = c(0, 1), xlab = "number of spatial components", ylab = "cumulative share",
     main = "spatial eigenvalues are spread out")
lines(cumsum(pve_psi)[1:min(30, length(pve_psi))], type = "b", pch = 16, cex = 0.6, col = "#d95f02")
legend("bottomright", c("rows (U)", "columns (V)"), col = c("#1f3f8f", "#d95f02"), lwd = 2, bty = "n")
dev.off()

# (b) principal eigen-videos
eigen_videos <- function(cols_use, npc) {
  Ftr <- scores_full$sparse$TRAIN[, cols_use, drop = FALSE]
  m <- colMeans(Ftr)
  sv <- svd(sweep(Ftr, 2, m), nu = 0, nv = npc)
  tot <- sum(apply(scores_full$sparse$TRAIN, 2, var))       # total predicted-video variance
  Zall <- sweep(do.call(rbind, scores_full$sparse)[, cols_use, drop = FALSE], 2, m) %*% sv$v
  lapply(seq_len(npc), function(j) {
    v <- numeric(ncol(scores_full$sparse$TRAIN)); v[which(cols_use)] <- sv$v[, j]
    sdj <- sd(Zall[seq_along(ids_by$TRAIN), j])
    # orient so that +2 SD corresponds to higher EF (purely cosmetic)
    sgn <- if (cor(Zall[, j], cov_dt$EF) < 0) -1 else 1
    vid <- score_video(sgn * v * sdj)
    dyn <- sum(sweep(vid, 1:2, apply(vid, 1:2, mean))^2) / sum(vid^2)
    list(z = sgn * Zall[, j], vid = vid, share = sdj^2 / tot, dyn = dyn)
  })
}
npc <- 6L
ev <- list(all = eigen_videos(rep(TRUE, length(prod_col)), npc),
           dynamic = eigen_videos(!static_cols, npc))
interp_tab <- rbindlist(lapply(names(ev), function(part) rbindlist(lapply(seq_len(npc), function(j) {
  e <- ev[[part]][[j]]
  rho <- vapply(cov_dt, function(v) cor(e$z, v, method = "spearman", use = "complete.obs"), 1)
  data.table(part = part, pc = j, variance_share = e$share, time_varying_share = e$dyn,
             as.data.table(as.list(setNames(rho, paste0("rho_", names(rho))))))
}))))
fwrite(interp_tab, file.path(args$output, "eigen_video_interpretation.csv"))
print(interp_tab[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)])

top_corr <- function(r) {
  rr <- unlist(r[, grep("^rho_", names(r)), with = FALSE])
  o <- order(-abs(rr))[1:3]
  paste(sprintf("%s %+.2f", sub("rho_", "", names(rr)[o]), rr[o]), collapse = ", ")
}
for (pt in names(ev)) {
  png(file.path(fig_dir, sprintf("eigen_videos_%s.png", pt)), 1700, 300 * npc, res = 110)
  par(mfrow = c(npc, 7), mar = c(0.3, 0.3, 2.4, 0.3), oma = c(0, 0, 2, 0))
  for (j in seq_len(npc)) {
    e <- ev[[pt]][[j]]
    lim_g <- range(mu[ed_es, mvec])
    for (g in ed_es) for (sg in c(-2, 2)) {
      z <- mu_arr[, , g] + sg * e$vid[, , g]; z[!mask] <- NA
      show(z, lim_g, gray_col, sprintf("PC%d %s2 SD, %s", j, if (sg < 0) "-" else "+",
                                       if (g == 1L) "ED" else "ES"))
    }
    dl <- c(-1, 1) * max(abs(e$vid[rep(mask, nt)]))
    for (g in ed_es) { z <- e$vid[, , g]; z[!mask] <- NA
      show(z, dl, div_col, sprintf("PC%d deviation at %s", j, if (g == 1L) "ED" else "ES")) }
    plot.new()
    r <- interp_tab[part == pt & pc == j]
    text(0, 0.8, sprintf("PC%d: %.1f%% of variance", j, 100 * r$variance_share), adj = 0, cex = 0.95)
    text(0, 0.6, sprintf("time-varying: %.0f%%", 100 * r$time_varying_share), adj = 0, cex = 0.95)
    text(0, 0.4, "strongest correlates:", adj = 0, cex = 0.85)
    text(0, 0.25, top_corr(r), adj = 0, cex = 0.8)
  }
  mtext(sprintf("Principal eigen-videos of the MatFPCA-predicted videos (%s); +2 SD = higher-EF side",
                if (pt == "all") "all parts" else "dynamic part only: time modes 2+"),
        outer = TRUE, cex = 1)
  dev.off()
}

# (c) EF groups (TEST): observed from all frames vs MatFPCA from sparse frames
grp <- cut(y$TEST, c(-Inf, 40, 55, Inf), right = FALSE, labels = c("reduced", "borderline", "normal"))
obs_edes <- lapply(c(ED = "ed", ES = "es"), function(w) t(vapply(ids_by$TEST, function(s) {
  x <- raw[[s]]; frame_mat(x[, , if (w == "ed") 1L else dim(x)[3], drop = FALSE])[1, ] }, numeric(sum(mvec)))))
gm <- list(obs = list(), pred = list())
for (g in c("reduced", "normal")) {
  gm$obs[[g]] <- lapply(obs_edes, function(M) colMeans(M[grp == g, , drop = FALSE]))
  vid <- score_video(colMeans(scores_full$sparse$TEST[grp == g, , drop = FALSE]), ed_es)
  gm$pred[[g]] <- list(ED = vid[, , 1][mvec] + mu[1, mvec], ES = vid[, , 2][mvec] + mu[nt, mvec])
}
chg <- function(src, g) gm[[src]][[g]]$ES - gm[[src]][[g]]$ED
png(file.path(fig_dir, "ef_groups.png"), 1700, 560, res = 110)
par(mfrow = c(2, 7), mar = c(0.3, 0.3, 2.4, 0.3), oma = c(0, 3, 2, 0))
lim_g <- range(unlist(gm$obs))
cl <- c(-1, 1) * max(abs(c(chg("obs", "reduced"), chg("obs", "normal"))))
cd <- c(-1, 1) * max(abs(chg("obs", "reduced") - chg("obs", "normal")))
agree <- cor(chg("obs", "reduced") - chg("obs", "normal"),
             chg("pred", "reduced") - chg("pred", "normal"))
for (src in c("obs", "pred")) {
  for (g in c("reduced", "normal")) for (ph in c("ED", "ES"))
    show(to_img(gm[[src]][[g]][[ph]]), lim_g, gray_col,
         sprintf("%s EF, %s (n=%d)", g, ph, sum(grp == g)))
  show(to_img(chg(src, "reduced")), cl, div_col, "ES - ED, reduced EF")
  show(to_img(chg(src, "normal")), cl, div_col, "ES - ED, normal EF")
  show(to_img(chg(src, "reduced") - chg(src, "normal")), cd, div_col, "change: reduced - normal")
  mtext(if (src == "obs") "observed\n(all frames)" else "MatFPCA\n(3-8 frames)",
        side = 2, outer = TRUE, at = if (src == "obs") 0.75 else 0.25, line = 0.5, cex = 0.8)
}
mtext(sprintf("TEST hearts by EF group: reduced EF < 40%%, normal EF >= 55%%; correlation of the two 'reduced - normal' change maps = %.2f",
              agree), outer = TRUE, cex = 0.95)
dev.off()
writeLines(sprintf("%.4f", agree), file.path(args$output, "ef_group_change_agreement.txt"))

# ---------------------------------------------------------------- summary ----
sel <- function(m) ef_sel[target == "EF" & method == m]
summary_lines <- c(
  sprintf("design: %s | train subjects: %d | train frames: %d | image %dx%d | mask %d px",
          args$design, length(train_ids), sum(lengths(train_pos)), p, q, sum(mask)),
  sprintf("fit time: %.0f s | ranks I x J x K = %s | sigma2 = %.4g",
          fit$elapsed, paste(ranks, collapse = " x "), C$sig2),
  sprintf("temporal variance in dynamic modes: %.1f%%  (go if >= 10%%)", 100 * dyn_share),
  sprintf("EF TEST R2 (k on VAL), MatFPCA video PCA, sparse frames: %.3f (k = %d)",
          sel("MatFPCA predicted video -> PCA (sparse frames)")$R2, sel("MatFPCA predicted video -> PCA (sparse frames)")$k),
  sprintf("EF TEST R2 (k on VAL), MatFPCA video PCA, dense frames:  %.3f (k = %d)",
          sel("MatFPCA predicted video -> PCA (dense frames)")$R2, sel("MatFPCA predicted video -> PCA (dense frames)")$k),
  sprintf("EF TEST R2 (k on VAL), sparse interpolation PCA:         %.3f",
          sel("PCA: Stacked sparse video (linear interpolation, 11 phases)")$R2),
  sprintf("EF TEST R2 (k on VAL), static ED-frame PCA:              %.3f",
          sel("PCA: ED frame only (static)")$R2),
  sprintf("held-out frame RMSE, MatFPCA / mean + static PCA: %.3f / %.3f",
          rec_sum[frames == "held-out" & method == "MatFPCA", rmse],
          rec_sum[frames == "held-out" & method == "Mean + static PCA (20)", rmse]),
  sprintf("score implementation check (relative diff): %.2g", score_check))
writeLines(summary_lines, file.path(args$output, "run_summary.txt"))
cat(summary_lines, sep = "\n")
log_message("Done.")
