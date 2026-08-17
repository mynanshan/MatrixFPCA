source(file.path("R", "matfpca.R"))
source(file.path("simulation", "gendata-simple.R"))

set.seed(20260625)

generated <- gen_data_simple(
  n = 20,
  p = 6,
  q = 5,
  npc_row = 6,
  npc_col = 5,
  npc_time = 3,
  nt = 21,
  sampling_ratio = c(0.20, 0.35),
  sig = 0.5,
  eigvecSet = 2,
  scores_distr = "normal",
  score_dependence = "common_amplitude",
  dependence_param = 1
)

stopifnot(
  identical(dim(generated$scores), c(20L, 6L, 5L, 3L)),
  abs(mean(generated$subject_amplitude^2) - 1) < 0.75,
  !"maxIter" %in% names(formals(fit_matfcov)),
  "temporal_method" %in% names(formals(fit_matfcov)),
  "nfolds" %in% names(formals(fit_matfcov)),
  "crossfit_diagonal" %in% names(formals(fit_matfcov)),
  identical(eval(formals(fit_matfcov)$nfolds), 10L),
  identical(
    eval(formals(fit_matfcov)$crossfit_diagonal),
    c("refined", "identity")
  )
)

fit_identity <- fit_matfcov(
  generated$data,
  wScheme = "identity",
  noise = TRUE,
  nbasis = "theory"
)
fit_refined <- fit_matfcov(
  generated$data,
  wScheme = "gaussian_refined",
  noise = TRUE,
  nbasis = "theory",
  opwIter = 2,
  nfolds = 2
)
fit_refined_legacy <- fit_matfcov(
  generated$data,
  wScheme = "gaussian_refined",
  noise = TRUE,
  nbasis = "theory",
  opwIter = 2,
  nfolds = 1
)
fit_refined_identity_diagonal <- fit_matfcov(
  generated$data,
  wScheme = "gaussian_refined",
  noise = TRUE,
  nbasis = "theory",
  opwIter = 2,
  nfolds = 2,
  crossfit_diagonal = "identity"
)
fit_joint_identity <- fit_matfcov(
  generated$data,
  wScheme = "identity",
  noise = TRUE,
  nbasis = "theory",
  temporal_method = "joint"
)
fit_joint_refined <- fit_matfcov(
  generated$data,
  wScheme = "gaussian_refined",
  noise = TRUE,
  nbasis = "theory",
  opwIter = 2,
  temporal_method = "joint",
  nfolds = 2
)
fit_joint_refined_identity_diagonal <- fit_matfcov(
  generated$data,
  wScheme = "gaussian_refined",
  noise = TRUE,
  nbasis = "theory",
  opwIter = 2,
  temporal_method = "joint",
  nfolds = 2,
  crossfit_diagonal = "identity"
)

for (fit in list(
  fit_identity,
  fit_refined,
  fit_refined_legacy,
  fit_refined_identity_diagonal,
  fit_joint_identity,
  fit_joint_refined,
  fit_joint_refined_identity_diagonal
)) {
  covariance <- fit$C
  stopifnot(
    abs(tr(covariance$U) - 6) < 1e-8,
    abs(tr(covariance$V) - 5) < 1e-8,
    min(eigen(covariance$U, symmetric = TRUE, only.values = TRUE)$values) > -1e-8,
    min(eigen(covariance$V, symmetric = TRUE, only.values = TRUE)$values) > -1e-8,
    min(eigen(covariance$K, symmetric = TRUE, only.values = TRUE)$values) > -1e-8,
    covariance$sig2 >= 0,
    !is.null(fit$raw$Gamma_U),
    !is.null(fit$raw$Gamma_V),
    !is.null(fit$raw$U_unprojected),
    !is.null(fit$raw$V_unprojected),
    identical(names(fit$temporal_fit$fit), c("K", "G")),
    identical(names(fit$temporal_fit$nbasis), c("K", "G")),
    fit$temporal_fit$temporal_method %in% c("two_step", "joint")
  )
}

stopifnot(
  is.null(fit_identity$W),
  inherits(fit_refined$W, "mode_weights"),
  inherits(fit_refined$W, "gaussian_refined"),
  inherits(fit_refined$W, "crossfit_mode_weights"),
  length(fit_refined$pilot) == 2L,
  length(fit_refined$beta_K) == 20L,
  all(is.finite(fit_refined$beta_K)),
  all(is.finite(fit_refined$alpha_K)),
  all(fit_refined$alpha_K > 0),
  fit_refined$temporal_fit$crossfit_scale > 0,
  length(fit_refined$fold_id) == 20L,
  length(unique(fit_refined$fold_id)) == 2L,
  all(vapply(get_W(fit_refined$W, 1L, "K"), is.matrix, logical(1))),
  all(vapply(get_W(fit_refined$W, 2L, "K"), is.matrix, logical(1))),
  all(vapply(get_W(fit_refined$W, 3L, "U"), is.matrix, logical(1))),
  all(vapply(get_W(fit_refined$W, 3L, "V"), is.matrix, logical(1))),
  !inherits(fit_refined_legacy$W, "crossfit_mode_weights"),
  is(fit_refined_legacy$pilot, "matfcov"),
  identical(fit_refined$temporal_fit$crossfit_diagonal, "refined"),
  identical(
    fit_refined_identity_diagonal$temporal_fit$crossfit_diagonal,
    "identity"
  ),
  identical(
    fit_joint_refined$temporal_fit$crossfit_diagonal,
    "refined"
  ),
  isTRUE(all.equal(
    fit_refined$C$sig2,
    fit_refined$temporal_fit$sig2
  )),
  is.null(fit_identity$temporal_fit$fit$G),
  is.null(fit_refined$temporal_fit$fit$G),
  "nugget_multiplier" %in%
    names(stats::coef(fit_joint_refined$temporal_fit$fit$K))
)

# The two-step estimator uses the manuscript's direct weighted average of
# diagonal residuals and does not fit a separate diagonal smoother.
identity_diagonal <- do.call(rbind, lapply(seq_len(count_subj(generated$data)), function(i) {
  Xi <- get_subj_tnsr(generated$data, i)
  tid <- generated$data$t_idx_list[[i]]
  data.frame(
    t = generated$data$t_grid[tid],
    Kval = apply(Xi@data^2, 3L, sum) / prod(generated$data$vdim)
  )
}))
identity_residual <- identity_diagonal$Kval - as.numeric(stats::predict(
  fit_identity$temporal_fit$fit$K,
  newdata = data.table::data.table(
    t1 = identity_diagonal$t,
    t2 = identity_diagonal$t
  )
))
stopifnot(isTRUE(all.equal(
  fit_identity$C$sig2,
  max(0, mean(identity_residual)),
  tolerance = 1e-10
)))

refined_diagonal <- do.call(rbind, lapply(seq_len(count_subj(generated$data)), function(i) {
  Xi <- get_subj_tnsr(generated$data, i)
  tid <- generated$data$t_idx_list[[i]]
  data.frame(
    t = generated$data$t_grid[tid],
    Kval = diag(get_working_K(
      Xi,
      get_W(fit_refined$W, 1L, "K")[[i]],
      get_W(fit_refined$W, 2L, "K")[[i]]
    )) / (prod(generated$data$vdim) * fit_refined$beta_K[[i]]),
    rho = fit_refined$alpha_K[[i]] / fit_refined$beta_K[[i]]
  )
}))
refined_residual <- refined_diagonal$Kval - as.numeric(stats::predict(
  fit_refined$temporal_fit$fit$K,
  newdata = data.table::data.table(
    t1 = refined_diagonal$t,
    t2 = refined_diagonal$t
  )
))
refined_sig2 <- max(
  0,
  sum(refined_diagonal$rho * refined_residual) /
    sum(refined_diagonal$rho^2)
)
stopifnot(isTRUE(all.equal(
  fit_refined$C$sig2,
  refined_sig2,
  tolerance = 1e-10
)))

for (fold in sort(unique(fit_refined$fold_id))) {
  held_out <- which(fit_refined$fold_id == fold)
  training <- fit_refined$training_subjects[[paste0("fold", fold)]]
  stopifnot(!length(intersect(held_out, training)))
}

# A held-out response must not affect any weight subsequently applied to it.
modified_data <- generated$data
modified_X1 <- get_subj_tnsr(modified_data, 1L)
modified_data$X_list[[1L]] <- rTensor::as.tensor(modified_X1@data * 100)
modified_crossfit <- .compute_crossfit_weights(
  modified_data,
  fold_id = fit_refined$fold_id,
  wScheme = "gaussian_refined",
  noise = TRUE,
  nbasis = "theory",
  opwIter = 2,
  temporal_method = "two_step"
)
for (target_mode in list(
  c("K", "W1"), c("K", "W2"),
  c("U", "W2"), c("U", "W3"),
  c("V", "W1"), c("V", "W3")
)) {
  target <- target_mode[[1L]]
  mode <- target_mode[[2L]]
  stopifnot(isTRUE(all.equal(
    fit_refined$W[[target]][[mode]][[1L]],
    modified_crossfit$W[[target]][[mode]][[1L]]
  )))
}

stopifnot(
  identical(.make_weight_folds(6L, 3L), c(1L, 2L, 3L, 1L, 2L, 3L)),
  identical(
    .make_weight_folds(6L, 2L, c("a", "a", "b", "b", "a", "b")),
    c(1L, 1L, 2L, 2L, 1L, 2L)
  ),
  inherits(
    try(.make_weight_folds(4L, 5L), silent = TRUE),
    "try-error"
  )
)

# At exactly zero noise, the Gaussian-refined formula must reduce to the
# plug-in inverse-covariance scheme.
U_test <- diag(c(2, 1))
V_test <- diag(c(3, 1.5, 0.75))
LK_test <- list(diag(c(2, 1)), diag(c(1.5, 0.5, 0.25)))
weights_refined_zero <- compute_W(
  U_test, V_test, LK_test,
  sig2 = 0,
  wScheme = "gaussian_refined"
)
weights_inverse <- compute_W(
  U_test, V_test, LK_test,
  sig2 = 0,
  wScheme = "invcov"
)
stopifnot(
  inherits(weights_refined_zero, "shared"),
  isTRUE(all.equal(weights_refined_zero, weights_inverse))
)

singular_weight <- .gaussian_block_weight(
  diag(c(3, 1, 0)),
  A = 2,
  B = 0.1,
  C = 0.01,
  ridge = 0
)
stopifnot(abs(singular_weight[3, 3]) < 1e-12)

cat("Temporal estimator checks passed.\n")
