#!/usr/bin/env Rscript

# Benchmark MatFPCA against nonseparable methods for sparse longitudinal images.
# Run from the repository root; see README-method-comparison.md for rationale.

parser <- argparse::ArgumentParser(
  description = "Compare MatFPCA with sparse nonseparable image/functional PCA methods."
)
parser$add_argument("-n", "--n_subjects", type = "integer", default = 50L)
parser$add_argument("-p", "--d_row", type = "integer", default = 4L)
parser$add_argument("-q", "--d_col", type = "integer", default = 4L)
parser$add_argument("-nt", "--n_time_points", type = "integer", default = 51L)
parser$add_argument("-s1", "--s_obs_min", type = "integer", default = 2L,
                    help = "Minimum sampling ratio in units of 0.05.")
parser$add_argument("-s2", "--s_obs_max", type = "integer", default = 4L,
                    help = "Maximum sampling ratio in units of 0.05.")
parser$add_argument("-sig", "--noise_lvl", type = "double", default = 0.1)
parser$add_argument("-rr", "--rank_row", type = "integer", default = 0L)
parser$add_argument("-rc", "--rank_col", type = "integer", default = 0L)
parser$add_argument("-rt", "--rank_time", type = "integer", default = 5L)
parser$add_argument("-set", "--eigset", type = "integer", default = 2L)
parser$add_argument("-eval", "--eigvalset", type = "integer", default = 1L)
parser$add_argument("--eigensystem", choices = c("ar1", "spectral"), default = "ar1")
parser$add_argument("-scd", "--scores_distr",
                    choices = c("normal", "uniform", "beta32", "t5"),
                    default = "normal")
parser$add_argument("-dep", "--score_dependence",
                    choices = c("independent", "common_amplitude", "common_signed", "heavy_amplitude"),
                    default = "independent")
parser$add_argument("-dpar", "--dependence_param", type = "double", default = 1)
parser$add_argument("-nb", "--nbasis", type = "integer", default = 0L,
                    help = "MatFPCA basis count; 0 uses the theory-guided choice.")
parser$add_argument("--spline_nbasis", type = "integer", default = 7L)
parser$add_argument("--spline_smoothness", type = "double", default = 1e-4)
parser$add_argument("--happ_npc", type = "integer", default = 3L,
                    help = "PACE components retained per pixel for HAPP-MFPCA.")
parser$add_argument("--face_knots", type = "integer", default = 5L)
parser$add_argument("--face_pve", type = "double", default = 0.999)
parser$add_argument("--mfaces_knots", type = "integer", default = 5L)
parser$add_argument("--mfaces_pve", type = "double", default = 0.999)
parser$add_argument("--mfaces_search_length", type = "integer", default = 10L)
parser$add_argument("--mfaces_max_outcomes", type = "integer", default = 25L)
parser$add_argument("--mgs_rank_candidates", default = "3",
                    help = "Comma-separated per-pixel ranks; multiple values invoke AIC selection.")
parser$add_argument("--mgs_basis_candidates", default = "7",
                    help = "Comma-separated per-pixel spline basis sizes.")
parser$add_argument("--mgs_mu_nbasis", type = "integer", default = 7L)
parser$add_argument("--mgs_maxit", type = "integer", default = 300L)
parser$add_argument("--product_spatial_npc", type = "integer", default = 10L)
parser$add_argument("--product_temporal_npc", type = "integer", default = 5L)
parser$add_argument("--shi_ncomp", type = "integer", default = 5L)
parser$add_argument("--shi_nbasis_time", type = "integer", default = 7L)
parser$add_argument("--shi_nbasis_row", type = "integer", default = 5L)
parser$add_argument("--shi_nbasis_col", type = "integer", default = 5L)
parser$add_argument("--shi_maxit", type = "integer", default = 100L)
parser$add_argument("--shi_tolerance", type = "double", default = 1e-5)
parser$add_argument("--shi_smoothness", type = "double", default = 1e-2)
parser$add_argument("--shi_score_ridge", default = "auto",
                    help = "Shi score ridge multiplier, or 'auto' for visit-block CV.")
parser$add_argument("--spaco_rank", type = "integer", default = 5L)
parser$add_argument("--spaco_nbasis", type = "integer", default = 7L)
parser$add_argument("--spaco_smoothness", type = "double", default = 1e-2)
parser$add_argument("--spaco_ridge", type = "double", default = 1e-8,
                    help = "SPACO numerical diagonal-loading multiplier.")
parser$add_argument("--spaco_score_ridge", default = "auto",
                    help = "SPACO score ridge multiplier, or 'auto' for visit-block CV.")
parser$add_argument("--spaco_maxit", type = "integer", default = 100L)
parser$add_argument("--allow_slow_competitors", type = "integer", default = 0L,
                    help = "Allow mFACEs beyond --mfaces_max_outcomes.")
parser$add_argument("--eval_ncomp", type = "integer", default = 5L,
                    help = "Leading global eigenspace dimension to evaluate.")
parser$add_argument("--cov_ncomp", type = "integer", default = 0L,
                    help = paste(
                      "Common covariance rank for every method, including MatFPCA;",
                      "0 retains each method's full numerical rank."
                    ))
parser$add_argument("--methods", type = "character",
                    default = paste(c(
                      "matfpca_refined", "matfpca_identity", "happ_pace",
                      "happ_face", "mfaces", "mgsfpca", "shi_mdfpca",
                      "spaco_minus", "sparse_product_fpca", "spline_pca"
                    ), collapse = ","))
parser$add_argument("--weight_folds", type = "integer", default = 5L)
parser$add_argument("--temporal_method", choices = c("two_step", "joint"), default = "two_step")
parser$add_argument("--crossfit_diagonal", choices = c("refined", "identity"), default = "refined")
parser$add_argument("--ncores", type = "integer", default = 1L)
parser$add_argument("--blas_threads", type = "integer", default = 0L,
                    help = "BLAS/OpenMP threads; 0 uses 1 when method-level parallelism is active.")
parser$add_argument("--backend", choices = c("multisession", "multicore"), default = "multisession")
parser$add_argument("--seed", type = "integer", default = 1L)
parser$add_argument("--save", type = "integer", default = ifelse(interactive(), 0L, 1L))
parser$add_argument("--save_fits", type = "integer", default = 0L)
parser$add_argument("--results_dir", default = file.path("simulation", "results-methods"))
parser$add_argument("--result_tag", default = "",
                    help = "Optional scheduler/job tag appended to the result filename.")
args <- parser$parse_args()

required_packages <- c("fda", "rTensor", "data.table", "Matrix")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1L), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

source(file.path("R", "matfdata.R"))
source(file.path("R", "matfpca.R"))
source(file.path("simulation", "gendata-simple.R"))
source(file.path("simulation", "competitors", "core.R"))

n <- max(2L, args$n_subjects)
p <- max(1L, args$d_row)
q <- max(1L, args$d_col)
nt <- max(11L, args$n_time_points)
sig <- max(0, args$noise_lvl)
s1 <- max(1L, args$s_obs_min)
s2 <- min(max(s1 + 1L, args$s_obs_max), 20L)
sampling_ratio <- 0.05 * c(s1, s2)
npc_row <- if (args$rank_row <= 0L) p else min(p, max(1L, args$rank_row))
npc_col <- if (args$rank_col <= 0L) q else min(q, max(1L, args$rank_col))
npc_time <- min(nt, max(1L, args$rank_time))
nbasis <- if (args$nbasis <= 0L) "theory" else max(3L, min(nt, args$nbasis))
eval_ncomp <- max(1L, args$eval_ncomp)
cov_ncomp <- if (args$cov_ncomp <= 0L) NULL else max(eval_ncomp, args$cov_ncomp)
if (!is.finite(args$spaco_ridge) || args$spaco_ridge < 0) {
  stop("spaco_ridge must be finite and nonnegative.")
}
spaco_score_ridge <- if (tolower(args$spaco_score_ridge) == "auto") {
  "auto"
} else {
  parsed <- suppressWarnings(as.numeric(args$spaco_score_ridge))
  if (length(parsed) != 1L || !is.finite(parsed) || parsed < 0) {
    stop("spaco_score_ridge must be 'auto' or a finite nonnegative number.")
  }
  parsed
}
shi_score_ridge <- if (tolower(args$shi_score_ridge) == "auto") {
  "auto"
} else {
  parsed <- suppressWarnings(as.numeric(args$shi_score_ridge))
  if (length(parsed) != 1L || !is.finite(parsed) || parsed < 0) {
    stop("shi_score_ridge must be 'auto' or a finite nonnegative number.")
  }
  parsed
}
ncores <- max(1L, args$ncores)
blas_threads <- if (args$blas_threads <= 0L) {
  if (ncores > 1L) 1L else ncores
} else {
  max(1L, args$blas_threads)
}
weight_folds <- min(n, max(1L, args$weight_folds))
methods <- unique(trimws(strsplit(args$methods, ",", fixed = TRUE)[[1L]]))
parse_integer_list <- function(value, name) {
  parsed <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  if (!length(parsed) || anyNA(parsed) || any(parsed < 1L)) {
    stop(name, " must be a comma-separated list of positive integers.")
  }
  unique(parsed)
}
mgs_rank_candidates <- parse_integer_list(args$mgs_rank_candidates, "mgs_rank_candidates")
mgs_basis_candidates <- parse_integer_list(args$mgs_basis_candidates, "mgs_basis_candidates")
valid_methods <- c(
  "matfpca_refined", "matfpca_identity", "happ_pace", "happ_face",
  "mfaces", "mgsfpca", "shi_mdfpca", "spaco_minus",
  "sparse_product_fpca", "spline_pca"
)
unknown_methods <- setdiff(methods, valid_methods)
if (length(unknown_methods)) {
  stop("Unknown method(s): ", paste(unknown_methods, collapse = ", "),
       ". Valid choices are: ", paste(valid_methods, collapse = ", "))
}

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(blas_threads)
  RhpcBLASctl::omp_set_num_threads(blas_threads)
}

cat(
  "Simulation settings\n",
  "  seed: ", args$seed, "\n",
  "  n, p, q, nt: ", paste(c(n, p, q, nt), collapse = ", "), "\n",
  "  sampling ratio: ", paste(sampling_ratio, collapse = "--"), "\n",
  "  noise SD: ", sig, "\n",
  "  methods: ", paste(methods, collapse = ", "), "\n",
  "  method workers, BLAS threads: ", ncores, ", ", blas_threads, "\n",
  "  evaluation eigenspace rank: ", eval_ncomp, "\n",
  "  common covariance rank: ", if (is.null(cov_ncomp)) "full" else cov_ncomp, "\n",
  sep = ""
)

set.seed(args$seed)
generated <- gen_data_simple(
  n = n,
  p = p,
  q = q,
  npc_row = npc_row,
  npc_col = npc_col,
  npc_time = npc_time,
  nt = nt,
  t_range = c(0, 1),
  sampling_ratio = sampling_ratio,
  sig = sig,
  eigvecSet = args$eigset,
  eigvalSet = args$eigvalset,
  spatial_model = args$eigensystem,
  scores_distr = args$scores_distr,
  score_dependence = args$score_dependence,
  dependence_param = args$dependence_param
)
data <- generated$data
truth <- generated$cov

fit_one <- function(method) {
  cat("\n==== ", method, " ====\n", sep = "")
  fit <- NULL
  elapsed <- system.time({
    fit <- switch(
      method,
      matfpca_refined = fit_matfcov(
        data,
        wScheme = "gaussian_refined",
        noise = sig > 0,
        nbasis = nbasis,
        opwIter = 3L,
        temporal_method = args$temporal_method,
        nfolds = weight_folds,
        crossfit_diagonal = args$crossfit_diagonal
      ),
      matfpca_identity = fit_matfcov(
        data,
        wScheme = "identity",
        noise = sig > 0,
        nbasis = nbasis,
        temporal_method = args$temporal_method
      ),
      happ_pace = fit_happ_pace(
        data,
        npc = args$happ_npc,
        ncomp = cov_ncomp,
        ncores = ncores,
        backend = args$backend
      ),
      happ_face = fit_happ_face(
        data,
        npc = args$happ_npc,
        ncomp = cov_ncomp,
        knots = args$face_knots,
        pve = args$face_pve,
        ncores = ncores,
        backend = args$backend
      ),
      mfaces = fit_mfaces(
        data,
        ncomp = cov_ncomp,
        knots = args$mfaces_knots,
        pve = args$mfaces_pve,
        search_length = args$mfaces_search_length,
        max_outcomes = args$mfaces_max_outcomes,
        allow_large = as.logical(args$allow_slow_competitors)
      ),
      mgsfpca = fit_mgsfpca(
        data,
        ncomp = cov_ncomp,
        rank_candidates = mgs_rank_candidates,
        basis_candidates = mgs_basis_candidates,
        mu_nbasis = args$mgs_mu_nbasis,
        maxit = args$mgs_maxit
      ),
      shi_mdfpca = fit_shi_mdfpca(
        data,
        ncomp = if (is.null(cov_ncomp)) {
          max(eval_ncomp, args$shi_ncomp)
        } else {
          cov_ncomp
        },
        nbasis_time = args$shi_nbasis_time,
        nbasis_row = args$shi_nbasis_row,
        nbasis_col = args$shi_nbasis_col,
        maxit = args$shi_maxit,
        tolerance = args$shi_tolerance,
        smoothness = args$shi_smoothness,
        score_ridge = shi_score_ridge
      ),
      spaco_minus = fit_spaco_minus(
        data,
        ncomp = if (is.null(cov_ncomp)) {
          max(eval_ncomp, args$spaco_rank)
        } else {
          cov_ncomp
        },
        nbasis = args$spaco_nbasis,
        smoothness = args$spaco_smoothness,
        ridge = args$spaco_ridge,
        score_ridge = spaco_score_ridge,
        maxit = args$spaco_maxit
      ),
      sparse_product_fpca = fit_sparse_product_fpca(
        data,
        ncomp = cov_ncomp,
        spatial_npc = args$product_spatial_npc,
        temporal_npc = args$product_temporal_npc
      ),
      spline_pca = fit_spline_pca(
        data,
        nbasis = args$spline_nbasis,
        ncomp = cov_ncomp,
        smoothness = args$spline_smoothness
      )
    )
  })[["elapsed"]]
  metrics <- if (startsWith(method, "matfpca_")) {
    evaluate_product_fit(
      fit, truth, eval_ncomp = eval_ncomp, covariance_ncomp = cov_ncomp
    )
  } else {
    evaluate_full_eigenfit(fit, truth, eval_ncomp = eval_ncomp)
  }
  cat(
    "elapsed seconds: ", format(elapsed, digits = 5), "\n",
    "relative squared Hilbert--Schmidt covariance error: ",
    format(metrics[["rel_hs_cov"]], digits = 5), "\n",
    "mean sin^2 leading-subspace angle: ",
    format(metrics[["mean_sin2"]], digits = 5), "\n",
    "relative leading-eigenvalue SSE: ",
    format(metrics[["rel_eigenvalue_sse"]], digits = 5), "\n",
    sep = ""
  )
  list(fit = fit, elapsed = unname(elapsed), metrics = metrics)
}

results <- lapply(methods, function(method) {
  tryCatch(
    fit_one(method),
    error = function(error) {
      message(method, " failed: ", conditionMessage(error))
      list(fit = NULL, elapsed = NA_real_, metrics = NULL, error = conditionMessage(error))
    }
  )
})
names(results) <- methods

rows <- data.table::rbindlist(lapply(methods, function(method) {
  result <- results[[method]]
  metrics <- result$metrics
  data.table::data.table(
    method = method,
    status = if (is.null(result$error)) "ok" else "error",
    error = if (is.null(result$error)) NA_character_ else result$error,
    seed = args$seed,
    n = n,
    p = p,
    q = q,
    nt = nt,
    sampling_min = sampling_ratio[1L],
    sampling_max = sampling_ratio[2L],
    sig = sig,
    rank_row = npc_row,
    rank_col = npc_col,
    rank_time = npc_time,
    eigensystem = args$eigensystem,
    eigvec_set = args$eigset,
    eigval_set = args$eigvalset,
    score_distribution = args$scores_distr,
    score_dependence = args$score_dependence,
    dependence_param = generated$dependence_param,
    result_tag = if (nzchar(args$result_tag)) args$result_tag else NA_character_,
    eval_ncomp = eval_ncomp,
    cov_ncomp = if (is.null(cov_ncomp)) NA_integer_ else cov_ncomp,
    happ_npc = args$happ_npc,
    product_spatial_npc = args$product_spatial_npc,
    product_temporal_npc = args$product_temporal_npc,
    mgs_rank_candidates = paste(mgs_rank_candidates, collapse = ","),
    mgs_basis_candidates = paste(mgs_basis_candidates, collapse = ","),
    shi_nbasis_time = args$shi_nbasis_time,
    shi_nbasis_row = args$shi_nbasis_row,
    shi_nbasis_col = args$shi_nbasis_col,
    shi_smoothness = args$shi_smoothness,
    shi_score_ridge = as.character(shi_score_ridge),
    shi_selected_score_ridge = if (method == "shi_mdfpca" && !is.null(result$fit)) {
      result$fit$diagnostics$selected_score_ridge
    } else {
      NA_real_
    },
    spaco_smoothness = args$spaco_smoothness,
    spaco_ridge = args$spaco_ridge,
    spaco_score_ridge = as.character(spaco_score_ridge),
    spaco_selected_score_ridge = if (method == "spaco_minus" && !is.null(result$fit)) {
      result$fit$diagnostics$selected_score_ridge
    } else {
      NA_real_
    },
    spaco_minimum_cv_score_ridge = if (
      method == "spaco_minus" && !is.null(result$fit)
    ) {
      result$fit$diagnostics$minimum_cv_score_ridge
    } else {
      NA_real_
    },
    spaco_score_ridge_candidates = if (
      method == "spaco_minus" && !is.null(result$fit)
    ) {
      paste(result$fit$diagnostics$score_ridge_candidates, collapse = ",")
    } else {
      NA_character_
    },
    spaco_score_ridge_cv_error = if (
      method == "spaco_minus" && !is.null(result$fit)
    ) {
      paste(signif(result$fit$diagnostics$score_ridge_cv_error, 8L), collapse = ",")
    } else {
      NA_character_
    },
    spaco_score_ridge_cv_se = if (
      method == "spaco_minus" && !is.null(result$fit)
    ) {
      paste(signif(result$fit$diagnostics$score_ridge_cv_se, 8L), collapse = ",")
    } else {
      NA_character_
    },
    spaco_maxit = args$spaco_maxit,
    spline_nbasis = args$spline_nbasis,
    spline_smoothness = args$spline_smoothness,
    observed_visits_min = min(generated$visit_counts),
    observed_visits_mean = mean(generated$visit_counts),
    observed_visits_max = max(generated$visit_counts),
    failed_pixels = if (method %in% c("happ_pace", "happ_face") && !is.null(result$fit)) {
      pixel_status <- result$fit$diagnostics$pixel_status
      as.integer(sum(pixel_status[names(pixel_status) != "ok"]))
    } else {
      NA_integer_
    },
    rel_hs_cov = if (is.null(metrics)) NA_real_ else metrics[["rel_hs_cov"]],
    oracle_rank_rel_hs = if (is.null(metrics)) NA_real_ else metrics[["oracle_rank_rel_hs"]],
    excess_rel_hs_cov = if (is.null(metrics)) NA_real_ else metrics[["excess_rel_hs_cov"]],
    mean_sin2 = if (is.null(metrics)) NA_real_ else metrics[["mean_sin2"]],
    max_sin2 = if (is.null(metrics)) NA_real_ else metrics[["max_sin2"]],
    rel_eigenvalue_sse = if (is.null(metrics)) NA_real_ else metrics[["rel_eigenvalue_sse"]],
    evaluated_components = if (is.null(metrics)) NA_integer_ else as.integer(metrics[["evaluated_components"]]),
    estimated_rank = if (is.null(metrics)) NA_integer_ else as.integer(metrics[["estimated_rank"]]),
    sig2_hat = if (is.null(metrics)) NA_real_ else metrics[["sig2_hat"]],
    runtime = result$elapsed
  )
}), fill = TRUE)

print(rows[, .(
  method, status, rel_hs_cov, oracle_rank_rel_hs, excess_rel_hs_cov,
  mean_sin2, max_sin2,
  rel_eigenvalue_sse, estimated_rank, failed_pixels, runtime
)])

if (as.logical(args$save)) {
  dir.create(args$results_dir, recursive = TRUE, showWarnings = FALSE)
  token <- function(x) gsub("[^A-Za-z0-9]+", "p", format(x, trim = TRUE))
  file_id <- paste0(
    "n", n, "_p", p, "_q", q, "_nt", nt,
    "_s", s1, "-", s2, "_sig", token(sig),
    "_r", npc_row, "-", npc_col, "-", npc_time,
    "_eig", token(args$eigensystem),
    "_dist", token(args$scores_distr),
    "_dep", token(args$score_dependence),
    "_dp", token(generated$dependence_param),
    "_sd", args$seed
  )
  if (nzchar(args$result_tag)) {
    file_id <- paste0(file_id, "_tag", token(args$result_tag))
  }
  csv_path <- file.path(args$results_dir, paste0("methods_", file_id, ".csv"))
  temp_path <- tempfile(pattern = ".methods-", tmpdir = args$results_dir)
  data.table::fwrite(rows, temp_path)
  if (!file.rename(temp_path, csv_path)) {
    unlink(temp_path)
    stop("Could not publish result file: ", csv_path)
  }
  cat("Saved metrics to ", csv_path, "\n", sep = "")

  if (as.logical(args$save_fits)) {
    fit_path <- file.path(args$results_dir, paste0("fits_", file_id, ".rds"))
    saveRDS(list(results = results, generated = generated, arguments = args), fit_path)
    cat("Saved fits to ", fit_path, "\n", sep = "")
  }
}
