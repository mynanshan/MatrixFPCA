#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript

# Read command line arguments
parser <- argparse::ArgumentParser(
  description = "Determine the settings for this simple simulation."
)

# Add arguments
parser$add_argument("-n", "--n_subjects",
  type = "integer", default = 20,
  help = "Number of subjects"
)
parser$add_argument("-nt", "--n_time_points",
  type = "integer", default = 51,
  help = "Number of grid points for time"
)
parser$add_argument("-t0", "--t0",
  type = "numeric", default = 0,
  help = "Start of the time domain"
)
parser$add_argument("-t1", "--t1",
  type = "numeric", default = 1,
  help = "End of the time domain"
)
parser$add_argument("-s1", "--s_obs_min",
  type = "integer", default = 1,
  help = "Min sampling ratio, 1 unit for 0.05."
)
parser$add_argument("-s2", "--s_obs_max",
  type = "integer", default = 3,
  help = "Max sampling ratio, 1 unit for 0.05."
)
parser$add_argument("-p", "--d_row",
  type = "integer", default = 30,
  help = "Dimension of rows"
)
parser$add_argument("-q", "--d_col",
  type = "integer", default = 30,
  help = "Dimension of columns"
)
parser$add_argument("-sig", "--noise_lvl",
  type = "numeric", default = 0,
  help = "Noise level sigma"
)
parser$add_argument("-rr", "--rank_row",
  type = "integer", default = 0,
  help = "Rank of the row dimension (0: full rank)"
)
parser$add_argument("-rc", "--rank_col",
  type = "integer", default = 0,
  help = "Rank of the column dimension (0: full rank)"
)
parser$add_argument("-rt", "--rank_time",
  type = "integer", default = 5,
  help = "Rank of the time dimension"
)
parser$add_argument("-nb", "--nbasis",
  type = "integer", default = 0,
  help = "Marginal basis dimension (0: theory-guided automatic choice)."
)
parser$add_argument("-set", "--eigset",
  type = "integer", default = 2,
  help = "Which data generation settings."
)
parser$add_argument("-scd", "--scores_distr",
  type = "character", default = "normal",
  choices = c("normal", "uniform", "beta32", "t5"),
  help = "Distribution of the FPC scores."
)
parser$add_argument("-dep", "--score_dependence",
  type = "character", default = "independent",
  choices = c(
    "independent",
    "common_amplitude",
    "common_signed",
    "heavy_amplitude"
  ),
  help = "Dependence structure of standardized FPC scores."
)
parser$add_argument("-dpar", "--dependence_param",
  type = "numeric", default = 1,
  help = paste(
    "Var(A_i^2) for common_amplitude;",
    "Pareto tail index for heavy_amplitude."
  )
)
parser$add_argument("--seed",
  type = "integer", default = -1,
  help = "The base seed"
)
parser$add_argument("--save",
  type = "integer", default = ifelse(interactive(), 0, 1),
  help = "Whether to save the results"
)
parser$add_argument("--results_dir",
  type = "character", default = file.path("experiment", "results"),
  help = "Directory for one-file-per-configuration simulation results."
)
parser$add_argument("--manual",
  type = "integer", default = ifelse(interactive(), 1, 0),
  help = "Whether allow manual change of settings."
)

args <- parser$parse_args()

n <- max(args[["n_subjects"]], 1)
p <- max(args[["d_row"]], 1)
q <- max(args[["d_col"]], 1)
nt <- max(args[["n_time_points"]], 11)
t0 <- args[["t0"]]
t1 <- max(args[["t1"]], t0 + 1e-3)
t_range <- c(t0, t1)

s_unit <- 0.05
s1 <- max(args[["s_obs_min"]], 1)
s2 <- min(max(args[["s_obs_max"]], s1 + 1), floor(1 / s_unit))
sr1 <- s1 * s_unit
sr2 <- s2 * s_unit
sampling_ratio = c(sr1, sr2)  # (Pmin, Pmax)

sig <- max(args[["noise_lvl"]], 0)

rr <- args[["rank_row"]]
rc <- args[["rank_col"]]
rt <- args[["rank_time"]]

npc_row <- if (rr == 0) p else min(max(rr, 1), p)
npc_col <- if (rc == 0) q else min(max(rc, 1), q)
npc_time <- min(max(rt, 1), nt)

nbasis <- if (args[["nbasis"]] <= 0) "theory" else
  min(max(args[["nbasis"]], 3), nt)
eigset <- args[["eigset"]]
scores_distr <- args[["scores_distr"]]
score_dependence <- args[["score_dependence"]]
dependence_param <- args[["dependence_param"]]
if (score_dependence %in% c("independent", "common_signed")) {
  dependence_param <- 0
}
seed <- args[["seed"]]
seed <- ifelse(seed < 0, ceiling(runif(1)*9999), seed)
ifsave <- as.logical(args[["save"]])
results_dir <- args[["results_dir"]]

manual_mode <- as.logical(args$manual)

if (manual_mode) {
  scores_distr <- "beta32"
  score_dependence <- "common_signed"
  dependence_param <- 0
}

cat(
  ">>>>>>>>>>> Settings <<<<<<<<<<< \n",
  "seed =", seed, "\n",
  "n =", n, "\n",
  "p =", p, "\n",
  "q =", q, "\n",
  "sig =", sig, "\n",
  "nt =", nt, "\n",
  "nbasis=", nbasis, "\n",
  "sr1 =", sr1, "\n",
  "sr2 =", sr2, "\n",
  "npc_row=", npc_row, "\n",
  "npc_col=", npc_col, "\n",
  "npc_time=", npc_time, "\n",
  "scores_distr=", scores_distr, "\n",
  "score_dependence=", score_dependence, "\n",
  "dependence_param=", dependence_param, "\n"
)

library(fda)
library(Matrix)
library(rTensor)
library(data.table)

# dev.off()
# pa = par(no.readonly = TRUE)

# Specify the number of cores to use
is_slurm <- nzchar(Sys.getenv("SLURM_JOB_ID"))
if (is_slurm) {
  ncores_all <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
} else {
  ncores_all <- RhpcBLASctl::get_num_procs()
}
RhpcBLASctl::blas_set_num_threads(ncores_all)

source("matfdata.R")
source("matfpca.R")
source("plot_utils.R")
source(file.path("experiment", "gendata-simple.R"))

# Data Generation ---------------------------------------------------------

set.seed(seed)

# data format:
# dimension 1: row
# dimension 2: column
# dimension 3: time

genData <- gen_data_simple(
  n, p, q,
  npc_row, npc_col, npc_time,
  nt, t_range,
  sampling_ratio = sampling_ratio,
  sig = sig,
  eigvecSet = eigset,
  scores_distr = scores_distr,
  score_dependence = score_dependence,
  dependence_param = dependence_param
)

data <- genData$data
cov_true <- genData$cov
eig_true <- genData$eig
# scores_true <- genData$scores
# X_true = recover_matfdata(eig_true$funcs, scores_true)

# Estimation --------------------------------------------------------------

# The heuristic inverse-covariance schemes remain available in fit_matfcov(),
# but paper simulations compare only the two estimators defined in the theory.
weight_schemes <- c("identity", "optimal")

evaluate_fit <- function(fit, truth) {
  estimate <- fit$C
  sseU <- matFNorm2(estimate$U - truth$U)
  sseV <- matFNorm2(estimate$V - truth$V)
  iseK <- fun2dNorm2(estimate$K - truth$K, truth$t_grid)
  siseC <- matfcovDiffNorm2(estimate, truth)
  c(
    sseU = sseU,
    sseV = sseV,
    iseK = iseK,
    siseC = siseC,
    relsseU = sseU / matFNorm2(truth$U),
    relsseV = sseV / matFNorm2(truth$V),
    reliseK = iseK / fun2dNorm2(truth$K, truth$t_grid),
    relsiseC = siseC / matfcovNorm2(truth)
  )
}

fit_results <- lapply(weight_schemes, function(scheme) {
  cat("==== Start Estimation:", scheme, "weights ====\n")
  fit <- fit_matfcov(
    data,
    wScheme = scheme,
    noise = TRUE,
    nbasis = nbasis,
    opwIter = 3
  )
  metrics <- evaluate_fit(fit, cov_true)
  cat(
    ">>>>> Weighting scheme: ", scheme, "\n",
    "Rel. SSE (U): ", metrics[["relsseU"]], "\n",
    "Rel. SSE (V): ", metrics[["relsseV"]], "\n",
    "Rel. ISE (K): ", metrics[["reliseK"]], "\n",
    "Rel. SISE (C): ", metrics[["relsiseC"]], "\n",
    sep = ""
  )
  list(fit = fit, metrics = metrics)
})
names(fit_results) <- weight_schemes


# par(mfrow=c(3,2), mar=c(2,1,1,1), oma=c(0,2,2,0))
# rangeK <- range(c(cov_true$K, cov_hat$K))
# tile_plot.matrix(cov_true$K, zlim=rangeK)
# mtext("K", side=2, line=1.5)
# mtext("True", side=3, line=1)
# tile_plot.matrix(cov_hat$K)
# mtext("Estimate", side=3, line=1)
# rangeU <- range(c(cov_true$U, cov_hat$U))
# tile_plot.matrix(cov_true$U, zlim=rangeU)
# mtext("U", side=2, line=1.5)
# tile_plot.matrix(cov_hat$U)
# rangeV <- range(c(cov_true$V, cov_hat$V))
# tile_plot.matrix(cov_true$V, zlim=rangeV)
# mtext("V", side=2, line=1.5)
# tile_plot.matrix(cov_hat$V)
# par(pa)

# evaluate eigen estimation
# eig_obj = eigen.matfcov(cov_hat)
# eig_obj = match_signs.matfeigen(eig_obj, eig_true)
# eig_obj_trunc = trunc.matfeigen(eig_obj)
# r_row = ncol(eig_obj_trunc$funcs$Phi)
# r_col = ncol(eig_obj_trunc$funcs$Psi)
# r_time = min(ncol(eig_obj_trunc$funcs$Theta), ncol(eig_true$funcs$Theta))

# par(mfrow=c(3,2), mar=c(2,1,1,1), oma=c(0,2,2,0))
# tile_plot.matrix(eig_true$funcs$Phi[,1:r_row,drop=F])
# mtext(expression(Phi), side=2, line=1.5)
# mtext("True", side=3, line=1)
# tile_plot.matrix(eig_obj_trunc$funcs$Phi, xlab="t")
# mtext("Estimate", side=3, line=1)
# tile_plot.matrix(eig_true$funcs$Psi[,1:r_col,drop=F])
# mtext(expression(Psi), side=2, line=1.5)
# tile_plot.matrix(eig_obj_trunc$funcs$Psi, xlab="t")
# tile_plot.matrix(eig_true$funcs$Theta[,1:r_time,drop=F])
# mtext(expression(Theta), side=2, line=1.5)
# tile_plot.matrix(eig_obj_trunc$funcs$Theta[,1:r_time,drop=F], xlab="t")
# par(pa)

# mean((eig_obj_trunc$funcs$Phi - eig_true$funcs$Phi[,1:r_row])^2)
# mean((eig_obj_trunc$funcs$Psi - eig_true$funcs$Psi[,1:r_col])^2)
# mean((eig_obj_trunc$funcs$Theta[,1:r_time] - eig_true$funcs$Theta[,1:r_time])^2)


# estimate scores and recover the process
# scores_hat = get_matfpca_scores(data, eig_obj_trunc)
# X_hat = recover_matfdata(eig_obj_trunc$funcs, scores_hat)
# mean((X_hat@data - X_true@data)^2) / mean(X_true@data^2)

# i = sample(1:n, 1)
# par(mfrow=c(2,1), mar=c(2,1,1,1), oma=c(0,2,2,0))
# tile_plot.matf(X_hat[i,,,t_idx_list[[i]]], paste0("t=", t_grid[t_idx_list[[i]]]))
# tile_plot.matf(X[i,,,t_idx_list[[i]]], paste0("t=", t_grid[t_idx_list[[i]]]))
# par(pa)

# par(mfrow=c(2,1), mar=c(2,1,1,1), oma=c(0,2,2,0))
# tile_plot.matf(X_hat[i,,,seq(1,nt,5)], paste0("t=", t_grid[seq(1,nt,5)]),
#   label_cex = 1.)
# tile_plot.matf(X[i,,,seq(1,nt,5)], paste0("t=", t_grid[seq(1,nt,5)]), label_cex = 1.)
# par(pa)

# j_list = seq(1, p, 2)
# k_list = seq(1, q, 2)
# i = sample(1:n, 1)
# par(mfrow=c(length(j_list),length(k_list)), mar=c(2,1,1,1), oma=c(0,2,2,0))
# for (j in j_list) {
#   for (k in k_list) {
#     plot(
#       t_grid, X[i,j,k,]@data, type="l",
#       ylim = range(c(X[i,j,k,]@data, X_hat[i,j,k,]@data))
#     )
#     lines(t_grid, X_hat[i,j,k,]@data, col=2, lty=2)
#     points(t_grid[data$t_idx_list[[i]]], data$X_list[[i]][j,k,]@data)
#   }
# }
# par(pa)


# Save results -----------------------------------------------------------

if (ifsave) {

  dir_path <- results_dir
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)

  file_token <- function(x) {
    gsub("[^A-Za-z0-9]+", "p", format(x, scientific = FALSE, trim = TRUE))
  }
  file_id <- paste0(
    "_n", n,
    "_p", p,
    "_sig", file_token(sig),
    "_d", scores_distr,
    "_dep", score_dependence,
    "_dp", file_token(dependence_param),
    "_s", s1, "-", s2,
    "_nt", nt,
    "_r", npc_row, "-", npc_col, "-", npc_time,
    "_set", eigset,
    "_sd", seed
  )
  file_name <- paste0("simu", file_id, ".csv")
  file_path <- file.path(dir_path, file_name)

  header <- c(
    "method", "seed", "weight",
    "n", "p", "q", "sig", "distr", "dependence", "dependence_param",
    "s1", "s2", "nt",
    "sseU", "sseV", "iseK", "siseC",
    "relsseU", "relsseV", "reliseK", "relsiseC"
  )

  formatter <- function(x, digits = 7, logarithm = TRUE) {
    if (logarithm) return(round(log(x + 10 * .Machine$double.eps), digits))
    round(x, digits)
  }

  output_rows <- lapply(weight_schemes, function(scheme) {
    metrics <- fit_results[[scheme]]$metrics
    as.list(stats::setNames(c(
      "MatFPCA", seed, scheme,
      n, p, q, sig, scores_distr, score_dependence, dependence_param,
      s1, s2, nt,
      formatter(metrics[["sseU"]]),
      formatter(metrics[["sseV"]]),
      formatter(metrics[["iseK"]]),
      formatter(metrics[["siseC"]]),
      formatter(metrics[["relsseU"]]),
      formatter(metrics[["relsseV"]]),
      formatter(metrics[["reliseK"]]),
      formatter(metrics[["relsiseC"]])
    ), header))
  })
  data.table::fwrite(
    data.table::rbindlist(output_rows),
    file = file_path
  )

}
