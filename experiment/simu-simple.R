#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=simu-simple
#SBATCH --output=logs/simu_simple_%j.out
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G

# Read command line arguments
parser <- argparse::ArgumentParser(
  description = "Determine the settings for this simple simulation."
)

# Add arguments
parser$add_argument("-n", "--n_subjects",
  type = "integer", default = 30,
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
  type = "integer", default = 50,
  help = "Dimension of rows"
)
parser$add_argument("-q", "--d_col",
  type = "integer", default = 50,
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
  type = "integer", default = 21,
  help = "Number of basis functions."
)
parser$add_argument("-set", "--eigset",
  type = "integer", default = 2,
  help = "Which data generation settings."
)
parser$add_argument("--seed",
  type = "integer", default = -1,
  help = "The base seed"
)
parser$add_argument("--save",
  type = "integer", default = ifelse(interactive(), 0, 1),
  help = "Whether to save the results"
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

nbasis <- min(max(args[["nbasis"]], 1), nt)
eigset <- args[["eigset"]]
seed <- args[["seed"]]
seed <- ifelse(seed < 0, ceiling(runif(1)*9999), seed)
ifsave <- as.logical(args[["save"]])

cat(
  ">>>>>>>>>>> Settings <<<<<<<<<<< \n",
  "seed=", seed, "\n",
  "n=", n, "\n",
  "p=", p, "\n",
  "q=", q, "\n",
  "npc_time=", npc_time, "\n"
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
  eigvecSet = eigset
)

data <- genData$data
cov_true <- genData$cov
eig_true <- genData$eig
# scores_true <- genData$scores
# X_true = recover_matfdata(eig_true$funcs, scores_true)

# Estimation --------------------------------------------------------------

cat("==== Start Estimation: Identity weights ====\n")

cobj_idt <- fit_matfcov(
  data, wScheme = "identity", noise = TRUE, nbasis = nbasis, maxIter = 1
)
cov_idt <- cobj_idt$C

cat("==== Start Estimation: Inv.cov. weights ====\n")

cobj_inv <- fit_matfcov(
  data, wScheme = "invcov", noise = TRUE, nbasis = nbasis,
  maxIter = 2
)
cov_inv <- cobj_inv$C

cat("==== Start Estimation: Optimal weights ====\n")

cobj_opt <- fit_matfcov(
  data, wScheme = "optimal", noise = TRUE, nbasis = nbasis,
  maxIter = 2, opwIter = 3
)
cov_opt <- cobj_opt$C

sseU_idt <- matFNorm2(cov_idt$U - cov_true$U)
sseV_idt <- matFNorm2(cov_idt$V - cov_true$V)
iseK_idt <- fun2dNorm2(cov_idt$K - cov_true$K, cov_true$t_grid)
siseC_idt <- matfcovDiffNorm2(cov_idt, cov_true)
relsseU_idt <- sseU_idt / matFNorm2(cov_true$U)
relsseV_idt <- sseV_idt / matFNorm2(cov_true$V)
reliseK_idt <- iseK_idt / fun2dNorm2(cov_true$K, cov_true$t_grid)
relsiseC_idt <- siseC_idt / matfcovNorm2(cov_true)

sseU_inv <- matFNorm2(cov_inv$U - cov_true$U)
sseV_inv <- matFNorm2(cov_inv$V - cov_true$V)
iseK_inv <- fun2dNorm2(cov_inv$K - cov_true$K, cov_true$t_grid)
siseC_inv <- matfcovDiffNorm2(cov_inv, cov_true)
relsseU_inv <- sseU_inv / matFNorm2(cov_true$U)
relsseV_inv <- sseV_inv / matFNorm2(cov_true$V)
reliseK_inv <- iseK_inv / fun2dNorm2(cov_true$K, cov_true$t_grid)
relsiseC_inv <- siseC_inv / matfcovNorm2(cov_true)

sseU_opt <- matFNorm2(cov_opt$U - cov_true$U)
sseV_opt <- matFNorm2(cov_opt$V - cov_true$V)
iseK_opt <- fun2dNorm2(cov_opt$K - cov_true$K, cov_true$t_grid)
siseC_opt <- matfcovDiffNorm2(cov_opt, cov_true)
relsseU_opt <- sseU_opt / matFNorm2(cov_true$U)
relsseV_opt <- sseV_opt / matFNorm2(cov_true$V)
reliseK_opt <- iseK_opt / fun2dNorm2(cov_true$K, cov_true$t_grid)
relsiseC_opt <- siseC_opt / matfcovNorm2(cov_true)

cat(
  ">>>>> Weighting scheme: identity\n",
  "Rel. SSE (U): ", relsseU_idt, "\n",
  "Rel. SSE (V): ", relsseV_idt, "\n",
  "Rel. ISE (K): ", reliseK_idt, "\n",
  "Rel. SISE (C): ", relsiseC_idt, "\n",
  sep = ""
)
cat(
  ">>>>> Weighting scheme: invcov\n",
  "Rel. SSE (U): ", relsseU_inv, "\n",
  "Rel. SSE (V): ", relsseV_inv, "\n",
  "Rel. ISE (K): ", reliseK_inv, "\n",
  "Rel. SISE (C): ", relsiseC_inv, "\n",
  sep = ""
)
cat(
  ">>>>> Weighting scheme: optimal\n",
  "Rel. SSE (U): ", relsseU_opt, "\n",
  "Rel. SSE (V): ", relsseV_opt, "\n",
  "Rel. ISE (K): ", reliseK_opt, "\n",
  "Rel. SISE (C): ", relsiseC_opt, "\n",
  sep = ""
)


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

  dir_path <- file.path("experiment", "results")
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
  file_id <- paste0(
    "_r", npc_row, "-", npc_col, "-", npc_time, "_set", eigset, "_sd", seed
  )
  file_name <- paste0("simu", file_id, ".txt")
  file_path <- file.path(dir_path, file_name)

  header <- c(
    "method", "seed", "weight",
    "n", "p", "q", "sig", "s1", "s2", "nt",
    "sseU", "sseV", "iseK", "siseC",
    "relsseU", "relsseV", "reliseK", "relsiseC"
  )

  file_info <- file.info(file_path)
    
  # If empty, write the predefined text
  if (is.na(file_info$size) || file_info$size == 0) {
    writeLines(paste0(header, collapse=","), con = file_path)
  }

  formatter <- function(x, digits = 7, logarithm = TRUE) {
    if (logarithm) return(round(log(x + 10 * .Machine$double.eps), digits))
    round(x, digits)
  }

  out_idt <- c(
    "MatFPCA", seed, "identity", 
    n, p, q, sig, s1, s2, nt,
    formatter(sseU_idt), formatter(sseV_idt), formatter(iseK_idt), formatter(siseC_idt),
    formatter(relsseU_idt), formatter(relsseV_idt), formatter(reliseK_idt), formatter(relsiseC_idt)
  )
  out_inv <- c(
    "MatFPCA", seed, "invcov", 
    n, p, q, sig, s1, s2, nt,
    formatter(sseU_inv), formatter(sseV_inv), formatter(iseK_inv), formatter(siseC_inv),
    formatter(relsseU_inv), formatter(relsseV_inv), formatter(reliseK_inv), formatter(relsiseC_inv)
  )
  out_opt <- c(
    "MatFPCA", seed, "optimal", 
    n, p, q, sig, s1, s2, nt,
    formatter(sseU_opt), formatter(sseV_opt), formatter(iseK_opt), formatter(siseC_opt),
    formatter(relsseU_opt), formatter(relsseV_opt), formatter(reliseK_opt), formatter(relsiseC_opt)
  )

  write(paste0(out_idt, collapse = ","), file = file_path, append = TRUE)
  write(paste0(out_inv, collapse = ","), file = file_path, append = TRUE)
  write(paste0(out_opt, collapse = ","), file = file_path, append = TRUE)

}
