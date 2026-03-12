library(fda)
library(Matrix)
library(rTensor)
library(data.table)

# dev.off()
pa = par(no.readonly = TRUE)

source("matfdata.R")
source("matfpca.R")
source("plot_utils.R")
source(file.path("experiment", "gendata-simple.R"))

# Specify the number of cores to use
# num_cores <- 4
# cl <- parallel::makeCluster(num_cores)
# # Register the cluster with doParallel
# doParallel::registerDoParallel(cl)

num_cores <- RhpcBLASctl::blas_get_num_procs()
RhpcBLASctl::blas_set_num_threads(6 - 2)

# Data Generation ---------------------------------------------------------

# Generate matrix-valued functional data with
# a small number of "rank-one" eigenfunctions

# dimension and temporal space
n = 30 # n sample
p = 50  # dim row
q = 50  # dim column
nt = 50 + 1 # num time points
t_range = c(0,1)

# rank settings
# NOTE: set npc == dimension to ensure full-rank covariances
npc_row = p
# npc_row = 5
npc_col = q
# npc_col = 5
npc_time = 3

# sparse irregular sampling
sampling_ratio = c(0.05, 0.1)  # (Pmin, Pmax)
# sampling_ratio = c(0.2, 0.3)  # (Pmin, Pmax)
# sampling_ratio = c(0.8, 0.9)  # (Pmin, Pmax)

# data format:
# dimension 1: row
# dimension 2: column
# dimension 3: time

genData <- gen_data_simple(
  n, p, q,
  npc_row, npc_col, npc_time,
  nt, t_range,
  sig = 1e-1,
  sampling_ratio,
  eigvecSet = 2
)

data <- genData$data
cov_true <- genData$cov
eig_true <- genData$eig
# scores_true <- genData$scores
# X_true = recover_matfdata(eig_true$funcs, scores_true)

# Estimation --------------------------------------------------------------

cobj_idt <- fit_matfcov(
  data, wScheme = "identity", noise = TRUE, nbasis = 21, maxIter = 1
)
cov_idt <- cobj_idt$C

cobj_inv <- fit_matfcov(
  data, wScheme = "invcov", noise = TRUE, nbasis = 21,
  maxIter = 2
)
cov_inv <- cobj_inv$C

cobj_opt <- fit_matfcov(
  data, wScheme = "optimal", noise = TRUE, nbasis = 21,
  maxIter = 2, opwIter = 3
)
cov_opt <- cobj_opt$C

cat(
  ">>>>> Weighting scheme: identity\n",
  "Rel. MSE (U): ", matFNorm2(cov_idt$U - cov_true$U) / matFNorm2(cov_true$U), "\n",
  "Rel. MSE (V): ", matFNorm2(cov_idt$V - cov_true$V) / matFNorm2(cov_true$V), "\n",
  "Rel. MSE (K): ", fun2dNorm2(cov_idt$K - cov_true$K, cov_true$t_grid) / fun2dNorm2(cov_true$K, cov_true$t_grid), "\n",
  "Rel. MSE (C): ", matfcovDiffNorm2(cov_idt, cov_true) / matfcovNorm2(cov_true),
  sep = ""
)
cat(
  ">>>>> Weighting scheme: inv. cov.\n",
  "Rel. MSE (U): ", matFNorm2(cov_inv$U - cov_true$U) / matFNorm2(cov_true$U), "\n",
  "Rel. MSE (V): ", matFNorm2(cov_inv$V - cov_true$V) / matFNorm2(cov_true$V), "\n",
  "Rel. MSE (K): ", fun2dNorm2(cov_inv$K - cov_true$K, cov_true$t_grid) / fun2dNorm2(cov_true$K, cov_true$t_grid), "\n",
  "Rel. MSE (C): ", matfcovDiffNorm2(cov_inv, cov_true) / matfcovNorm2(cov_true),
  sep = ""
)
cat(
  ">>>>> Weighting scheme: optimal\n",
  "Rel. MSE (U): ", matFNorm2(cov_opt$U - cov_true$U) / matFNorm2(cov_true$U), "\n",
  "Rel. MSE (V): ", matFNorm2(cov_opt$V - cov_true$V) / matFNorm2(cov_true$V), "\n",
  "Rel. MSE (K): ", fun2dNorm2(cov_opt$K - cov_true$K, cov_true$t_grid) / fun2dNorm2(cov_true$K, cov_true$t_grid), "\n",
  "Rel. MSE (C): ", matfcovDiffNorm2(cov_opt, cov_true) / matfcovNorm2(cov_true),
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

# # evaluate eigen estimation
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
