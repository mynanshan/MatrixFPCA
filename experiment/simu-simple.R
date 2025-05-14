library(fda)
library(Matrix)
library(rTensor)
library(data.table)

pa = par(no.readonly = TRUE)

source("matfdata.R")
source("matfpca.R")
source("plot_utils.R")
source("experiment/gendata-simple.R")

# Specify the number of cores to use
# num_cores <- 4
# cl <- parallel::makeCluster(num_cores)
# # Register the cluster with doParallel
# doParallel::registerDoParallel(cl)

plan(multisession(workers = 4))

# Data Generation ---------------------------------------------------------

# Generate matrix-valued functional data with
# a small number of "rank-one" eigenfunctions

# dimension and temporal space
n = 10 # n sample
p = 100  # dim row
q = 100  # dim column
nt = 50 + 1 # num time points
t_range = c(0,1)

# rank settings
# NOTE: set npc == dimension to ensure full-rank covariances
npc_row = p
npc_col = q
npc_time = 5

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
  sampling_ratio
)

data <- genData$data
cov_true <- genData$cov
eig_true <- genData$eig
scores_true <- genData$scores
X_true = recover_matfdata(eig_true$funcs, scores_true)

# Estimation --------------------------------------------------------------

cov_fit = fit_matfcov(
  data, max_iter = 10, tol=1e-8, moment_type = "basic")
cov_hat = cov_fit$cov
cat("Ran for", cov_fit$iter, "iterations.")

# evaluate cov estimation
# mean((evalf.matfcov(cov_hat) - evalf.matfcov(cov_true))^2)
mean((scale_mat_trace(cov_hat$U) - scale_mat_trace(cov_true$U))^2)
mean((scale_mat_trace(cov_hat$V) - scale_mat_trace(cov_true$V))^2)
mean((scale_mat_trace(cov_hat$K) - scale_mat_trace(cov_true$K))^2)


par(mfrow=c(3,2), mar=c(2,1,1,1), oma=c(0,2,2,0))
tile_plot.matrix(cov_true$K)
mtext("K", side=2, line=1.5)
mtext("True", side=3, line=1)
tile_plot.matrix(cov_hat$K)
mtext("Estimate", side=3, line=1)
tile_plot.matrix(cov_true$U)
mtext("U", side=2, line=1.5)
tile_plot.matrix(cov_hat$U)
tile_plot.matrix(cov_true$V)
mtext("V", side=2, line=1.5)
tile_plot.matrix(cov_hat$V)
par(pa)

# evaluate eigen estimation
eig_obj = eigen.matfcov(cov_hat)
eig_obj = match_signs.matfeigen(eig_obj, eig_true)
eig_obj_trunc = trunc.matfeigen(eig_obj)
r_row = ncol(eig_obj_trunc$funcs$Phi)
r_col = ncol(eig_obj_trunc$funcs$Psi)
r_time = min(ncol(eig_obj_trunc$funcs$Rho), ncol(eig_true$funcs$Rho))

par(mfrow=c(3,2), mar=c(2,1,1,1), oma=c(0,2,2,0))
tile_plot.matrix(eig_true$funcs$Phi[,1:r_row,drop=F])
mtext(expression(Phi), side=2, line=1.5)
mtext("True", side=3, line=1)
tile_plot.matrix(eig_obj_trunc$funcs$Phi, xlab="t")
mtext("Estimate", side=3, line=1)
tile_plot.matrix(eig_true$funcs$Psi[,1:r_col,drop=F])
mtext(expression(Psi), side=2, line=1.5)
tile_plot.matrix(eig_obj_trunc$funcs$Psi, xlab="t")
tile_plot.matrix(eig_true$funcs$Rho[,1:r_time,drop=F])
mtext(expression(Rho), side=2, line=1.5)
tile_plot.matrix(eig_obj_trunc$funcs$Rho[,1:r_time,drop=F], xlab="t")
par(pa)

mean((eig_obj_trunc$funcs$Phi - eig_true$funcs$Phi[,1:r_row])^2)
mean((eig_obj_trunc$funcs$Psi - eig_true$funcs$Psi[,1:r_col])^2)
mean((eig_obj_trunc$funcs$Rho[,1:r_time] - eig_true$funcs$Rho[,1:r_time])^2)


# estimate scores and recover the process
scores_hat = get_matfpca_scores(data, eig_obj_trunc)
X_hat = recover_matfdata(eig_obj_trunc$funcs, scores_hat)
mean((X_hat@data - X_true@data)^2) / mean(X_true@data^2)

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

# sim studies:
# - other methods (baseline) 
# - different sample size
# - misspecified?
