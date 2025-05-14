library(fda)
library(Matrix)
library(rTensor)
library(data.table)
library(lubridate)

pa = par(no.readonly = TRUE)

source("matfdata.R")
source("matfpca.R")
source("plot_utils.R")

data_dir = "../data/ADNI/ADNIpng/axial"
obs_names = list.dirs(data_dir, full.names = FALSE)
obs_names = obs_names[obs_names != ""]
obs_name_split = strsplit(obs_names, "_")
subjs = sapply(obs_name_split, \(sp) sp[1])
time_extractor = function(time_label) {
  if (time_label == "bl") {
    return(0)
  } else {
    return(as.integer(stringr::str_remove(time_label, "m")))
  }
}
time_labels = sapply(obs_name_split, \(sp) sp[2])
times = sapply(obs_name_split, \(sp) time_extractor(sp[2]))
subj_names = sort(as.integer(unique(subjs)))
subj_times = lapply(subj_names, \(subj) {
  return(times[subjs == as.character(subj)])
})
subj_time_labels = lapply(subj_names, \(subj) {
  return(time_labels[subjs == as.character(subj)])
})
for (i in seq_along(subj_names)) {
  ord = order(subj_times[[i]])
  subj_times[[i]] = subj_times[[i]][ord]
  subj_time_labels[[i]] = subj_time_labels[[i]][ord]
}
# get_file_path = function(subj, time_label) {
#   file.path(data_dir,
#     paste0(subj, "_", time_label),
#     "100.png")
# }
# save_subj_rds = function(subj) {
#   i = which(subj_names == subj)
#   times = subj_times[[i]]
#   tlabs = subj_time_labels[[i]]
#   m = length(tlabs)
#   dims = dim(png::readPNG(get_file_path(subj, tlabs[1])))
#   dat = sapply(tlabs, \(ti) png::readPNG(get_file_path(subj, ti))) |>
#     array(dim = c(dims, m))
#   saveRDS(times, file.path(data_dir, paste0(subj, "_", "times.rds")))
#   saveRDS(dat, file.path(data_dir, paste0(subj, "_", "img.rds")))
# }
# for (subj in subj_names){
#   print(subj)
#   save_subj_rds(subj)
# }


n = length(subj_names)
X_list = lapply(1:n, \(i) readRDS(file.path(
  data_dir, paste0(subj_names[i], "_img.rds"))))
t_list = lapply(1:n, \(i) readRDS(file.path(
  data_dir, paste0(subj_names[i], "_times.rds"))))
t_grid = sort(unique(unlist(t_list)))

# rescaling
for (i in seq_len(n)) {
  s = max(X_list[[i]]) - min(X_list[[i]])
  X_list[[i]] = (X_list[[i]] - min(X_list[[i]])) / s
}

# calculate mean
X_mean = array(0, dim = c(dim(X_list[[1]][,,1]), length(t_grid)))
for (j in seq_along(t_grid)) {
  nj = sum(sapply(t_list, \(ti) t_grid[j] %in% ti))
  for (i in seq_len(n)) {
    jid = which(t_list[[i]] == t_grid[j])
    if (length(jid) == 0) {
      next
    } else {
      X_mean[,,j] = X_mean[,,j] + X_list[[i]][,,jid]
    }
  }
  X_mean[,,j] = X_mean[,,j] / nj
}
tile_plot.matf(as.tensor(X_mean), paste0("t=", t_grid))

# de-mean
for (i in seq_len(n)) {
  jid = which(t_list[[i]] %in% t_grid)
  X_list[[i]] = X_list[[i]] - X_mean[,,jid,drop=F]
}
tile_plot.matf(as.tensor(X_list[[6]]), paste0("t=", t_list[[i]]))


# main algo -----------------------------------------------------

t_idx_list = lapply(
  seq_len(n), \(i) which(t_list[[i]] %in% t_grid)
)

data = create_matrix_funcdata_irreg(
  X_list = lapply(X_list, \(X) as.tensor(X)),
  t_grid = t_grid,
  t_idx_list = t_idx_list
)

rm(X_list)


p = dim(data$X_list[[1]])[1]
q = dim(data$X_list[[1]])[2]
nt = length(t_grid)

# cov_fit = fit_matfcov(
#   data, max_iter = 2, tol=1e-8, nugget = c(1e-4, 1e-4, 1e-4))
# saveRDS(cov_fit, "adni_covfit.rds")
cov_fit = readRDS("adni_covfit.rds")
cov_hat = cov_fit$cov

eig_obj = eigen.matfcov(cov_hat)
eig_obj_trunc = trunc.matfeigen(eig_obj)

scores_hat = readRDS("adni_scores.rds")
# scores_hat = get_matfpca_scores(data, eig_obj_trunc, nugget = 1e-1)
# saveRDS(scores_hat, "adni_scores.rds")
X_hat = readRDS("adni_Xhat.rds")
# X_hat = recover_matfdata(eig_obj_trunc$funcs, scores_hat)
# saveRDS(X_hat, "adni_Xhat.rds")

i = sample(1:n, 1)  # i = 2, 5
tile_plot.matf(data$X_list[[i]], paste0("t=", t_grid[t_idx_list[[i]]]))
tile_plot.matf(X_hat[i,,,t_idx_list[[i]]], paste0("t=", t_grid[t_idx_list[[i]]]))

eigF = evalf.matfeigf(
  .trunc.matfeigen(eig_obj_trunc, 3, 3, 3)$funcs,
  t_idx = t_idx_list[[i]])
for (r in 1:27) {
  tile_plot.matf(
    as.tensor(array(eigF[,r], dim=c(data$vdim, length(t_idx_list[[i]])))),
    paste0("t=", t_grid[t_idx_list[[i]]]))
}


# dynamic FPCA, pointwise algo, comparison
# not using low-rank info, would it be noisy?
# pixel-wise: FPCA
# time-point FPCA and then smooth
