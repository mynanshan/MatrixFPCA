#!<r-location>
#SBATCH --job-name=simu-simple        # Job name
#SBATCH --output=logs/simu_simple_%j.out         # Standard output file (%j expands to jobID)
#SBATCH --time=6:00:00                   # Maximum runtime (hh:mm:ss)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G                          # Memory allocation
library(fda)
library(Matrix)
library(rTensor)
library(argparse)
library(RhpcBLASctl)

RhpcBLASctl::blas_set_num_threads(4)

pa = par(no.readonly = TRUE)

source("matfdata.R")
source("matfpca.R")
source("plot_utils.R")
source("experiment/gendata-simple.R")

# Read command line arguments
parser <- argparse::ArgumentParser(
  description = "Determine the settings for this simple simulation."
)

# Add arguments
parser$add_argument("-n", "--n_subjects",
  type = "integer", default = 10,
  help = "Number of subjects"
)
parser$add_argument("-p", "--d_row",
  type = "integer", default = 10,
  help = "Dimension of rows"
)
parser$add_argument("-q", "--d_col",
  type = "integer", default = 10,
  help = "Dimension of columns"
)
parser$add_argument("-rt", "--rank_time",
  type = "integer", default = 5,
  help = "Rank of the time dimension"
)
parser$add_argument("-set", "--data_set",
  type = "integer", default = 1,
  help = "Which data generation settings."
)
parser$add_argument("-cm", "--compare_methods",
  type = "logical", default = FALSE,
  help = "Whether compare different FPCA methods"
)
parser$add_argument("--repeat",
  type = "integer", default = 100,
  help = "Number of repeat"
)

args <- parser$parse_args()

n <- args[["n_subjects"]]
p <- args[["d_row"]]
q <- args[["d_col"]]
npc_time = args[["rank_time"]]
eigset = args[["data_set"]]
compare_methods <- args[["compare_methods"]]
nrepeat <- args[["repeat"]]

dir_path <- file.path("experiment", "results")
file_id <- paste0(
  "_n",n,"_p",p,"_q",q,
  "_rt",npc_time,
  "_cm",as.integer(compare_methods)
)
file_name <- paste0("res-simu-simple", file_id, ".txt")
metafile_name <- paste0("meta-simu-simple", file_id, ".txt")

header <- c(
  "method",
  "n","p","q",
  "normU", "normV", "normK",
  "sumU", "sumV", "sumK",
  "maxU", "maxV", "maxK",
  "trU", "trV", "trK",
  "rrow", "rcol", "rtime",
  "rrow.trunc", "rcol.trunc", "rtime.trunc",
  "mseU", "mseV", "mseK",
  "relmseU", "relmseV", "relmseK",
  "maxdiffU", "maxdiffV", "maxdiffK",
  "msePhi", "msePsi", "mseRho",
  "mseX", "relmseX",
  "time"
)

cat(header, sep = "\n", file = file.path(dir_path, metafile_name))

# Data Generation ---------------------------------------------------------

for (ir in 1:nrepeat) {
cat("Repeat", ir, "\n")
# dimension and temporal space
# n, p, q
nt = 50 + 1 # num time points
t_range = c(0,1)

# rank settings
npc_row = p
npc_col = q
# npc_time

# sparse irregular sampling
sampling_ratio = c(0.05, 0.1)  # (Pmin, Pmax)

genData <- gen_data_simple(
  n, p, q,
  npc_row, npc_col, npc_time,
  nt, t_range,
  sampling_ratio,
  eigset = eigset
)

data <- genData$data
cov_true <- genData$cov
eig_true <- genData$eig
scores_true <- genData$scores
X_true = recover_matfdata(eig_true$funcs, scores_true)

# Mat-FPCA ----------------------------------------------------------------

time_start = Sys.time()

# fit covariances
cov_hat = fit_matfcov(data, wScheme = "OBS")

# select ranks
eig_obj = eigen.matfcov(cov_hat)
eig_obj = match_signs.matfeigen(eig_obj, eig_true)
eig_obj_trunc = trunc.matfeigen(eig_obj)
r_row = min(ncol(eig_obj_trunc$funcs$Phi), ncol(eig_true$funcs$Phi))
r_col = min(ncol(eig_obj_trunc$funcs$Psi), ncol(eig_true$funcs$Psi))
r_time = min(ncol(eig_obj_trunc$funcs$Rho), ncol(eig_true$funcs$Rho))

# estimates FPC scores
scores_hat = get_matfpca_scores(data, eig_obj_trunc)

# recover the continuous data
X_hat = recover_matfdata(eig_obj_trunc$funcs, scores_hat)

time_end = Sys.time()

mseU <- mean((cov_hat$U - cov_true$U)^2)
relmseU <- mseU / mean((cov_true$U)^2)
maxdiffU <- max(abs(cov_hat$U - cov_true$U))
mseV <- mean((cov_hat$V - cov_true$V)^2)
relmseV <- mseV / mean((cov_true$V)^2)
maxdiffV <- max(abs(cov_hat$V - cov_true$V))
mseK <- mean((cov_hat$K - cov_true$K)^2)
relmseK <- mseK / mean((cov_true$K)^2)
maxdiffK <- max(abs(cov_hat$K - cov_true$K))

output <- c(
  # settings
  "MatFPCA",
  n, p, q,
  sqrt(sum((cov_true$U)^2)),
  sqrt(sum((cov_true$V)^2)),
  sqrt(sum((cov_true$K)^2)),
  sum(cov_true$U),
  sum(cov_true$V),
  sum(cov_true$K),
  max(abs(cov_true$U)),
  max(abs(cov_true$V)),
  max(abs(cov_true$K)),
  sum(diag(cov_true$U)),
  sum(diag(cov_true$V)),
  sum(diag(cov_true$K)),
  npc_row, npc_col, npc_time,
  r_row, r_col, r_time,
  # model performance
  mseU |> round(5),
  mseV |> round(5),
  mseK |> round(5),
  relmseU |> round(5),
  relmseV |> round(5),
  relmseK |> round(5),
  maxdiffU |> round(5),
  maxdiffV |> round(5),
  maxdiffK |> round(5),
  mean((eig_obj_trunc$funcs$Phi - eig_true$funcs$Phi[,1:r_row])^2) |> round(5),
  mean((eig_obj_trunc$funcs$Psi - eig_true$funcs$Psi[,1:r_col])^2) |> round(5),
  mean((eig_obj_trunc$funcs$Rho[,1:r_time] - eig_true$funcs$Rho[,1:r_time])^2) |> round(5),
  mean((X_hat@data - X_true@data)^2) |> round(5),
  (mean((X_hat@data - X_true@data)^2) / mean(X_true@data^2)) |> round(5),
  as.numeric(time_end - time_start) |> round(3)
)

write(paste0(output, collapse = ","),
  file = file.path(dir_path, file_name), append = TRUE)

if (compare_methods) {
  # TODO
  
  # Chiou's Multivariate FPCA -----------------------------------------------
  # TODO: look for codes of this method
  
  
  # Happ's MFPCA ------------------------------------------------------------
  
  # convert data to Happ's funData format
  # t_list <- get_t_list(data)
  # indices <- expand.grid(i=1:p, j=1:q)
  # funDataList <- mapply(
  #   \(i, j) {
  #     funData::irregFunData(
  #       t_list, lapply(data$X_list, \(X) X@data[i,j,]))
  #   }, indices$i, indices$j
  # )
  # mfunData <- funData::multiFunData(funDataList)
  
  
  # Haixu's MFPC ------------------------------------------------------------
  
}
  
}

