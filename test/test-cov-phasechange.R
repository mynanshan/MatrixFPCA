#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=test-cov
#SBATCH --output=logs/test-cov_%j.out
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G

library(fda)
library(mgcv)
library(foreach)
library(tibble)

# Read command line arguments
parser <- argparse::ArgumentParser(
  description = "Determine the settings for this simple simulation."
)
parser$add_argument("-n", "--n_subjects",
  type = "integer", default = 50,
  help = "Number of subjects"
)
parser$add_argument("--seed",
  type = "integer", default = 1,
  help = "seed"
)
args <- parser$parse_args()
n <- args[["n_subjects"]]
seed <- args[["seed"]]

cat(
  ">>>>>>>>>>> Experiment start <<<<<<<<<<< \n",
  "n =", n, "\n",
  "seed =", seed, "\n"
)

genData1D <- function(n, K = 5) {

  t0 <- 0
  t1 <- 1
  nt <- 100 + 1
  tspan <- t1 - t0
  tgrid <- seq(t0, t1, length.out = nt)

  # Fourier basis as FPCs
  Phi <- fourier(tgrid, nbasis=K+1)[,2:(K+1)] / sqrt(tspan)
  lam <- 1 * (1:K)^(-2)
  covMat <- crossprod(t(sweep(Phi, 2, sqrt(lam), "*")))

  # functional data
  scores <- matrix(rnorm(n * K), nrow = n, ncol = K)
  scores <- sweep(scores, 2, sqrt(lam), "*")
  X <- scores %*% t(Phi)
  
  list(
    tgrid = tgrid,
    covMat = covMat,
    Phi = Phi,
    lam = lam,
    scores = scores,
    X = X
  )
}

sampleFuncData1D <- function(X, tgrid, m) {
  n <- nrow(X)
  nt <- ncol(X)
  stopifnot(nt == length(tgrid))

  Ltid <- lapply(seq_len(n), \(i) sort(sample(1:nt, m)))
  Lt <- lapply(Ltid, \(tid) tgrid[tid])
  Ly <- lapply(
    seq_len(n), \(i) {
      tid <- Ltid[[i]]
      X[i, tid]
    }
  )
  
  list(
    Ly = Ly,
    Ltid = Ltid,
    Lt = Lt
  )
}

getRawCovs <- function(Ly, Ltid, tgrid) {

  n <- length(Ly)
  raw_covs <- foreach(
    i = seq_len(n),
    .combine = rbind
  ) %do%
    {
      # get the working K
      Ci <- outer(Ly[[i]], Ly[[i]])

      # only keep the lower half due to symmetry
      low_tri_id <- lower.tri(Ci, diag = TRUE)
      Ci <- Ci[low_tri_id]

      # create a meshgrid of (t1, t2)
      t_cov_mesh <- pracma::meshgrid(tgrid[Ltid[[i]]])
      tid_mesh <- pracma::meshgrid(Ltid[[i]])

      # only keep the corresponding lower half (t2 > t1)
      t_cov_grid <- cbind(t_cov_mesh$X[low_tri_id], t_cov_mesh$Y[low_tri_id])
      tid_grid <- cbind(tid_mesh$X[low_tri_id], tid_mesh$Y[low_tri_id])

      # points weighting
      w <- rep(1, length(Ci))

      data.table::data.table(
        t1 = t_cov_grid[, 1],
        t2 = t_cov_grid[, 2],
        t1id = tid_grid[, 1],
        t2id = tid_grid[, 2],
        Cval = as.vector(Ci),
        w = w
      )
    }
  
  raw_covs
}

scaleCovVar <- function(rawcov, covMat, scaling = 1) {
  t1id <- rawcov$t1id
  t2id <- rawcov$t2id
  Cvals <- rawcov$Cval
  Ctrue <- covMat[cbind(t1id, t2id)]
  Cvals <- Ctrue + (Cvals - Ctrue) * scaling
  rawcov$Cval <- Cvals
  rawcov
}

estCov <- function(rawcov, tgrid, nb=7) {

  nt <- length(tgrid)
  rawcov$w <- rawcov$w / mean(rawcov$w)
  fit_C <- gam(
    Cval ~ te(t1, t2, k = c(nb, nb), bs = "ps"),
    data = rawcov,
    weights = rawcov$w
  )

  Cmat <- matrix(0, nt, nt)
  low_tri_id <- lower.tri(Cmat, diag = TRUE)
  t_cov_mesh <- pracma::meshgrid(tgrid)
  t_cov_grid <- cbind(t_cov_mesh$X[low_tri_id], t_cov_mesh$Y[low_tri_id])
  Cmat[low_tri_id] <- predict(
    fit_C,
    newdata = data.table::data.table(
      t1 = t_cov_grid[, 1],
      t2 = t_cov_grid[, 2]
    )
  )
  C <- Cmat + t(Cmat)
  diag(C) <- diag(Cmat)
  C
}


# Experiment -------------------------------------------------------------

set.seed(seed)
fulldat <- genData1D(n)

mList <- seq(4,20,4)
scalingList <- c(3,1,0.3,0.1,0.03,0.01,0.003,0.001)
nbList <- seq(5,17,4)

for (m in mList) {

  cat("----- m =", m, "------\n")

  set.seed(seed)
  dat <- sampleFuncData1D(fulldat$X, fulldat$tgrid, m)

  dir_path <- file.path("test", "results")
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
  file_name <- paste0("testcov_n",n,"_m",m,".csv")

  res = tibble(
    seed = numeric(),
    n = numeric(),
    m = numeric(),
    nb = numeric(),
    s = numeric(),
    mse = numeric()
  )

  rawcov <- getRawCovs(dat$Ly, dat$Ltid, fulldat$tgrid)

  for (scaling in scalingList) {

    for (nb in nbList) {

      Chat <- estCov(
        scaleCovVar(rawcov, fulldat$covMat, scaling),
        fulldat$tgrid, nb = nb
      )
      mse <- mean((Chat - fulldat$covMat)^2)

      res <- add_row(
        res,
        seed = seed, n = n, m = m, nb = nb, s = scaling, mse = mse
      )

    }

  }

  readr::write_csv(res, file = file.path(dir_path, file_name), append = TRUE)

}
