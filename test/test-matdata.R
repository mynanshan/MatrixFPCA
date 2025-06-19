library(rTensor)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RhpcBLASctl)

blas_set_num_threads(4)

res = tibble(
  n = numeric(),
  p = numeric(),
  q = numeric(),
  distr = character(),
  mseU = numeric(),
  sumkVarXX = numeric(),
  avgkVarXX = numeric(),
  maxkVarXX = numeric(),
  sumkCovXX = numeric(),
  avgkCovXX = numeric(),
  avgabskCovXX = numeric(),
  maxkCovXX = numeric()
)

for (irep in 1:20) {

for (p in c(10, 20, 40, 80, 120)) {
  
  n = 200
  q = p + 1
  sig = 0.05
  
  
  Phi = pracma::randortho(p)
  Psi = pracma::randortho(q)
  
  lamU = rev(seq(0.1, 1, length.out=p))
  lamV = rev(seq(0.1, 1, length.out=q))
  lamU = lamU / sum(lamU) * p
  lamV = lamV / sum(lamV) * q
  
  U = Phi %*% diag(lamU,p,p) %*% t(Phi)
  V = Psi %*% diag(lamV,q,q) %*% t(Psi)
  
  Z = rnorm(n * p * q) |>
    array(dim = c(n, p, q)) |> 
    as.tensor()
  
  scores = ttl(
    Z, list(
      diag(sqrt(lamU), p, p),
      diag(sqrt(lamV), q, q) 
    ), c(2:3)
  )
  
  X = ttl(scores, list(Phi, Psi), c(2:3))
  X = X@data
  
  for (distr in c("gaussian", "unif", "student")) {
    
    if (distr == "gaussian") {
      X = X + rnorm(n*p*q, 0, sqrt(sig))
    } else if (distr == "unif") {
      width = sqrt(12*sig)
      X = X + runif(n*p*q, -0.5*width, 0.5*width)
    } else if (distr == "student") {
      # nu = 4, var = nu / (nu - 2) = 2
      # scaling: sqrt(sig / 2)
      X = X + rt(n*p*q, 4) * sqrt(sig / 2)
    }
    
    Utildes = apply(X, 1, \(Xi) Xi %*% t(Xi)) |> 
      array(dim = c(p, p, n)) |> 
      aperm(c(3, 1, 2))  # (n, p, p)
    Uhat = apply(Utildes, c(2,3), mean) / q
    
    mseU = mean((Uhat - U)^2)
    print(mseU)
    
    covXX = c()
    jjid = c()
    for (j1 in seq(1, p, ceiling(p / 20))) {
      for (j2 in seq(j1, p, ceiling(p / 20))) {
        Xj1j2 = X[,j1,] * X[,j2,]  # (n, q)
        covXX = rbind(covXX, c(cov(Xj1j2)))
        if (j1 == j2) {
          jjid = c(jjid, nrow(covXX))
        }
      }
    }
    covXX = array(covXX, c(dim(covXX)[1], q, q))
    index = expand.grid(k1=1:q, k2=1:q)
    cat("Cov(Xjk, Xj'k) \n")
    print(summary(c(apply(covXX[jjid,,], 1, diag))))
    print(summary(c(apply(covXX[-jjid,,], 1, diag))))
    sumkVarXX = mean(apply(covXX, 1, \(C) sum(diag(C))))
    avgkVarXX = mean(apply(covXX, 1, \(C) mean(diag(C))))
    maxkVarXX = mean(apply(covXX, 1, \(C) max(diag(C))))
    cat("Cov(Xjk, Xj'k') \n")
    print(summary(c(apply(covXX[jjid,,], 1, \(C) C[index$k1 != index$k2]))))
    print(summary(c(apply(covXX[-jjid,,], 1, \(C) C[index$k1 != index$k2]))))
    sumkCovXX = mean(apply(covXX, 1, \(C) sum(C[index$k1 != index$k2])))
    avgkCovXX = mean(apply(covXX, 1, \(C) mean(C[index$k1 != index$k2])))
    avgabskCovXX = mean(apply(covXX, 1, \(C) mean(abs(C[index$k1 != index$k2]))))
    maxkCovXX = mean(apply(covXX, 1, \(C) max(abs(C[index$k1 != index$k2]))))
    
    res = add_row(
      res,
      n = n, p = p, q = q, distr = distr,
      mseU = mseU,
      sumkVarXX = sumkVarXX,
      avgkVarXX = avgkVarXX,
      maxkVarXX = maxkVarXX,
      sumkCovXX = sumkCovXX,
      avgkCovXX = avgkCovXX,
      avgabskCovXX = avgabskCovXX,
      maxkCovXX = maxkCovXX
    ) 
  }
}
  
}

readr::write_csv(res, "test/matdata-res.csv")

res = readr::read_csv("test/matdata-res.csv")

ggplot(res) +
  geom_point(aes(x = p, y = mseU * q)) + 
  facet_wrap(vars(distr))

ggplot(res) +
  geom_point(aes(x = p, y = sumkVarXX / q)) + 
  facet_wrap(vars(distr))

ggplot(res) +
  geom_point(aes(x = p, y = avgkVarXX)) + 
  facet_wrap(vars(distr))

ggplot(res) +
  geom_point(aes(x = p, y = maxkVarXX)) + 
  facet_wrap(vars(distr))

ggplot(res) +
  geom_point(aes(x = p, y = sumkCovXX / q)) + 
  facet_wrap(vars(distr))

ggplot(res) +
  geom_point(aes(x = p, y = avgabskCovXX)) + 
  facet_wrap(vars(distr))

ggplot(res) +
  geom_point(aes(x = p, y = maxkCovXX / log(q))) + 
  facet_wrap(vars(distr))

