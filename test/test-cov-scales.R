library(fda)
library(Matrix)
library(rTensor)
library(data.table)

# dev.off()
pa = par(no.readonly = TRUE)

source(file.path("R", "matfdata.R"))
source(file.path("R", "matfpca.R"))
source(file.path("simulation", "gendata-simple.R"))

for (p in c(10, 20, 50, 100, 150, 200)) {
  
  # dimension and temporal space
  n = 100 # n sample
  q = p  # dim column
  nt = 50 + 1 # num time points
  t_range = c(0,1)
  
  npc_row = p
  npc_col = q
  npc_time = 3
  
  sampling_ratio = c(0.05, 0.1)  # (Pmin, Pmax)
  
  genData <- gen_data_simple(
    n, p, q,
    npc_row, npc_col, npc_time,
    nt, t_range,
    sampling_ratio
  )
  
  data <- genData$data
  cov_true <- genData$cov
  
  index = expand.grid(j1 = 1:p, j2 = 1:p)
  
  cat("p =", p, "\n")
  cat("Diagonal elements:\n")
  print(summary(c(diag(cov_true$U))))
  cat("Off-diagonal elements:\n")
  print(summary(c(cov_true$U[index$j1 != index$j2])))
  
  UdiagVar = c()
  UoffdVar = c()
  for (j1 in seq(1, p, ceiling(p/20))) {
    for (j2 in seq(j1, q, ceiling(p/20))) {
      UjjtildeList = sapply(
          1:n, \(i) {
            get_working_U(data$X_list[[i]])[j1,j2]
          }
        )
      Ujjtilde = mean(UjjtildeList)
      ENi = mean(sampling_ratio) * nt
      EK = mean(diag(cov_true$K))
      Ujjtrue = cov_true$U[j1,j2]
      UjjM2Hat = mean((UjjtildeList/(q * EK * ENi))^2)
      Ujj2 = Ujjtrue^2
      UjjVarHat = (UjjM2Hat - Ujj2)
      if (j1 == j2) {
        UdiagVar = c(UdiagVar, UjjVarHat)
      } else {
        UoffdVar = c(UoffdVar, UjjVarHat)
      }
      
      Xj1 = lapply(data$X_list, \(X) X[j1,,])
      Xj2 = lapply(data$X_list, \(X) X[j2,,])
      
    }
  }
  cat("Estimated variance of Ujj:\n")
  print(summary(UdiagVar))
  cat("Estimated variance of Ujj':\n")
  print(summary(UoffdVar))
  cat("q * Estimated variance of Ujj':\n")
  print(summary(UoffdVar * q))
  cat("\n")
}
