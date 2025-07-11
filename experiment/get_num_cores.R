#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/r/4.5.0/bin/Rscript
# get_num_cores.R

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: get_cores.R <p>")
}
p <- as.numeric(args[1])

# ==== YOUR CORE‐CHOICE LOGIC HERE ====
# e.g. 1 core per 200 rows, minimum 1, maximum 32
num_cores <- min(32, max(1, 2^ceiling(log10(p)) + 1))

cat(num_cores)
