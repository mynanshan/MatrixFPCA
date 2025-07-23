#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/r/4.5.0/bin/Rscript
# get_mem_per_core.R

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: get_mem_per_core.R <p>")
}
p <- as.numeric(args[1])

# ==== YOUR MEM‐CHOICE LOGIC HERE ====
mem <- min(1, max(8, 2^ceiling(log10(p))))

cat(num_cores)
