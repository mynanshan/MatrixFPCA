#!/usr/bin/env Rscript

parser <- argparse::ArgumentParser(description = "Summarize method-comparison CSV files.")
parser$add_argument("--results_dir", default = file.path("experiment", "results-methods"))
parser$add_argument("--output", default = file.path("experiment", "results-methods", "summary.csv"))
parser$add_argument("--manifest_dir", default = file.path("logs", "methods", "manifests"))
parser$add_argument("--status_output",
                    default = file.path("experiment", "results-methods", "scheduler-status.csv"))
args <- parser$parse_args()
if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required.")

files <- list.files(args$results_dir, pattern = "^methods_.*\\.csv$", full.names = TRUE)
if (!length(files)) stop("No method result CSV files found in ", args$results_dir, ".")
results <- data.table::rbindlist(lapply(files, data.table::fread), fill = TRUE)
group <- c(
  "method", "n", "p", "q", "nt", "sampling_min", "sampling_max", "sig",
  "rank_row", "rank_col", "rank_time", "eigensystem", "score_distribution",
  "score_dependence", "dependence_param", "eigvec_set", "eigval_set", "eval_ncomp"
)
summary <- results[, .(
  jobs = .N,
  successful = sum(status == "ok"),
  success_rate = mean(status == "ok"),
  mean_rel_hs_cov = mean(rel_hs_cov, na.rm = TRUE),
  sd_rel_hs_cov = stats::sd(rel_hs_cov, na.rm = TRUE),
  mean_oracle_rank_rel_hs = mean(oracle_rank_rel_hs, na.rm = TRUE),
  mean_excess_rel_hs_cov = mean(excess_rel_hs_cov, na.rm = TRUE),
  mean_sin2 = mean(mean_sin2, na.rm = TRUE),
  mean_runtime = mean(runtime, na.rm = TRUE),
  median_runtime = stats::median(runtime, na.rm = TRUE),
  max_runtime = max(runtime, na.rm = TRUE)
), by = group]
numeric_columns <- names(summary)[vapply(summary, is.numeric, logical(1L))]
for (column in numeric_columns) {
  data.table::set(summary, which(!is.finite(summary[[column]])), column, NA_real_)
}
data.table::setorderv(summary, c("p", "n", "method"))
dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(summary, args$output)
print(summary)
cat("Wrote ", args$output, "\n", sep = "")

manifest_files <- list.files(
  args$manifest_dir, pattern = "^simu_methods_.*\\.tsv$", full.names = TRUE
)
if (length(manifest_files)) {
  manifests <- data.table::rbindlist(lapply(manifest_files, data.table::fread), fill = TRUE)
  manifests <- unique(manifests, by = "job_id", fromLast = TRUE)
  status_dir <- file.path(args$results_dir, "status")
  manifests[, scheduler_status := vapply(job_id, function(id) {
    path <- file.path(status_dir, paste0(id, ".status"))
    if (!file.exists(path)) return("pending")
    strsplit(readLines(path, n = 1L, warn = FALSE), "\\t")[[1L]][1L]
  }, character(1L))]
  status_summary <- manifests[, .(
    tasks = .N,
    ok = sum(scheduler_status == "ok"),
    timeout = sum(scheduler_status == "timeout"),
    error = sum(scheduler_status == "error"),
    pending = sum(scheduler_status == "pending")
  ), by = .(scenario, method_group, n, p, q, nt, s1, s2, sigma,
            distribution, dependence, dep_param, eigensystem)]
  data.table::setorderv(status_summary, c("p", "n", "scenario", "method_group"))
  dir.create(dirname(args$status_output), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(status_summary, args$status_output)
  cat("Wrote ", args$status_output, "\n", sep = "")
}
