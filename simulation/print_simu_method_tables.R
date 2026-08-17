#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

parser <- argparse::ArgumentParser(
  description = paste(
    "Print method-comparison tables from the numeric summary produced by",
    "analyze_simu_methods.R."
  )
)
parser$add_argument(
  "--input",
  default = file.path(
    "simulation", "summary-methods", "method_comparison_complete.csv"
  ),
  help = "Numeric comparison CSV from analyze_simu_methods.R (default: %(default)s)."
)
parser$add_argument(
  "--metrics",
  default = paste(
    c("rel_hs_cov", "excess_rel_hs_cov", "sin2", "runtime"),
    collapse = ","
  ),
  help = paste(
    "Comma-separated metrics to print. Available: rel_hs_cov,",
    "oracle_rank_rel_hs, excess_rel_hs_cov, sin2, eigenvalue_sse, runtime,",
    "or all (default: %(default)s)."
  )
)
parser$add_argument(
  "--scenarios", default = NULL,
  help = paste(
    "Optional comma-separated scenario IDs or labels. By default all",
    "scenarios are printed in analyzer order."
  )
)
parser$add_argument(
  "--digits", type = "integer", default = 4L,
  help = "Digits after the decimal for accuracy metrics (default: %(default)s)."
)
parser$add_argument(
  "--runtime_digits", type = "integer", default = 2L,
  help = "Digits after the decimal for runtime (default: %(default)s)."
)
args <- parser$parse_args()

if (!file.exists(args$input)) {
  stop(
    "Comparison summary does not exist: ", args$input, "\n",
    "Run Rscript simulation/analyze_simu_methods.R first."
  )
}
if (args$digits < 0L || args$digits > 10L) {
  stop("--digits must be between 0 and 10.")
}
if (args$runtime_digits < 0L || args$runtime_digits > 10L) {
  stop("--runtime_digits must be between 0 and 10.")
}

metric_definitions <- data.table(
  metric = c(
    "rel_hs_cov", "oracle_rank_rel_hs", "excess_rel_hs_cov", "sin2",
    "eigenvalue_sse", "runtime"
  ),
  mean_column = c(
    "mean_rel_hs_cov", "mean_oracle_rank_rel_hs",
    "mean_excess_rel_hs_cov", "mean_sin2", "mean_rel_eigenvalue_sse",
    "mean_runtime"
  ),
  sd_column = c(
    "sd_rel_hs_cov", "sd_oracle_rank_rel_hs", "sd_excess_rel_hs_cov",
    "sd_mean_sin2", "sd_rel_eigenvalue_sse", "sd_runtime"
  ),
  label = c(
    "Relative HS error", "Oracle-rank error", "Excess relative HS",
    "Mean sin^2 angle", "Relative eigenvalue SSE", "Runtime (s)"
  )
)

split_argument <- function(value) {
  if (is.null(value) || !nzchar(trimws(value))) return(character())
  unique(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
}

requested_metrics <- split_argument(args$metrics)
if (identical(requested_metrics, "all")) {
  requested_metrics <- metric_definitions$metric
}
unknown_metrics <- setdiff(requested_metrics, metric_definitions$metric)
if (length(unknown_metrics)) {
  stop("Unknown metric(s): ", paste(unknown_metrics, collapse = ", "))
}
if (!length(requested_metrics)) stop("At least one metric must be requested.")
metric_definitions <- metric_definitions[match(requested_metrics, metric)]

comparison <- fread(args$input)
options(width = max(getOption("width"), 200L))
identifier_columns <- c(
  "scenario_order", "scenario", "scenario_label", "method_order", "method",
  "method_label", "expected", "successful"
)
required_columns <- c(
  identifier_columns,
  metric_definitions$mean_column,
  metric_definitions$sd_column
)
missing_columns <- setdiff(required_columns, names(comparison))
if (length(missing_columns)) {
  stop(
    "Comparison summary lacks required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

requested_scenarios <- split_argument(args$scenarios)
if (length(requested_scenarios)) {
  known_scenarios <- unique(comparison[, c("scenario", "scenario_label")])
  unknown_scenarios <- setdiff(
    requested_scenarios,
    c(known_scenarios$scenario, known_scenarios$scenario_label)
  )
  if (length(unknown_scenarios)) {
    stop("Unknown scenario(s): ", paste(unknown_scenarios, collapse = ", "))
  }
  comparison <- comparison[
    scenario %in% requested_scenarios | scenario_label %in% requested_scenarios
  ]
}
setorder(comparison, scenario_order, method_order)

format_mean_sd <- function(mean_value, sd_value, digits) {
  out <- rep("--", length(mean_value))
  has_mean <- is.finite(mean_value)
  has_sd <- has_mean & is.finite(sd_value)
  number_format <- paste0("%.", digits, "f")
  mean_sd_format <- paste0(number_format, " (", number_format, ")")
  out[has_mean] <- sprintf(number_format, mean_value[has_mean])
  out[has_sd] <- sprintf(
    mean_sd_format, mean_value[has_sd], sd_value[has_sd]
  )
  out
}

overview_path <- file.path(dirname(args$input), "completion_overview.csv")
if (file.exists(overview_path)) {
  overview <- fread(overview_path)
  overview_columns <- intersect(
    c(
      "profile", "study_version", "expected_repetitions",
      "successful_method_fits", "expected_method_fits", "failed_method_fits",
      "missing_method_fits", "complete"
    ),
    names(overview)
  )
  cat("Completion overview\n")
  print(overview[, ..overview_columns], row.names = FALSE)
}

for (scenario_id in unique(comparison$scenario)) {
  section <- comparison[scenario == scenario_id]
  output <- section[, .(
    Method = method_label,
    Successful = paste0(successful, "/", expected)
  )]

  for (metric_index in seq_len(nrow(metric_definitions))) {
    definition <- metric_definitions[metric_index]
    digits <- if (definition$metric == "runtime") {
      args$runtime_digits
    } else {
      args$digits
    }
    output[, (definition$label) := format_mean_sd(
      section[[definition$mean_column]],
      section[[definition$sd_column]],
      digits
    )]
  }

  setcolorder(output, c("Method", metric_definitions$label, "Successful"))
  cat("\n", section$scenario_label[[1L]], " [", scenario_id, "]\n", sep = "")
  print(output, row.names = FALSE)
}

cat(
  "\nEntries are Monte Carlo mean (SD) over successful fits; -- denotes an ",
  "unavailable metric.\n",
  sep = ""
)
