#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

parser <- argparse::ArgumentParser(
  description = paste(
    "Audit and summarize method-comparison simulations, write complete",
    "accuracy/runtime/completion tables, and print publication-style tables."
  )
)
parser$add_argument(
  "--results_dir", default = file.path("simulation", "results-methods"),
  help = "Directory containing methods_*.csv and status/ (default: %(default)s)."
)
parser$add_argument(
  "--output_dir", default = file.path("simulation", "summary-methods"),
  help = "Directory for generated audit and comparison tables (default: %(default)s)."
)
parser$add_argument(
  "--profile", choices = c("main", "scaling", "pilot"), default = "main",
  help = "Launcher scenario bundle to analyze (default: %(default)s)."
)
parser$add_argument(
  "--study_version", default = "v3",
  help = "Version tag emitted by run_simu_methods.sh (default: %(default)s)."
)
parser$add_argument(
  "--expected_repetitions", type = "integer",
  help = "Expected seeds starting at one (default: 100 for main/scaling; 3 for pilot)."
)
parser$add_argument(
  "--digits", type = "integer", default = 4L,
  help = "Digits after the decimal in formatted accuracy tables (default: %(default)s)."
)
parser$add_argument(
  "--runtime_digits", type = "integer", default = 2L,
  help = "Digits after the decimal for formatted runtimes (default: %(default)s)."
)
parser$add_argument(
  "--metrics",
  default = paste(c("rel_hs_cov", "sin2", "runtime"),
                  collapse = ","),
  help = paste(
    "Comma-separated metrics to print: rel_hs_cov, oracle_rank_rel_hs,",
    "excess_rel_hs_cov, sin2, eigenvalue_sse, runtime, or all",
    "(default: %(default)s)."
  )
)
parser$add_argument(
  "--scenarios", default = NULL,
  help = "Optional comma-separated scenario IDs or labels to print."
)
parser$add_argument(
  "--quiet", action = "store_true",
  help = "Write all output files but do not print scenario comparison tables."
)
parser$add_argument(
  "--fail_on_incomplete", action = "store_true",
  help = "Exit with status 1 after writing outputs if any expected fit is missing or failed."
)
args <- parser$parse_args()

if (!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", args$study_version)) {
  stop("--study_version must contain only letters, digits, underscores, or hyphens.")
}
if (!dir.exists(args$results_dir)) {
  stop("Results directory does not exist: ", args$results_dir)
}
if (args$digits < 0L || args$digits > 10L) {
  stop("--digits must be between 0 and 10.")
}
if (args$runtime_digits < 0L || args$runtime_digits > 10L) {
  stop("--runtime_digits must be between 0 and 10.")
}
expected_repetitions <- if (is.null(args$expected_repetitions)) {
  if (args$profile == "pilot") 3L else 100L
} else {
  args$expected_repetitions
}
if (expected_repetitions < 1L) stop("--expected_repetitions must be positive.")

result_files <- list.files(
  args$results_dir, pattern = "^methods_.*[.]csv$", full.names = TRUE
)
if (!length(result_files)) {
  stop("No methods_*.csv files found in ", args$results_dir, ".")
}

read_result <- function(path) {
  out <- tryCatch(
    fread(path),
    error = function(error) {
      warning("Could not read ", path, ": ", conditionMessage(error))
      NULL
    }
  )
  if (!is.null(out)) out[, source_file := basename(path)]
  out
}

results <- rbindlist(lapply(result_files, read_result), use.names = TRUE, fill = TRUE)
required_columns <- c(
  "method", "status", "error", "seed", "n", "p", "q", "nt",
  "sampling_min", "sampling_max", "sig", "rank_row", "rank_col",
  "rank_time", "eigensystem", "eigvec_set", "eigval_set",
  "score_distribution", "score_dependence", "dependence_param",
  "result_tag", "rel_hs_cov", "oracle_rank_rel_hs",
  "excess_rel_hs_cov", "mean_sin2", "rel_eigenvalue_sse", "runtime",
  "cov_ncomp", "happ_npc", "product_spatial_npc", "product_temporal_npc",
  "mgs_rank_candidates", "mgs_basis_candidates", "shi_nbasis_time",
  "shi_nbasis_row", "shi_nbasis_col", "shi_smoothness", "shi_score_ridge",
  "spaco_smoothness", "spaco_ridge", "spaco_score_ridge",
  "spaco_selected_score_ridge", "spaco_minimum_cv_score_ridge",
  "spaco_score_ridge_candidates", "spaco_score_ridge_cv_error",
  "spaco_score_ridge_cv_se", "spaco_maxit",
  "observed_visits_min", "observed_visits_mean", "observed_visits_max"
)
missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns)) {
  stop("Result files lack required columns: ", paste(missing_columns, collapse = ", "))
}

# The result tag is the launcher's stable task identifier:
# <study-version>_<scenario>_sd<seed>_<method-group>. Versioning prevents
# --no-overwrite from mistaking results from an older tuning configuration for
# completed current tasks. Legacy unversioned tags remain parseable for audits.
all_groups <- c(
  "proposal", "spline", "product", "spaco", "happ_pace", "happ_face",
  "mfaces", "mgsfpca", "shi"
)
group_pattern <- paste(all_groups, collapse = "|")
legacy_tag_pattern <- paste0("^(.*)_sd([0-9]+)_(", group_pattern, ")$")
requested_tag_pattern <- paste0(
  "^", args$study_version, "_(.*)_sd([0-9]+)_(", group_pattern, ")$"
)
requested_tag <- grepl(requested_tag_pattern, results$result_tag)
legacy_tag <- !requested_tag & grepl(legacy_tag_pattern, results$result_tag)
valid_tag <- requested_tag | legacy_tag
results[, `:=`(
  parsed_study_version = fifelse(
    requested_tag, args$study_version, fifelse(legacy_tag, "legacy", NA_character_)
  ),
  parsed_scenario = fifelse(
    requested_tag,
    sub(requested_tag_pattern, "\\1", result_tag),
    fifelse(legacy_tag, sub(legacy_tag_pattern, "\\1", result_tag), NA_character_)
  ),
  tag_seed = suppressWarnings(as.integer(ifelse(
    requested_tag, sub(requested_tag_pattern, "\\2", result_tag),
    ifelse(legacy_tag, sub(legacy_tag_pattern, "\\2", result_tag), NA_character_)
  ))),
  method_group = fifelse(
    requested_tag,
    sub(requested_tag_pattern, "\\3", result_tag),
    fifelse(legacy_tag, sub(legacy_tag_pattern, "\\3", result_tag), NA_character_)
  )
)]

scenario_definitions <- rbindlist(list(
  data.table(
    scenario = c("dim_p8", "dim_p16", "dim_p32"),
    scenario_label = c(
      "Baseline dimension 8 x 8", "Baseline dimension 16 x 16",
      "Baseline dimension 32 x 32"
    ),
    n = 50L, p = c(8L, 16L, 32L), q = c(8L, 16L, 32L),
    nt = 51L, sampling_min = 0.1, sampling_max = 0.2, sig = 0.1,
    score_distribution = "normal", score_dependence = "independent",
    dependence_param = 0, eigensystem = "ar1", eigvec_set = 2L,
    eigval_set = 1L, rank_time = 5L
  ),
  data.table(
    scenario = c(
      "n25_p8", "n100_p8", "visits_20_40_p8", "noise_05_p8",
      "spectral_p8", "t5_p8", "dependent_p8"
    ),
    scenario_label = c(
      "Sample size n = 25", "Sample size n = 100",
      "11-20 visits (20-40% grid)", "Noise sigma = 0.5",
      "Alternative spectrum", "t(5) scores", "Dependent scores"
    ),
    n = c(25L, 100L, rep(50L, 5L)), p = 8L, q = 8L, nt = 51L,
    sampling_min = c(0.1, 0.1, 0.2, rep(0.1, 4L)),
    sampling_max = c(0.2, 0.2, 0.4, rep(0.2, 4L)),
    sig = c(rep(0.1, 3L), 0.5, 0.1, 0.1, 0.1),
    score_distribution = c(rep("normal", 5L), "t5", "normal"),
    score_dependence = c(rep("independent", 6L), "common_amplitude"),
    dependence_param = c(rep(0, 6L), 1),
    eigensystem = c(rep("ar1", 4L), "spectral", "ar1", "ar1"),
    eigvec_set = 2L, eigval_set = c(rep(1L, 4L), 2L, 1L, 1L),
    rank_time = c(rep(5L, 4L), 9L, 5L, 5L)
  ),
  data.table(
    scenario = "baseline_p8",
    scenario_label = "Pilot baseline 8 x 8",
    n = 50L, p = 8L, q = 8L, nt = 51L,
    sampling_min = 0.1, sampling_max = 0.2, sig = 0.1,
    score_distribution = "normal", score_dependence = "independent",
    dependence_param = 0, eigensystem = "ar1", eigvec_set = 2L,
    eigval_set = 1L, rank_time = 5L
  ),
  data.table(
    scenario = c("rank6_p8", "rank12_p8"),
    scenario_label = c(
      "Covariance rank K = 6 (8 x 8 baseline)",
      "Covariance rank K = 12 (8 x 8 baseline)"
    ),
    n = 50L, p = 8L, q = 8L, nt = 51L,
    sampling_min = 0.1, sampling_max = 0.2, sig = 0.1,
    score_distribution = "normal", score_dependence = "independent",
    dependence_param = 0, eigensystem = "ar1", eigvec_set = 2L,
    eigval_set = 1L, rank_time = 5L
  )
), use.names = TRUE)
scenario_definitions[, `:=`(
  cov_ncomp = 9L,
  happ_npc = 3L,
  product_spatial_npc = 10L,
  product_temporal_npc = 5L,
  mgs_rank_candidates = "3",
  mgs_basis_candidates = "7",
  shi_nbasis_time = 7L,
  shi_nbasis_row = 5L,
  shi_nbasis_col = 5L,
  shi_smoothness = 1e-2,
  shi_score_ridge = "auto",
  spaco_smoothness = 1e-2,
  spaco_ridge = 1e-8,
  spaco_score_ridge = "auto",
  spaco_maxit = 100L
)]
scenario_definitions[scenario == "rank6_p8", cov_ncomp := 6L]
scenario_definitions[scenario == "rank12_p8", cov_ncomp := 12L]

profile_scenarios <- switch(
  args$profile,
  main = c(
    "dim_p8", "rank6_p8", "rank12_p8", "dim_p16", "dim_p32",
    "n25_p8", "n100_p8",
    "visits_20_40_p8", "noise_05_p8", "spectral_p8", "t5_p8",
    "dependent_p8"
  ),
  scaling = c("dim_p8", "dim_p16", "dim_p32"),
  pilot = "baseline_p8"
)
scenarios <- scenario_definitions[match(profile_scenarios, scenario)]
scenarios[, scenario_order := seq_len(.N)]

method_map <- data.table(
  method_group = all_groups,
  group_order = seq_along(all_groups),
  methods = c(
    "matfpca_refined,matfpca_identity", "spline_pca", "sparse_product_fpca",
    "spaco_minus", "happ_pace", "happ_face", "mfaces", "mgsfpca",
    "shi_mdfpca"
  )
)
method_labels <- data.table(
  method = c(
    "matfpca_refined", "matfpca_identity", "sparse_product_fpca",
    "spline_pca", "spaco_minus", "happ_pace", "happ_face", "mfaces",
    "mgsfpca", "shi_mdfpca"
  ),
  method_label = c(
    "MatFPCA refined", "MatFPCA identity", "Sparse Product FPCA",
    "Spline PCA", "SPACO-minus", "HAPP-PACE", "HAPP-FACE", "mFACEs",
    "mGSFPCA", "Shi MD-FPCA"
  ),
  method_order = seq_len(10L)
)

# Keep these eligibility rules synchronized with run_simu_methods.sh. They are
# encoded here instead of inferred from observed files so an entirely missing
# method group is still detected by the audit.
max_dimension <- c(
  proposal = 32L, spline = 32L, product = 32L, spaco = 32L,
  happ_pace = 32L, happ_face = 32L, mfaces = 8L, mgsfpca = 32L, shi = 32L
)
scenario_allowed <- function(scenario, group) {
  TRUE
}

expected_tasks <- scenarios[
  , .(seed = seq_len(expected_repetitions)),
  by = names(scenarios)
][
  , .(method_group = all_groups),
  by = c(names(scenarios), "seed")
]
expected_tasks <- expected_tasks[
  mapply(scenario_allowed, scenario, method_group) &
    p <= unname(max_dimension[method_group]) &
    q <= unname(max_dimension[method_group])
]
expected_tasks[, job_id := paste0(
  args$study_version, "_", scenario, "_sd", seed, "_", method_group
)]
expected_tasks <- merge(expected_tasks, method_map, by = "method_group", sort = FALSE)
setorder(expected_tasks, scenario_order, seed, group_order)

expected_methods <- expected_tasks[
  , .(method = unlist(strsplit(methods, ",", fixed = TRUE))),
  by = setdiff(names(expected_tasks), "methods")
]
expected_methods <- merge(expected_methods, method_labels, by = "method", sort = FALSE)
setorder(expected_methods, scenario_order, method_order, seed)

actual_tasks <- results[
  , .(
    result_files = uniqueN(source_file),
    result_rows = .N,
    observed_methods = uniqueN(method)
  ),
  by = .(result_tag)
]
task_audit <- merge(
  expected_tasks[, !"methods"], actual_tasks,
  by.x = "job_id", by.y = "result_tag", all.x = TRUE, sort = FALSE
)
task_audit[, result_file_present := !is.na(result_files)]
task_audit[is.na(result_files), `:=`(
  result_files = 0L, result_rows = 0L, observed_methods = 0L
)]
status_dir <- file.path(args$results_dir, "status")
task_audit[, scheduler_status := vapply(job_id, function(id) {
  path <- file.path(status_dir, paste0(id, ".status"))
  if (!file.exists(path)) return("missing")
  line <- readLines(path, n = 1L, warn = FALSE)
  if (!length(line)) return("empty")
  strsplit(line, "\t", fixed = TRUE)[[1L]][1L]
}, character(1L))]
task_audit[, scheduler_ok := scheduler_status == "ok"]
setorder(task_audit, scenario_order, seed, group_order)

actual_method_rows <- results[
  , .(
    rows_found = .N,
    observed_status = paste(sort(unique(status)), collapse = ";"),
    error = paste(sort(unique(error[!is.na(error) & nzchar(error)])), collapse = " | "),
    source_files = paste(sort(unique(source_file)), collapse = ";")
  ),
  by = .(result_tag, method)
]
method_audit <- merge(
  expected_methods, actual_method_rows,
  by.x = c("job_id", "method"), by.y = c("result_tag", "method"),
  all.x = TRUE, sort = FALSE
)
method_audit[is.na(rows_found), rows_found := 0L]
method_audit[is.na(observed_status), observed_status := "missing"]
method_audit[is.na(error), error := ""]
method_audit[, method_complete := rows_found == 1L & observed_status == "ok"]
setorder(method_audit, scenario_order, method_order, seed)

# Verify that a correctly named file also contains the settings promised by
# its task tag. This catches copied, renamed, or partially overwritten results
# that a filename/status-only completion check would miss.
config_fields <- c(
  "seed", "n", "p", "q", "nt", "sampling_min", "sampling_max", "sig",
  "rank_row", "rank_col", "rank_time", "eigensystem", "eigvec_set", "eigval_set",
  "score_distribution", "score_dependence", "dependence_param",
  "cov_ncomp", "happ_npc", "product_spatial_npc", "product_temporal_npc",
  "mgs_rank_candidates", "mgs_basis_candidates", "shi_nbasis_time",
  "shi_nbasis_row", "shi_nbasis_col", "shi_smoothness", "shi_score_ridge",
  "spaco_smoothness", "spaco_ridge", "spaco_score_ridge", "spaco_maxit"
)
actual_config <- results[
  result_tag %in% expected_methods$job_id,
  c("result_tag", "method", "source_file", config_fields), with = FALSE
]
setnames(actual_config, config_fields, paste0("actual_", config_fields))
expected_config <- copy(expected_methods)
expected_config[, `:=`(rank_row = p, rank_col = q)]
expected_config <- expected_config[
  , c("job_id", "method", config_fields), with = FALSE
]
setnames(expected_config, "job_id", "result_tag")
setnames(expected_config, config_fields, paste0("expected_", config_fields))
configuration_audit <- merge(
  expected_config, actual_config,
  by = c("result_tag", "method"), all.x = TRUE, sort = FALSE
)
different_value <- function(actual, expected) {
  missing_difference <- xor(is.na(actual), is.na(expected))
  both_present <- !is.na(actual) & !is.na(expected)
  value_difference <- if (is.numeric(actual) || is.integer(actual)) {
    both_present & abs(as.numeric(actual) - as.numeric(expected)) > 1e-10
  } else {
    both_present & as.character(actual) != as.character(expected)
  }
  missing_difference | value_difference
}
configuration_audit[, mismatch_fields := ""]
for (field in config_fields) {
  mismatch <- different_value(
    configuration_audit[[paste0("actual_", field)]],
    configuration_audit[[paste0("expected_", field)]]
  )
  configuration_audit[
    mismatch,
    mismatch_fields := paste0(
      mismatch_fields, fifelse(nzchar(mismatch_fields), ";", ""), field
    )
  ]
}
configuration_audit[, configuration_matches := nzchar(mismatch_fields) == FALSE]
configuration_mismatches <- configuration_audit[
  !is.na(source_file) & configuration_matches == FALSE
]

expected_tags <- expected_tasks$job_id
unexpected_result_rows <- results[!result_tag %in% expected_tags]
malformed_result_rows <- results[!valid_tag | tag_seed != seed]
duplicate_method_rows <- results[
  , .N, by = .(result_tag, method)
][N > 1L]
missing_tasks <- task_audit[!result_file_present | !scheduler_ok]
missing_method_rows <- method_audit[rows_found == 0L]
method_failures <- method_audit[rows_found > 0L & observed_status != "ok"]

selected_results <- merge(
  results,
  expected_methods[, .(result_tag = job_id, method, scenario_order, scenario_label,
                       method_label, method_order)],
  by = c("result_tag", "method"), all = FALSE, sort = FALSE
)

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
sd_or_na <- function(x) {
  if (sum(!is.na(x)) < 2L) NA_real_ else stats::sd(x, na.rm = TRUE)
}
median_or_na <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)

successful_results <- selected_results[status == "ok"]
metric_summary <- successful_results[
  , .(
    successful = .N,
    mean_rel_hs_cov = mean_or_na(rel_hs_cov),
    sd_rel_hs_cov = sd_or_na(rel_hs_cov),
    median_rel_hs_cov = median_or_na(rel_hs_cov),
    mean_oracle_rank_rel_hs = mean_or_na(oracle_rank_rel_hs),
    sd_oracle_rank_rel_hs = sd_or_na(oracle_rank_rel_hs),
    mean_excess_rel_hs_cov = mean_or_na(excess_rel_hs_cov),
    sd_excess_rel_hs_cov = sd_or_na(excess_rel_hs_cov),
    mean_sin2 = mean_or_na(mean_sin2),
    sd_mean_sin2 = sd_or_na(mean_sin2),
    mean_rel_eigenvalue_sse = mean_or_na(rel_eigenvalue_sse),
    sd_rel_eigenvalue_sse = sd_or_na(rel_eigenvalue_sse),
    mean_runtime = mean_or_na(runtime),
    sd_runtime = sd_or_na(runtime),
    median_runtime = median_or_na(runtime),
    max_runtime = if (all(is.na(runtime))) NA_real_ else max(runtime, na.rm = TRUE),
    mean_selected_spaco_score_ridge = mean_or_na(spaco_selected_score_ridge),
    median_selected_spaco_score_ridge = median_or_na(
      spaco_selected_score_ridge
    )
  ),
  by = .(scenario_order, parsed_scenario, scenario_label, method_order,
         method, method_label)
]
setnames(metric_summary, "parsed_scenario", "scenario")

completion_by_method <- method_audit[
  , .(
    expected = .N,
    observed = sum(rows_found > 0L),
    successful_audit = sum(method_complete),
    failed = sum(rows_found > 0L & observed_status != "ok"),
    missing = sum(rows_found == 0L)
  ),
  by = .(scenario_order, scenario, scenario_label, method_order, method, method_label)
]
comparison <- merge(
  completion_by_method, metric_summary,
  by = c("scenario_order", "scenario", "scenario_label", "method_order",
         "method", "method_label"),
  all.x = TRUE, sort = FALSE
)
comparison[, successful := successful_audit]
comparison[, successful_audit := NULL]
comparison[, success_rate := successful / expected]
comparison[, mc_se_rel_hs_cov := sd_rel_hs_cov / sqrt(successful)]
comparison[, rank_rel_hs_cov := frank(mean_rel_hs_cov, ties.method = "min", na.last = "keep"),
           by = scenario]
setorder(comparison, scenario_order, method_order)

format_mean_sd <- function(
  mean_value, sd_value, digits = args$digits, bold = rep(FALSE, length(mean_value))
) {
  out <- rep("--", length(mean_value))
  has_mean <- is.finite(mean_value)
  has_both <- has_mean & is.finite(sd_value)
  out[has_mean] <- sprintf(paste0("%.", digits, "f"), mean_value[has_mean])
  out[has_both] <- sprintf(
    paste0("%.", digits, "f (%.", digits, "f)"),
    mean_value[has_both], sd_value[has_both]
  )
  bold <- rep_len(bold, length(out)) & has_mean
  out[bold] <- paste0("\\textbf{", out[bold], "}")
  out
}

# Paper-facing results contain only the two principal accuracy criteria and
# fitting time. Audit counts remain available in the dedicated audit files,
# but are deliberately absent here. Smooth-then-PCA is a diagnostic weak
# baseline and is not part of the paper comparison.
paper_method_order <- data.table(
  method = c(
    "matfpca_refined", "matfpca_identity", "shi_mdfpca", "mfaces",
    "happ_face", "happ_pace", "mgsfpca", "spaco_minus",
    "sparse_product_fpca"
  ),
  method_short = c(
    "MatFPCA (refined)", "MatFPCA (identity)", "Shi MD-FPCA", "mFACEs",
    "HAPP-FACE", "HAPP-PACE", "mGSFPCA", "SPACO-minus", "Product FPCA"
  ),
  paper_method_order = seq_len(9L)
)
paper_scenario_order <- data.table(
  scenario = c(
    "dim_p8", "dim_p16", "dim_p32", "rank6_p8", "rank12_p8",
    "n25_p8", "n100_p8", "visits_20_40_p8", "noise_05_p8",
    "spectral_p8", "t5_p8", "dependent_p8"
  ),
  paper_scenario_label = c(
    "Baseline", "Dimension 16 x 16", "Dimension 32 x 32", "Rank K = 6",
    "Rank K = 12", "Sample size n = 25", "Sample size n = 100",
    "11--20 visits", "Noise sigma = 0.5", "Alternative spectrum",
    "t(5) scores", "Dependent scores"
  ),
  paper_scenario_latex = c(
    "Baseline", "Dimension $16 \\times 16$", "Dimension $32 \\times 32$",
    "Rank $K=6$", "Rank $K=12$", "Sample size $n=25$",
    "Sample size $n=100$", "11--20 visits", "Noise $\\sigma=0.5$",
    "Alternative spectrum", "$t_5$ scores", "Dependent scores"
  ),
  paper_scenario_order = seq_len(12L)
)

paper_results <- merge(
  comparison[method != "spline_pca"], paper_method_order,
  by = "method", all = FALSE, sort = FALSE
)
paper_results <- merge(
  paper_results, paper_scenario_order,
  by = "scenario", all = FALSE, sort = FALSE
)
setorder(paper_results, paper_scenario_order, paper_method_order)
paper_results <- paper_results[, .(
  scenario, scenario_label = paper_scenario_label, paper_scenario_latex,
  scenario_order = paper_scenario_order, method, method_label,
  method_short, method_order = paper_method_order,
  mean_relative_hs_error = mean_rel_hs_cov,
  sd_relative_hs_error = sd_rel_hs_cov,
  mean_subspace_loss = mean_sin2,
  sd_subspace_loss = sd_mean_sin2,
  mean_runtime_seconds = mean_runtime,
  sd_runtime_seconds = sd_runtime
)]

publication_table <- paper_results[, .(
  scenario_order, scenario, scenario_label, method_order, method, method_label,
  relative_hs_error = format_mean_sd(
    mean_relative_hs_error, sd_relative_hs_error
  ),
  mean_sin2_angle = format_mean_sd(mean_subspace_loss, sd_subspace_loss),
  runtime_seconds = format_mean_sd(
    mean_runtime_seconds, sd_runtime_seconds, digits = args$runtime_digits
  )
)]

primary_wide <- dcast(
  publication_table,
  scenario_order + scenario + scenario_label ~ method_label,
  value.var = "relative_hs_error"
)
setorder(primary_wide, scenario_order)

completion_by_scenario_group <- task_audit[
  , .(
    expected_tasks = .N,
    result_files = sum(result_file_present),
    scheduler_ok = sum(scheduler_ok),
    scheduler_failed = sum(!scheduler_ok)
  ),
  by = .(scenario_order, scenario, scenario_label, group_order, method_group)
]
method_group_completion <- method_audit[
  , .(
    expected_method_fits = .N,
    successful_method_fits = sum(method_complete),
    failed_method_fits = sum(rows_found > 0L & observed_status != "ok"),
    missing_method_fits = sum(rows_found == 0L)
  ),
  by = .(scenario, method_group)
]
completion_by_scenario_group <- merge(
  completion_by_scenario_group, method_group_completion,
  by = c("scenario", "method_group"), all.x = TRUE, sort = FALSE
)
setorder(completion_by_scenario_group, scenario_order, group_order)

completion_overview <- data.table(
  profile = args$profile,
  study_version = args$study_version,
  expected_repetitions = expected_repetitions,
  expected_tasks = nrow(expected_tasks),
  result_files_present = sum(task_audit$result_file_present),
  scheduler_ok = sum(task_audit$scheduler_ok),
  scheduler_not_ok = sum(!task_audit$scheduler_ok),
  expected_method_fits = nrow(expected_methods),
  successful_method_fits = sum(method_audit$method_complete),
  failed_method_fits = nrow(method_failures),
  missing_method_fits = nrow(missing_method_rows),
  duplicate_method_rows = nrow(duplicate_method_rows),
  configuration_mismatches = nrow(configuration_mismatches),
  malformed_result_rows = nrow(malformed_result_rows),
  result_rows_outside_profile = nrow(unexpected_result_rows)
)
completion_overview[, complete :=
  scheduler_not_ok == 0L & failed_method_fits == 0L &
    missing_method_fits == 0L & duplicate_method_rows == 0L &
    configuration_mismatches == 0L & malformed_result_rows == 0L]

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(completion_overview, file.path(args$output_dir, "completion_overview.csv"))
fwrite(task_audit, file.path(args$output_dir, "task_audit.csv"))
fwrite(method_audit, file.path(args$output_dir, "method_audit.csv"))
fwrite(completion_by_scenario_group,
       file.path(args$output_dir, "completion_by_scenario_group.csv"))
fwrite(missing_tasks, file.path(args$output_dir, "missing_or_failed_tasks.csv"))
fwrite(missing_method_rows, file.path(args$output_dir, "missing_method_rows.csv"))
fwrite(method_failures, file.path(args$output_dir, "method_failures.csv"))
fwrite(duplicate_method_rows, file.path(args$output_dir, "duplicate_method_rows.csv"))
fwrite(configuration_mismatches,
       file.path(args$output_dir, "configuration_mismatches.csv"))
fwrite(malformed_result_rows, file.path(args$output_dir, "malformed_result_rows.csv"))
fwrite(unexpected_result_rows,
       file.path(args$output_dir, "result_rows_outside_profile.csv"))
fwrite(comparison, file.path(args$output_dir, "method_comparison_complete.csv"))
fwrite(publication_table, file.path(args$output_dir, "method_comparison_formatted.csv"))
fwrite(primary_wide, file.path(args$output_dir, "relative_hs_error_wide.csv"))
fwrite(paper_results, file.path(args$output_dir, "paper_results_long.csv"))

# Paper tables put methods in rows and use metric-by-setting columns. Separate
# experiment series stay in separate tables so all values remain readable at
# ordinary page width. Standard deviations remain in paper_results_long.csv.
format_table_mean <- function(value, digits, bold = FALSE) {
  out <- rep("--", length(value))
  available <- is.finite(value)
  out[available] <- sprintf(paste0("%.", digits, "f"), value[available])
  emphasize <- rep_len(bold, length(value)) & available
  out[emphasize] <- paste0("\\textbf{", out[emphasize], "}")
  out
}

source(file.path("simulation", "latex_table_helper.R"))
build_paper_series_table <- function(
  table_stem, scenario_ids, setting_labels, setting_latex,
  caption, label, include_time = TRUE
) {
  stopifnot(
    length(scenario_ids) == length(setting_labels),
    length(scenario_ids) == length(setting_latex)
  )
  metrics <- data.table(
    metric = c("Cov", "Subspace", "Time"),
    value_column = c(
      "mean_relative_hs_error", "mean_subspace_loss", "mean_runtime_seconds"
    ),
    digits = c(min(args$digits, 3L), min(args$digits, 3L), args$runtime_digits),
    metric_order = 1:3
  )
  if (!include_time) metrics <- metrics[metric != "Time"]
  settings <- data.table(
    scenario = scenario_ids,
    setting = setting_labels,
    setting_latex = setting_latex,
    setting_order = seq_along(scenario_ids)
  )
  series <- merge(
    paper_results[scenario %in% scenario_ids], settings,
    by = "scenario", all = FALSE, sort = FALSE
  )
  long <- rbindlist(lapply(seq_len(nrow(metrics)), function(index) {
    specification <- metrics[index]
    series[, .(
      method_order, Method = method_short, setting, setting_latex,
      setting_order,
      metric = specification$metric,
      metric_order = specification$metric_order,
      digits = specification$digits,
      value = get(specification$value_column)
    )]
  }))
  long[, best := {
    candidates <- value[is.finite(value)]
    is.finite(value) & value <= min(candidates) + sqrt(.Machine$double.eps)
  }, by = .(metric_order, setting_order)]
  long[, `:=`(
    display = format_table_mean(value, digits),
    display_latex = format_table_mean(value, digits, bold = best)
  )]
  long[, column_id := paste0("m", metric_order, "s", setting_order)]

  column_order <- long[
    order(metric_order, setting_order), unique(column_id)
  ]
  wide <- dcast(
    long, method_order + Method ~ column_id,
    value.var = "display_latex", fill = "--"
  )
  csv_wide <- dcast(
    long, method_order + Method ~ column_id,
    value.var = "display", fill = "--"
  )
  setcolorder(wide, c("method_order", "Method", column_order))
  setcolorder(csv_wide, c("method_order", "Method", column_order))
  setorder(wide, method_order)
  setorder(csv_wide, method_order)

  csv_table <- copy(csv_wide[, !"method_order"])
  flat_names <- unlist(lapply(metrics$metric, function(metric_name) {
    paste(metric_name, setting_labels, sep = ": ")
  }))
  setnames(csv_table, c("Method", flat_names))
  fwrite(csv_table, file.path(args$output_dir, paste0(table_stem, ".csv")))

  metric_groups <- setNames(
    rep(length(scenario_ids), nrow(metrics)), metrics$metric
  )
  latex_table <- latex_result_table(
    wide[, !"method_order"],
    settings_cols = "Method",
    metric_groups = metric_groups,
    column_names = c("Method", rep(setting_latex, times = nrow(metrics))),
    caption = caption,
    label = label,
    align = c("l", rep("c", length(column_order))),
    scale_down = TRUE,
    font_size = 9,
    collapse_cols = 1L,
    latex_hline = "none"
  )
  write_latex_result_table(
    latex_table, file.path(args$output_dir, paste0(table_stem, ".txt"))
  )
  invisible(list(data = csv_table, latex = latex_table))
}

dimension_table <- build_paper_series_table(
  "simulation_dimension_table",
  c("dim_p8", "dim_p16", "dim_p32"),
  c("8 x 8", "16 x 16", "32 x 32"),
  c("$8 \\times 8$", "$16 \\times 16$", "$32 \\times 32$"),
  paste(
    "Performance as image dimension varies. Cov denotes relative squared",
    "Hilbert--Schmidt covariance error, Subspace denotes mean squared sine of",
    "the principal angles, and Time is fitting time in seconds. Entries are",
    "Monte Carlo means; smaller is better, bold is best in column, and --",
    "denotes a method not run."
  ),
  "tab:simulation-dimension"
)
rank_table <- build_paper_series_table(
  "simulation_rank_table",
  c("rank6_p8", "dim_p8", "rank12_p8"),
  c("K = 6", "K = 9", "K = 12"),
  c("$K=6$", "$K=9$", "$K=12$"),
  paste(
    "Performance as covariance reconstruction rank varies. Cov denotes",
    "relative squared Hilbert--Schmidt covariance error, Subspace denotes mean",
    "squared sine of the principal angles, and Time is fitting time in seconds.",
    "Entries are Monte Carlo means; smaller is better and bold is best in",
    "column."
  ),
  "tab:simulation-rank"
)
sample_size_table <- build_paper_series_table(
  "simulation_sample_size_table",
  c("n25_p8", "dim_p8", "n100_p8"),
  c("n = 25", "n = 50", "n = 100"),
  c("$n=25$", "$n=50$", "$n=100$"),
  paste(
    "Performance as sample size varies. Cov denotes relative squared",
    "Hilbert--Schmidt covariance error, Subspace denotes mean squared sine of",
    "the principal angles, and Time is fitting time in seconds. Entries are",
    "Monte Carlo means; smaller is better and bold is best in column."
  ),
  "tab:simulation-sample-size"
)
sensitivity_table <- build_paper_series_table(
  "simulation_sensitivity_table",
  c(
    "dim_p8", "visits_20_40_p8", "noise_05_p8", "spectral_p8", "t5_p8",
    "dependent_p8"
  ),
  c("Base", "Visits", "Noise", "Spectrum", "t5", "Dependent"),
  c("Base", "Visits", "Noise", "Spectrum", "$t_5$", "Dep."),
  paste(
    "Sensitivity to observation density and distributional settings. Base is",
    "the primary setting; Visits uses 11--20 rather than 6--10 visits; Noise",
    "uses $\\sigma=0.5$; Spectrum uses the alternative spectrum; $t_5$ uses",
    "heavy-tailed scores; and Dep. uses dependent scores. Cov denotes relative",
    "squared Hilbert--Schmidt covariance error and Subspace denotes mean squared",
    "sine of the principal angles. Entries are Monte Carlo means; smaller is",
    "better and bold is best in column."
  ),
  "tab:simulation-sensitivity",
  include_time = FALSE
)

# Remove outputs from the superseded methods-as-columns layout.
unlink(file.path(
  args$output_dir,
  c(
    "paper_accuracy_table.csv", "paper_runtime_table.csv",
    "simulation_accuracy_table.txt", "simulation_runtime_table.txt"
  )
))

options(width = max(getOption("width"), 200L))
print(completion_overview)
if (nrow(method_failures)) {
  cat("\nMethod-level failures:\n")
  print(method_failures[, .(scenario, seed, method_label, observed_status, error)])
}
cat("\nWrote audit and comparison tables to ", args$output_dir, "\n", sep = "")

# Generate console tables from the same in-memory comparison used for the CSV
# outputs so the two representations cannot drift apart.
split_argument <- function(value) {
  if (is.null(value) || !nzchar(trimws(value))) return(character())
  unique(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
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

printed_comparison <- comparison[method != "spline_pca"]
printed_comparison[, console_method_order := match(method, paper_method_order$method)]
setorder(printed_comparison, scenario_order, console_method_order)
requested_scenarios <- split_argument(args$scenarios)
if (length(requested_scenarios)) {
  unknown_scenarios <- setdiff(
    requested_scenarios,
    c(unique(comparison$scenario), unique(comparison$scenario_label))
  )
  if (length(unknown_scenarios)) {
    stop("Unknown scenario(s): ", paste(unknown_scenarios, collapse = ", "))
  }
  printed_comparison <- printed_comparison[
    scenario %in% requested_scenarios | scenario_label %in% requested_scenarios
  ]
}

if (!isTRUE(args$quiet)) {
  for (scenario_id in unique(printed_comparison$scenario)) {
    section <- printed_comparison[scenario == scenario_id]
    output <- section[, .(
      Method = method_label
    )]
    for (metric_index in seq_len(nrow(metric_definitions))) {
      definition <- metric_definitions[metric_index]
      metric_digits <- if (definition$metric == "runtime") {
        args$runtime_digits
      } else {
        args$digits
      }
      output[, (definition$label) := format_mean_sd(
        section[[definition$mean_column]],
        section[[definition$sd_column]],
        metric_digits
      )]
    }
    setcolorder(output, c("Method", metric_definitions$label))
    cat("\n", section$scenario_label[[1L]], " [", scenario_id, "]\n", sep = "")
    print(output, row.names = FALSE)
  }
  cat(
    "\nEntries are Monte Carlo mean (SD) over successful fits; -- denotes ",
    "an unavailable metric.\n",
    sep = ""
  )
}

if (isTRUE(args$fail_on_incomplete) && !completion_overview$complete) {
  quit(save = "no", status = 1L)
}
