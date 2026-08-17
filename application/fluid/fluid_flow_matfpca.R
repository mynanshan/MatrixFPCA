#!/usr/bin/env Rscript

parser <- argparse::ArgumentParser(
  description = paste(
    "Simulate iid 2D Navier--Stokes subjects at irregular times and assess",
    "subject-level dynamic recovery with MatFPCA."
  )
)
parser$add_argument("--n-subjects", type = "integer", default = 40L)
parser$add_argument("--grid-size", type = "integer", default = 24L)
parser$add_argument("--n-time", type = "integer", default = 41L)
parser$add_argument("--t-end", type = "double", default = 6)
parser$add_argument("--solver-dt", type = "double", default = 0.01)
parser$add_argument("--visits-min", type = "integer", default = 8L)
parser$add_argument("--visits-max", type = "integer", default = 14L)
parser$add_argument(
  "--noise-fraction", type = "double", default = 0.05,
  help = "Observation-noise SD as a fraction of the latent vorticity SD."
)
parser$add_argument("--viscosity-min", type = "double", default = 0.0025)
parser$add_argument("--viscosity-max", type = "double", default = 0.006)
parser$add_argument("--shear-min", type = "double", default = 0.16)
parser$add_argument("--shear-max", type = "double", default = 0.24)
parser$add_argument("--perturbation-min", type = "double", default = 0.04)
parser$add_argument("--perturbation-max", type = "double", default = 0.08)
parser$add_argument("--random-modes", type = "integer", default = 4L)
parser$add_argument("--random-vorticity-sd", type = "double", default = 0.08)
parser$add_argument(
  "--mean-df", type = "integer", default = 5L,
  help = "Common temporal spline dimension for training-only mean centering."
)
parser$add_argument(
  "--weight-scheme", type = "character", default = "gaussian_refined",
  choices = c("identity", "gaussian_refined")
)
parser$add_argument("--weight-folds", type = "integer", default = 3L)
parser$add_argument("--opw-iter", type = "integer", default = 3L)
parser$add_argument(
  "--nbasis", type = "character", default = "theory",
  help = "Temporal covariance basis dimension or 'theory'."
)
parser$add_argument("--mode-pve", type = "double", default = 0.95)
parser$add_argument("--max-spatial-rank", type = "integer", default = 12L)
parser$add_argument("--max-temporal-rank", type = "integer", default = 10L)
parser$add_argument("--seed", type = "integer", default = 20260817L)
parser$add_argument("--example-subjects", type = "integer", default = 4L)
parser$add_argument("--save-model", type = "integer", default = 1L)
parser$add_argument("--save-truth", type = "integer", default = 0L)
parser$add_argument(
  "--output-dir", type = "character",
  default = file.path("application", "fluid", "results")
)
args <- parser$parse_args()

required_packages <- c(
  "data.table", "fda", "Matrix", "rTensor", "RhpcBLASctl"
)
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1L), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

if (args$n_subjects < 4L) stop("--n-subjects must be at least 4.")
if (args$n_time < 9L) stop("--n-time must be at least 9.")
if (args$t_end <= 0 || args$solver_dt <= 0) {
  stop("--t-end and --solver-dt must be positive.")
}
if (args$viscosity_min <= 0 || args$viscosity_max < args$viscosity_min) {
  stop("Invalid viscosity range.")
}
if (args$shear_min <= 0 || args$shear_max < args$shear_min) {
  stop("Invalid shear-thickness range.")
}
if (args$perturbation_min < 0 ||
    args$perturbation_max < args$perturbation_min) {
  stop("Invalid perturbation range.")
}
if (args$noise_fraction < 0) stop("--noise-fraction cannot be negative.")
if (args$random_modes < 0L || args$random_vorticity_sd < 0) {
  stop("--random-modes and --random-vorticity-sd cannot be negative.")
}
if (args$mode_pve <= 0 || args$mode_pve > 1) {
  stop("--mode-pve must lie in (0, 1].")
}

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path("R", "matfdata.R"))
source(file.path("R", "matfpca.R"))
source(file.path("application", "fluid", "fluid_flow_utils.R"))

available_cores <- if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK"))) {
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
} else {
  min(4L, RhpcBLASctl::get_num_procs())
}
RhpcBLASctl::blas_set_num_threads(max(1L, available_cores))

log_message <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
}

time_grid <- seq(0, args$t_end, length.out = args$n_time)
log_message(
  "Simulating ", args$n_subjects, " iid Navier--Stokes subjects on a ",
  args$grid_size, " x ", args$grid_size, " grid..."
)
simulation_timing <- system.time({
  cohort <- simulate_fluid_cohort(
    n_subjects = args$n_subjects,
    grid_size = args$grid_size,
    output_times = time_grid,
    solver_dt = args$solver_dt,
    seed = args$seed,
    viscosity_range = c(args$viscosity_min, args$viscosity_max),
    shear_range = c(args$shear_min, args$shear_max),
    perturbation_range = c(
      args$perturbation_min, args$perturbation_max
    ),
    random_modes = args$random_modes,
    random_vorticity_sd = args$random_vorticity_sd,
    progress = TRUE
  )
})
latent <- cohort$latent
signal_sd <- stats::sd(as.numeric(latent))
noise_sd <- args$noise_fraction * signal_sd

visit_indices <- sample_fluid_visits(
  args$n_subjects, args$n_time,
  args$visits_min, args$visits_max,
  seed = args$seed + 1L
)
set.seed(args$seed + 2L)
observed_fields <- lapply(seq_len(args$n_subjects), function(subject_index) {
  indices <- visit_indices[[subject_index]]
  field <- latent[subject_index, , , indices, drop = FALSE]
  dim(field) <- c(args$grid_size, args$grid_size, length(indices))
  if (noise_sd > 0) {
    field <- field + stats::rnorm(length(field), 0, noise_sd)
  }
  field
})
raw_mat_data <- create_matrix_funcdata_irreg(
  X_list = lapply(observed_fields, rTensor::as.tensor),
  t_grid = time_grid,
  t_idx_list = lapply(visit_indices, as.integer)
)

log_message("Estimating and removing the common temporal mean from observed visits...")
mean_field <- estimate_fluid_mean(
  observed_fields, visit_indices, time_grid, mean_df = args$mean_df
)
centered_observations <- lapply(seq_len(args$n_subjects), function(subject_index) {
  indices <- visit_indices[[subject_index]]
  centered <- observed_fields[[subject_index]] - mean_field[, , indices, drop = FALSE]
  rTensor::as.tensor(centered)
})
mat_data <- create_matrix_funcdata_irreg(
  X_list = centered_observations,
  t_grid = time_grid,
  t_idx_list = lapply(visit_indices, as.integer)
)

fit_nbasis <- if (grepl("^[0-9]+$", args$nbasis)) {
  as.integer(args$nbasis)
} else {
  args$nbasis
}
n_folds <- if (args$weight_scheme == "identity") {
  1L
} else {
  min(max(2L, args$weight_folds), args$n_subjects)
}
estimate_noise <- noise_sd > 0
log_message(
  "Fitting MatFPCA with ", args$weight_scheme, " weights and ",
  n_folds, " subject fold(s)..."
)
fit_timing <- system.time({
  covariance_fit <- fit_matfcov(
    mat_data,
    wScheme = args$weight_scheme,
    noise = estimate_noise,
    nbasis = fit_nbasis,
    opwIter = args$opw_iter,
    temporal_method = "two_step",
    nfolds = n_folds,
    crossfit_diagonal = "refined"
  )
})

rank_for_pve <- function(values, pve, rank_cap) {
  values <- pmax(as.numeric(values), 0)
  if (!length(values) || sum(values) <= 0) {
    stop("A fitted covariance mode has no positive eigenvalues.")
  }
  selected <- which(cumsum(values) / sum(values) >= pve)[[1L]]
  min(selected, max(1L, as.integer(rank_cap)), length(values))
}

eig_full <- eigen.matfcov(covariance_fit$C)
ranks <- c(
  row = rank_for_pve(
    eig_full$values$lam_phi, args$mode_pve, args$max_spatial_rank
  ),
  column = rank_for_pve(
    eig_full$values$lam_psi, args$mode_pve, args$max_spatial_rank
  ),
  time = rank_for_pve(
    eig_full$values$lam_theta, args$mode_pve, args$max_temporal_rank
  )
)
eig <- trunc.matfeigen(
  eig_full,
  I = ranks[["row"]], J = ranks[["column"]], K = ranks[["time"]],
  verbose = FALSE
)
log_message(
  "Estimating subject scores at row/column/time ranks ",
  paste(ranks, collapse = "/"), "..."
)
score_timing <- system.time({
  scores <- get_matfpca_scores(
    mat_data,
    eig,
    sig2 = covariance_fit$C$sig2,
    K_grid = covariance_fit$C$K
  )
  centered_reconstruction <- recover_matfdata(eig$funcs, scores)@data
})
reconstruction <- sweep(
  centered_reconstruction, c(2L, 3L, 4L), mean_field, "+"
)

safe_correlation <- function(x, y) {
  if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x, y))
}

log_message("Evaluating latent-field and physics recovery at every time point...")
time_rows <- vector("list", args$n_subjects * args$n_time)
row_index <- 0L
for (subject_index in seq_len(args$n_subjects)) {
  observed <- seq_len(args$n_time) %in% visit_indices[[subject_index]]
  true_diagnostics <- cohort$diagnostics[[subject_index]]
  for (time_index in seq_len(args$n_time)) {
    row_index <- row_index + 1L
    truth <- latent[subject_index, , , time_index]
    estimate <- reconstruction[subject_index, , , time_index]
    error <- estimate - truth
    estimated_diagnostics <- fluid_diagnostics(
      estimate, cohort$spectral_grid
    )
    truth_norm <- sum(truth^2)
    time_rows[[row_index]] <- data.table::data.table(
      subject = subject_index,
      time_index = time_index,
      time = time_grid[[time_index]],
      split = if (observed[[time_index]]) "observed" else "unobserved",
      rmse = sqrt(mean(error^2)),
      relative_rmse = if (truth_norm > 0) {
        sqrt(sum(error^2) / truth_norm)
      } else {
        NA_real_
      },
      correlation = safe_correlation(as.numeric(truth), as.numeric(estimate)),
      kinetic_energy_true = true_diagnostics[time_index, "kinetic_energy"],
      kinetic_energy_estimate = estimated_diagnostics[["kinetic_energy"]],
      enstrophy_true = true_diagnostics[time_index, "enstrophy"],
      enstrophy_estimate = estimated_diagnostics[["enstrophy"]]
    )
  }
}
recovery_by_time <- data.table::rbindlist(time_rows)

subject_rows <- vector("list", args$n_subjects * 3L)
row_index <- 0L
for (subject_index in seq_len(args$n_subjects)) {
  split_indices <- list(
    observed = visit_indices[[subject_index]],
    unobserved = setdiff(seq_len(args$n_time), visit_indices[[subject_index]]),
    all = seq_len(args$n_time)
  )
  for (split_name in names(split_indices)) {
    row_index <- row_index + 1L
    indices <- split_indices[[split_name]]
    truth <- latent[subject_index, , , indices, drop = FALSE]
    estimate <- reconstruction[subject_index, , , indices, drop = FALSE]
    error <- estimate - truth
    truth_vector <- as.numeric(truth)
    estimate_vector <- as.numeric(estimate)
    physics <- recovery_by_time[
      subject == subject_index & time_index %in% indices
    ]
    subject_rows[[row_index]] <- data.table::data.table(
      subject = subject_index,
      split = split_name,
      n_times = length(indices),
      rmse = sqrt(mean(error^2)),
      normalized_rmse = sqrt(mean(error^2)) / signal_sd,
      relative_rmse = sqrt(sum(error^2) / sum(truth^2)),
      correlation = safe_correlation(truth_vector, estimate_vector),
      kinetic_energy_relative_rmse = sqrt(
        sum((physics$kinetic_energy_estimate - physics$kinetic_energy_true)^2) /
          sum(physics$kinetic_energy_true^2)
      ),
      enstrophy_relative_rmse = sqrt(
        sum((physics$enstrophy_estimate - physics$enstrophy_true)^2) /
          sum(physics$enstrophy_true^2)
      )
    )
  }
}
reconstruction_metrics <- data.table::rbindlist(subject_rows)
reconstruction_metrics <- merge(
  reconstruction_metrics,
  data.table::as.data.table(cohort$parameters),
  by = "subject", sort = FALSE
)
reconstruction_summary <- reconstruction_metrics[, .(
  subjects = .N,
  rmse_mean = mean(rmse),
  normalized_rmse_mean = mean(normalized_rmse),
  relative_rmse_mean = mean(relative_rmse),
  relative_rmse_median = stats::median(relative_rmse),
  relative_rmse_q25 = stats::quantile(relative_rmse, 0.25),
  relative_rmse_q75 = stats::quantile(relative_rmse, 0.75),
  correlation_mean = mean(correlation, na.rm = TRUE),
  kinetic_energy_relative_rmse_mean = mean(kinetic_energy_relative_rmse),
  enstrophy_relative_rmse_mean = mean(enstrophy_relative_rmse)
), by = split]

schedule <- data.table::CJ(
  subject = seq_len(args$n_subjects),
  time_index = seq_len(args$n_time)
)
schedule[, `:=`(
  time = time_grid[time_index],
  observed = mapply(
    function(subject_id, index) index %in% visit_indices[[subject_id]],
    subject, time_index
  )
)]

data.table::fwrite(
  data.table::as.data.table(cohort$parameters),
  file.path(args$output_dir, "subject_parameters.csv")
)
data.table::fwrite(
  schedule, file.path(args$output_dir, "observation_schedule.csv")
)
data.table::fwrite(
  recovery_by_time, file.path(args$output_dir, "recovery_by_time.csv")
)
data.table::fwrite(
  reconstruction_metrics,
  file.path(args$output_dir, "reconstruction_metrics.csv")
)
data.table::fwrite(
  reconstruction_summary,
  file.path(args$output_dir, "reconstruction_summary.csv")
)
saveRDS(
  raw_mat_data,
  file.path(args$output_dir, "fluid_observed_data.rds"),
  compress = FALSE
)

vorticity_palette <- grDevices::colorRampPalette(
  c("#253494", "#2C7FB8", "#F7F7F7", "#F03B20", "#99000D")
)(101L)
error_palette <- grDevices::colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(101L)

draw_fluid_field <- function(field, value_limits, title, palette) {
  graphics::image(
    cohort$spectral_grid$x, cohort$spectral_grid$y, field,
    zlim = value_limits, col = palette, axes = FALSE,
    xlab = "", ylab = "", asp = 1, useRaster = TRUE
  )
  graphics::box()
  graphics::title(main = title, cex.main = 0.78)
}

snapshot_indices <- unique(as.integer(round(seq(
  1L, args$n_time, length.out = 4L
))))
draw_subject_recovery <- function(subject_index) {
  observed <- visit_indices[[subject_index]]
  truth_slice <- latent[subject_index, , , snapshot_indices, drop = FALSE]
  reconstruction_slice <- reconstruction[
    subject_index, , , snapshot_indices, drop = FALSE
  ]
  field_limit <- max(abs(c(truth_slice, reconstruction_slice)))
  error_limit <- max(abs(reconstruction_slice - truth_slice))
  if (error_limit == 0) error_limit <- 1
  field_range <- c(-field_limit, field_limit)
  error_range <- c(-error_limit, error_limit)

  for (time_index in snapshot_indices) {
    visit_label <- if (time_index %in% observed) "observed" else "held out"
    draw_fluid_field(
      latent[subject_index, , , time_index], field_range,
      sprintf("Truth: t = %.2f (%s)", time_grid[[time_index]], visit_label),
      vorticity_palette
    )
  }
  for (time_index in snapshot_indices) {
    draw_fluid_field(
      reconstruction[subject_index, , , time_index], field_range,
      sprintf("MatFPCA: t = %.2f", time_grid[[time_index]]),
      vorticity_palette
    )
  }
  for (time_index in snapshot_indices) {
    draw_fluid_field(
      reconstruction[subject_index, , , time_index] -
        latent[subject_index, , , time_index],
      error_range, sprintf("Error: t = %.2f", time_grid[[time_index]]),
      error_palette
    )
  }
  unobserved_error <- reconstruction_metrics[
    subject == subject_index & split == "unobserved", relative_rmse
  ]
  parameter <- cohort$parameters[subject_index, ]
  graphics::mtext(
    sprintf(
      "Subject %d: Re = %.0f, nu = %.4f, held-out relative RMSE = %.3f",
      subject_index, parameter$reynolds_number, parameter$viscosity,
      unobserved_error
    ),
    side = 3L, outer = TRUE, line = 1, cex = 1.05
  )
}

all_subject_pdf <- file.path(args$output_dir, "all_subject_dynamics_recovery.pdf")
grDevices::pdf(all_subject_pdf, width = 11, height = 8.5, onefile = TRUE)
for (subject_index in seq_len(args$n_subjects)) {
  graphics::par(mfrow = c(3L, 4L), mar = c(0.8, 0.8, 2.2, 0.4),
                oma = c(0, 0, 3, 0))
  draw_subject_recovery(subject_index)
}
grDevices::dev.off()

ordered_subjects <- reconstruction_metrics[
  split == "unobserved"
][order(relative_rmse), subject]
example_count <- min(max(1L, args$example_subjects), args$n_subjects)
example_positions <- unique(as.integer(round(seq(
  1L, length(ordered_subjects), length.out = example_count
))))
example_subjects <- ordered_subjects[example_positions]
for (subject_index in example_subjects) {
  grDevices::png(
    file.path(
      args$output_dir,
      sprintf("subject_%03d_dynamics_recovery.png", subject_index)
    ),
    width = 1800, height = 1300, res = 160
  )
  graphics::par(mfrow = c(3L, 4L), mar = c(0.8, 0.8, 2.2, 0.4),
                oma = c(0, 0, 3, 0))
  draw_subject_recovery(subject_index)
  grDevices::dev.off()
}

grDevices::png(
  file.path(args$output_dir, "irregular_observation_schedule.png"),
  width = 1400, height = 900, res = 160
)
schedule_matrix <- matrix(FALSE, args$n_subjects, args$n_time)
for (subject_index in seq_len(args$n_subjects)) {
  schedule_matrix[subject_index, visit_indices[[subject_index]]] <- TRUE
}
graphics::par(mar = c(4.5, 4.5, 2, 1))
graphics::image(
  time_grid, seq_len(args$n_subjects), t(schedule_matrix),
  col = c("#F2F2F2", "#2166AC"), xlab = "Time", ylab = "Subject",
  main = "Irregular snapshots supplied to MatFPCA", useRaster = TRUE
)
grDevices::dev.off()

grDevices::png(
  file.path(args$output_dir, "recovery_error_by_time.png"),
  width = 1400, height = 900, res = 160
)
time_summary <- recovery_by_time[, .(
  median = stats::median(relative_rmse),
  lower = stats::quantile(relative_rmse, 0.25),
  upper = stats::quantile(relative_rmse, 0.75)
), by = time]
graphics::par(mar = c(4.5, 4.8, 2, 1))
graphics::plot(
  time_summary$time, time_summary$median, type = "n",
  ylim = range(c(time_summary$lower, time_summary$upper)),
  xlab = "Time", ylab = "Subject relative RMSE",
  main = "Recovery over the dense evaluation grid"
)
graphics::polygon(
  c(time_summary$time, rev(time_summary$time)),
  c(time_summary$lower, rev(time_summary$upper)),
  col = grDevices::adjustcolor("#2C7FB8", alpha.f = 0.22), border = NA
)
graphics::lines(time_summary$time, time_summary$median, col = "#045A8D", lwd = 2)
graphics::rug(unlist(lapply(visit_indices, function(index) time_grid[index])),
              col = grDevices::adjustcolor("black", alpha.f = 0.15))
grDevices::dev.off()

representative_subject <- ordered_subjects[[ceiling(length(ordered_subjects) / 2)]]
representative_physics <- recovery_by_time[subject == representative_subject]
grDevices::png(
  file.path(args$output_dir, "representative_physics_recovery.png"),
  width = 1500, height = 750, res = 160
)
graphics::par(mfrow = c(1L, 2L), mar = c(4.2, 4.5, 2.5, 1))
for (quantity in c("kinetic_energy", "enstrophy")) {
  truth_name <- paste0(quantity, "_true")
  estimate_name <- paste0(quantity, "_estimate")
  value_range <- range(representative_physics[[truth_name]],
                       representative_physics[[estimate_name]])
  graphics::plot(
    representative_physics$time, representative_physics[[truth_name]],
    type = "l", lwd = 2, col = "black", ylim = value_range,
    xlab = "Time", ylab = gsub("_", " ", quantity),
    main = paste("Subject", representative_subject)
  )
  graphics::lines(
    representative_physics$time, representative_physics[[estimate_name]],
    lwd = 2, col = "#D7301F", lty = 2
  )
  observed_rows <- representative_physics$split == "observed"
  graphics::points(
    representative_physics$time[observed_rows],
    representative_physics[[truth_name]][observed_rows],
    pch = 16, col = "#2166AC"
  )
  graphics::legend(
    "topright", c("latent", "MatFPCA", "observed visit"),
    col = c("black", "#D7301F", "#2166AC"),
    lty = c(1, 2, NA), pch = c(NA, NA, 16), bty = "n", cex = 0.85
  )
}
grDevices::dev.off()

if (as.logical(args$save_model)) {
  model <- list(
    covariance_fit = covariance_fit,
    eig = eig,
    scores = scores,
    ranks = ranks,
    mean_field = mean_field,
    time_grid = time_grid,
    x_grid = cohort$spectral_grid$x,
    y_grid = cohort$spectral_grid$y,
    visit_indices = visit_indices,
    subject_parameters = cohort$parameters,
    signal_sd = signal_sd,
    noise_sd = noise_sd,
    arguments = args
  )
  saveRDS(
    model, file.path(args$output_dir, "fluid_matfpca_model.rds"),
    compress = FALSE
  )
}
if (as.logical(args$save_truth)) {
  saveRDS(
    list(latent = latent, reconstruction = reconstruction, cohort = cohort),
    file.path(args$output_dir, "fluid_latent_and_reconstruction.rds"),
    compress = FALSE
  )
}

unobserved_summary <- reconstruction_summary[split == "unobserved"]
summary_lines <- c(
  "2D incompressible Navier--Stokes MatFPCA application",
  "Equation: omega_t + u omega_x + v omega_y = nu Delta omega; omega = -Delta psi",
  sprintf("Subjects: %d", args$n_subjects),
  sprintf("Spatial grid: %d x %d", args$grid_size, args$grid_size),
  sprintf("Dense evaluation grid: %d times on [0, %.3f]", args$n_time, args$t_end),
  sprintf("Observed visits per subject: %d--%d", args$visits_min, args$visits_max),
  sprintf("Latent SD / observation noise SD: %.6f / %.6f", signal_sd, noise_sd),
  sprintf("Weight scheme / folds: %s / %d", args$weight_scheme, n_folds),
  sprintf("Retained row/column/time ranks: %s", paste(ranks, collapse = "/")),
  sprintf("Estimated noise variance: %.8f", covariance_fit$C$sig2),
  sprintf("Simulation elapsed seconds: %.3f", simulation_timing[["elapsed"]]),
  sprintf("MatFPCA fit elapsed seconds: %.3f", fit_timing[["elapsed"]]),
  sprintf("Score/reconstruction elapsed seconds: %.3f", score_timing[["elapsed"]]),
  sprintf(
    "Held-out relative RMSE, mean / median / IQR: %.4f / %.4f / [%.4f, %.4f]",
    unobserved_summary$relative_rmse_mean,
    unobserved_summary$relative_rmse_median,
    unobserved_summary$relative_rmse_q25,
    unobserved_summary$relative_rmse_q75
  ),
  sprintf("Held-out field correlation, mean: %.4f", unobserved_summary$correlation_mean),
  sprintf(
    "Held-out energy / enstrophy relative RMSE, mean: %.4f / %.4f",
    unobserved_summary$kinetic_energy_relative_rmse_mean,
    unobserved_summary$enstrophy_relative_rmse_mean
  )
)
writeLines(summary_lines, file.path(args$output_dir, "run_summary.txt"))
capture_session_info <- function() {
  previous_timezone <- Sys.getenv("TZ", unset = NA_character_)
  on.exit({
    if (is.na(previous_timezone)) {
      Sys.unsetenv("TZ")
    } else {
      Sys.setenv(TZ = previous_timezone)
    }
  })
  Sys.setenv(TZ = "UTC")
  capture.output(sessionInfo())
}
writeLines(
  capture_session_info(), file.path(args$output_dir, "session_info.txt")
)

cat(paste(summary_lines, collapse = "\n"), "\n")
log_message("Results written to ", normalizePath(args$output_dir))
