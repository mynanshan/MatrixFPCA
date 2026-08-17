source(file.path("application", "fluid", "fluid_flow_utils.R"))

grid <- fluid_spectral_grid(16L)
stopifnot(
  length(grid$x) == 16L,
  all(dim(grid$k2) == c(16L, 16L)),
  all(grid$dealias[abs(grid$kx) > floor(16L / 3L)] == FALSE)
)

subject <- draw_fluid_subject(42L)
times <- seq(0, 0.5, length.out = 6L)
solution <- simulate_fluid_subject(
  subject, output_times = times, solver_dt = 0.01, spectral_grid = grid
)
stopifnot(
  identical(dim(solution$vorticity), c(16L, 16L, 6L)),
  all(is.finite(solution$vorticity)),
  max(abs(solution$diagnostics[, "circulation"])) < 1e-10,
  all(diff(solution$diagnostics[, "enstrophy"]) <= 1e-8),
  solution$diagnostics[6L, "kinetic_energy"] <
    solution$diagnostics[1L, "kinetic_energy"]
)

cohort_one <- simulate_fluid_cohort(
  4L, 12L, seq(0, 0.4, length.out = 5L), 0.02, seed = 99L,
  progress = FALSE
)
cohort_two <- simulate_fluid_cohort(
  4L, 12L, seq(0, 0.4, length.out = 5L), 0.02, seed = 99L,
  progress = FALSE
)
stopifnot(
  identical(cohort_one$latent, cohort_two$latent),
  identical(cohort_one$parameters, cohort_two$parameters),
  length(unique(cohort_one$parameters$seed)) == 4L
)

visits <- sample_fluid_visits(
  n_subjects = 20L, n_time = 21L,
  visits_min = 5L, visits_max = 8L, seed = 100L
)
stopifnot(
  all(lengths(visits) >= 5L),
  all(lengths(visits) <= 8L),
  all(vapply(visits, function(index) min(index) <= ceiling(0.2 * 21), logical(1L))),
  all(vapply(visits, function(index) max(index) >= floor(0.8 * 21), logical(1L)))
)

observed_fields <- lapply(seq_len(4L), function(subject_index) {
  index <- visits[[subject_index]]
  field <- cohort_one$latent[subject_index, , , pmin(index, 5L), drop = FALSE]
  dim(field) <- c(12L, 12L, length(index))
  field
})
mean_field <- estimate_fluid_mean(
  observed_fields, visits[seq_len(4L)], seq(0, 1, length.out = 21L),
  mean_df = 4L
)
stopifnot(
  identical(dim(mean_field), c(12L, 12L, 21L)),
  all(is.finite(mean_field))
)

cat("Fluid-flow solver and sampling checks passed.\n")
