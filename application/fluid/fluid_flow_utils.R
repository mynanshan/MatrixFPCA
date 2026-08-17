# Utilities for the 2D incompressible Navier--Stokes MatFPCA application.
# The spatial domain is the periodic square [0, 2*pi)^2.  Vorticity and the
# stream function follow omega = -Delta(psi), u = psi_y, v = -psi_x.

fluid_spectral_grid <- function(grid_size) {
  grid_size <- as.integer(grid_size)
  if (!is.finite(grid_size) || grid_size < 8L || grid_size %% 2L != 0L) {
    stop("grid_size must be an even integer of at least 8.")
  }

  coordinates <- 2 * pi * (0:(grid_size - 1L)) / grid_size
  modes <- c(0:(grid_size / 2L),
             (-grid_size / 2L + 1L):-1L)
  kx <- outer(modes, rep(1, grid_size))
  ky <- outer(rep(1, grid_size), modes)
  k2 <- kx^2 + ky^2
  dealias_cutoff <- floor(grid_size / 3L)

  list(
    n = grid_size,
    x = coordinates,
    y = coordinates,
    kx = kx,
    ky = ky,
    k2 = k2,
    dealias = abs(kx) <= dealias_cutoff & abs(ky) <= dealias_cutoff
  )
}

fluid_fft2 <- function(field) {
  stats::fft(field)
}

fluid_ifft2_real <- function(coefficients) {
  Re(stats::fft(coefficients, inverse = TRUE) / length(coefficients))
}

fluid_velocity_from_vorticity_hat <- function(omega_hat, spectral_grid) {
  psi_hat <- matrix(0 + 0i, nrow(omega_hat), ncol(omega_hat))
  nonzero <- spectral_grid$k2 > 0
  psi_hat[nonzero] <- omega_hat[nonzero] / spectral_grid$k2[nonzero]

  list(
    u = fluid_ifft2_real(1i * spectral_grid$ky * psi_hat),
    v = fluid_ifft2_real(-1i * spectral_grid$kx * psi_hat)
  )
}

fluid_rhs_hat <- function(omega_hat, viscosity, spectral_grid) {
  omega_hat <- omega_hat * spectral_grid$dealias
  velocity <- fluid_velocity_from_vorticity_hat(omega_hat, spectral_grid)
  omega_x <- fluid_ifft2_real(1i * spectral_grid$kx * omega_hat)
  omega_y <- fluid_ifft2_real(1i * spectral_grid$ky * omega_hat)
  advection <- velocity$u * omega_x + velocity$v * omega_y
  nonlinear_hat <- fluid_fft2(-advection) * spectral_grid$dealias

  nonlinear_hat - viscosity * spectral_grid$k2 * omega_hat
}

fluid_rk4_step <- function(omega_hat, step_size, viscosity, spectral_grid) {
  k1 <- fluid_rhs_hat(omega_hat, viscosity, spectral_grid)
  k2 <- fluid_rhs_hat(
    omega_hat + 0.5 * step_size * k1, viscosity, spectral_grid
  )
  k3 <- fluid_rhs_hat(
    omega_hat + 0.5 * step_size * k2, viscosity, spectral_grid
  )
  k4 <- fluid_rhs_hat(omega_hat + step_size * k3, viscosity, spectral_grid)
  updated <- omega_hat + step_size * (k1 + 2 * k2 + 2 * k3 + k4) / 6
  updated * spectral_grid$dealias
}

fluid_diagnostics <- function(vorticity, spectral_grid) {
  omega_hat <- fluid_fft2(vorticity) * spectral_grid$dealias
  velocity <- fluid_velocity_from_vorticity_hat(omega_hat, spectral_grid)
  c(
    kinetic_energy = 0.5 * mean(velocity$u^2 + velocity$v^2),
    enstrophy = 0.5 * mean(vorticity^2),
    circulation = mean(vorticity)
  )
}

draw_fluid_subject <- function(
  seed,
  viscosity_range = c(0.0025, 0.006),
  shear_range = c(0.16, 0.24),
  perturbation_range = c(0.04, 0.08),
  random_modes = 4L,
  random_vorticity_sd = 0.08
) {
  set.seed(seed)
  random_modes <- as.integer(random_modes)
  viscosity <- exp(stats::runif(
    1L, log(viscosity_range[[1L]]), log(viscosity_range[[2L]])
  ))
  shear_thickness <- stats::runif(1L, shear_range[[1L]], shear_range[[2L]])
  perturbation <- stats::runif(
    1L, perturbation_range[[1L]], perturbation_range[[2L]]
  )

  list(
    seed = as.integer(seed),
    viscosity = viscosity,
    reynolds_number = 1 / viscosity,
    shear_thickness = shear_thickness,
    perturbation = perturbation,
    x_shift = stats::runif(1L, 0, 2 * pi),
    y_shift = stats::runif(1L, 0, 2 * pi),
    mode_kx = sample.int(3L, random_modes, replace = TRUE),
    mode_ky = sample.int(3L, random_modes, replace = TRUE),
    mode_phase_x = stats::runif(random_modes, 0, 2 * pi),
    mode_phase_y = stats::runif(random_modes, 0, 2 * pi),
    mode_amplitude = stats::rnorm(random_modes, 0, random_vorticity_sd)
  )
}

double_shear_vorticity <- function(subject, spectral_grid) {
  x_matrix <- outer(spectral_grid$x, rep(1, spectral_grid$n))
  y_matrix <- outer(rep(1, spectral_grid$n), spectral_grid$y)
  shifted_x <- (x_matrix - subject$x_shift) %% (2 * pi)
  shifted_y <- (y_matrix - subject$y_shift) %% (2 * pi)
  lower_layer <- shifted_y <= pi

  first_argument <- (shifted_y - pi / 2) / subject$shear_thickness
  second_argument <- (3 * pi / 2 - shifted_y) / subject$shear_thickness
  du_dy <- ifelse(
    lower_layer,
    1 / (subject$shear_thickness * cosh(first_argument)^2),
    -1 / (subject$shear_thickness * cosh(second_argument)^2)
  )
  dv_dx <- subject$perturbation * cos(shifted_x)
  vorticity <- dv_dx - du_dy

  if (length(subject$mode_amplitude)) {
    for (component in seq_along(subject$mode_amplitude)) {
      vorticity <- vorticity + subject$mode_amplitude[[component]] *
        sin(subject$mode_kx[[component]] * shifted_x +
              subject$mode_phase_x[[component]]) *
        sin(subject$mode_ky[[component]] * shifted_y +
              subject$mode_phase_y[[component]])
    }
  }

  # Remove roundoff in the zero Fourier mode.  Periodic vorticity must have
  # zero spatial integral for the stream-function inversion used here.
  vorticity - mean(vorticity)
}

simulate_fluid_subject <- function(
  subject,
  output_times,
  solver_dt,
  spectral_grid
) {
  output_times <- as.numeric(output_times)
  if (length(output_times) < 2L || output_times[[1L]] != 0 ||
      is.unsorted(output_times, strictly = TRUE)) {
    stop("output_times must be strictly increasing and start at zero.")
  }
  if (!is.finite(solver_dt) || solver_dt <= 0) {
    stop("solver_dt must be positive and finite.")
  }

  n_time <- length(output_times)
  vorticity <- array(NA_real_, c(spectral_grid$n, spectral_grid$n, n_time))
  diagnostics <- matrix(
    NA_real_, n_time, 3L,
    dimnames = list(NULL, c("kinetic_energy", "enstrophy", "circulation"))
  )

  initial <- double_shear_vorticity(subject, spectral_grid)
  omega_hat <- fluid_fft2(initial) * spectral_grid$dealias
  vorticity[, , 1L] <- fluid_ifft2_real(omega_hat)
  diagnostics[1L, ] <- fluid_diagnostics(vorticity[, , 1L], spectral_grid)
  current_time <- 0

  for (time_index in 2:n_time) {
    target_time <- output_times[[time_index]]
    while (current_time < target_time - 10 * .Machine$double.eps) {
      step_size <- min(solver_dt, target_time - current_time)
      omega_hat <- fluid_rk4_step(
        omega_hat, step_size, subject$viscosity, spectral_grid
      )
      current_time <- current_time + step_size
    }
    vorticity[, , time_index] <- fluid_ifft2_real(omega_hat)
    diagnostics[time_index, ] <- fluid_diagnostics(
      vorticity[, , time_index], spectral_grid
    )
  }

  if (any(!is.finite(vorticity))) {
    stop("The Navier--Stokes solver produced non-finite vorticity values.")
  }
  list(vorticity = vorticity, diagnostics = diagnostics)
}

simulate_fluid_cohort <- function(
  n_subjects,
  grid_size,
  output_times,
  solver_dt,
  seed,
  viscosity_range = c(0.0025, 0.006),
  shear_range = c(0.16, 0.24),
  perturbation_range = c(0.04, 0.08),
  random_modes = 4L,
  random_vorticity_sd = 0.08,
  progress = interactive()
) {
  n_subjects <- as.integer(n_subjects)
  if (!is.finite(n_subjects) || n_subjects < 2L) {
    stop("n_subjects must be at least two.")
  }
  spectral_grid <- fluid_spectral_grid(grid_size)
  set.seed(seed)
  subject_seeds <- sample.int(.Machine$integer.max, n_subjects)
  subjects <- lapply(subject_seeds, draw_fluid_subject,
    viscosity_range = viscosity_range,
    shear_range = shear_range,
    perturbation_range = perturbation_range,
    random_modes = random_modes,
    random_vorticity_sd = random_vorticity_sd
  )

  latent <- array(
    NA_real_, c(n_subjects, grid_size, grid_size, length(output_times))
  )
  diagnostics <- vector("list", n_subjects)
  for (subject_index in seq_len(n_subjects)) {
    if (progress) {
      message("Simulating fluid subject ", subject_index, "/", n_subjects, "...")
    }
    solution <- simulate_fluid_subject(
      subjects[[subject_index]], output_times, solver_dt, spectral_grid
    )
    latent[subject_index, , , ] <- solution$vorticity
    diagnostics[[subject_index]] <- solution$diagnostics
  }

  parameter_table <- data.frame(
    subject = seq_len(n_subjects),
    seed = vapply(subjects, `[[`, integer(1L), "seed"),
    viscosity = vapply(subjects, `[[`, numeric(1L), "viscosity"),
    reynolds_number = vapply(subjects, `[[`, numeric(1L), "reynolds_number"),
    shear_thickness = vapply(subjects, `[[`, numeric(1L), "shear_thickness"),
    perturbation = vapply(subjects, `[[`, numeric(1L), "perturbation"),
    x_shift = vapply(subjects, `[[`, numeric(1L), "x_shift"),
    y_shift = vapply(subjects, `[[`, numeric(1L), "y_shift"),
    mode_kx = vapply(
      subjects, function(subject) paste(subject$mode_kx, collapse = ";"),
      character(1L)
    ),
    mode_ky = vapply(
      subjects, function(subject) paste(subject$mode_ky, collapse = ";"),
      character(1L)
    ),
    mode_phase_x = vapply(
      subjects, function(subject) paste(subject$mode_phase_x, collapse = ";"),
      character(1L)
    ),
    mode_phase_y = vapply(
      subjects, function(subject) paste(subject$mode_phase_y, collapse = ";"),
      character(1L)
    ),
    mode_amplitude = vapply(
      subjects, function(subject) paste(subject$mode_amplitude, collapse = ";"),
      character(1L)
    )
  )

  list(
    latent = latent,
    diagnostics = diagnostics,
    parameters = parameter_table,
    subjects = subjects,
    spectral_grid = spectral_grid,
    output_times = output_times
  )
}

sample_fluid_visits <- function(
  n_subjects,
  n_time,
  visits_min,
  visits_max,
  seed
) {
  visits_min <- as.integer(visits_min)
  visits_max <- as.integer(visits_max)
  if (visits_min < 3L || visits_max < visits_min || visits_max >= n_time) {
    stop("Require 3 <= visits_min <= visits_max < n_time.")
  }
  set.seed(seed)
  early_pool <- seq_len(max(1L, ceiling(0.2 * n_time)))
  late_pool <- seq.int(min(n_time, floor(0.8 * n_time)), n_time)

  lapply(seq_len(n_subjects), function(subject_index) {
    visit_count <- sample(visits_min:visits_max, 1L)
    anchors <- unique(c(sample(early_pool, 1L), sample(late_pool, 1L)))
    remaining <- setdiff(seq_len(n_time), anchors)
    sort(as.integer(c(
      anchors,
      sample(remaining, visit_count - length(anchors), replace = FALSE)
    )))
  })
}

estimate_fluid_mean <- function(
  observed_fields,
  visit_indices,
  time_grid,
  mean_df = 5L,
  ridge = 1e-8
) {
  n_subjects <- length(observed_fields)
  if (length(visit_indices) != n_subjects) {
    stop("observed_fields and visit_indices must have the same length.")
  }
  spatial_dim <- dim(observed_fields[[1L]])[1:2]
  mean_df <- min(max(as.integer(mean_df), 4L), length(time_grid) - 1L)
  basis_grid <- splines::bs(
    time_grid, df = mean_df, degree = 3L, intercept = TRUE,
    Boundary.knots = range(time_grid)
  )
  stacked_indices <- unlist(visit_indices, use.names = FALSE)
  basis_observed <- basis_grid[stacked_indices, , drop = FALSE]
  subject_weights <- unlist(lapply(
    visit_indices, function(index) rep(1 / length(index), length(index))
  ), use.names = FALSE)
  sqrt_weights <- sqrt(subject_weights)

  response <- do.call(rbind, lapply(observed_fields, function(field) {
    t(matrix(field, nrow = prod(spatial_dim)))
  }))
  weighted_basis <- basis_observed * sqrt_weights
  weighted_response <- response * sqrt_weights
  penalty <- ridge * max(sum(diag(crossprod(weighted_basis))), 1)
  coefficients <- solve(
    crossprod(weighted_basis) + penalty * diag(ncol(weighted_basis)),
    crossprod(weighted_basis, weighted_response)
  )
  fitted_grid <- basis_grid %*% coefficients

  array(t(fitted_grid), dim = c(spatial_dim, length(time_grid)))
}
