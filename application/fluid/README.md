# 2D fluid-flow MatFPCA application

[`fluid_flow_matfpca.R`](fluid_flow_matfpca.R) creates independent fluid-flow
subjects, hides most of each subject's time snapshots, fits MatFPCA to the
remaining irregular observations, and evaluates recovery against the complete
latent trajectory.

## PDE and numerical experiment

The simulated field is vorticity on the periodic square
\(\Omega=[0,2\pi)^2\). Each subject solves the unforced, incompressible 2D
Navier--Stokes equation

\[
\partial_t\omega + u\,\partial_x\omega + v\,\partial_y\omega
=\nu\Delta\omega,\qquad
\omega=-\Delta\psi,\quad u=\partial_y\psi,\quad v=-\partial_x\psi.
\]

This vorticity-snapshot setup is common in dynamic mode decomposition studies,
while the double-shear-layer initial condition produces more genuinely
nonlinear evolution than an exact single-mode solution such as the
Taylor--Green vortex. The solver uses Fourier pseudo-spectral derivatives,
two-thirds dealiasing, and fourth-order Runge--Kutta time stepping.

Subject \(i\) independently receives:

- log-uniform viscosity \(\nu_i\), reported together with the dimensionless
  convention \(Re_i=1/\nu_i\);
- a random shear-layer thickness and transverse perturbation amplitude;
- independent uniform translations in both periodic spatial directions;
- independent random low-frequency vorticity perturbations.

The translations make the population vorticity mean spatially homogeneous,
but the realized finite cohort is still centered using a subject-weighted
temporal spline estimated from observed snapshots only. Gaussian pixel noise
is added only to the observed data. Each subject is forced to have at least one
early and one late observation, with all remaining visits sampled without
replacement. The unsupplied times are never used in mean estimation, covariance
estimation, rank selection, or score estimation.

MatFPCA is fitted to complete vorticity matrices at those subject-specific time
indices. Marginal ranks are selected by PVE and bounded by configurable caps.
Subject scores inferred from the irregular noisy visits reconstruct the full
dense trajectory. This assesses interpolation/denoising of a random process;
it is not a numerical PDE forecaster and does not give MatFPCA the viscosity or
initial-condition parameters.

## Running it

Run the default experiment from the repository root:

```bash
Rscript application/fluid/fluid_flow_matfpca.R
```

The defaults use 40 subjects, a 24 by 24 grid, 41 dense evaluation times, and
8--14 observed visits per subject. A quick end-to-end check is:

```bash
Rscript application/fluid/fluid_flow_matfpca.R \
  --n-subjects 8 \
  --grid-size 8 \
  --n-time 15 \
  --t-end 1 \
  --solver-dt 0.02 \
  --visits-min 4 \
  --visits-max 6 \
  --weight-scheme identity \
  --nbasis 5 \
  --max-spatial-rank 4 \
  --max-temporal-rank 4 \
  --output-dir /tmp/matfpca_fluid_smoke
```

Use `--help` for all controls. `--save-truth 1` saves the dense latent and
reconstructed arrays for custom downstream analysis; it is off by default to
keep the output compact. The default result directory is
`application/fluid/results`.

## Recovery outputs

The result directory contains:

- `all_subject_dynamics_recovery.pdf`: one page per subject, comparing latent
  vorticity, MatFPCA recovery, and signed error at four times. Panel titles say
  whether a time was observed or held out;
- `subject_*_dynamics_recovery.png`: examples spanning the empirical range of
  held-out recovery error;
- `reconstruction_metrics.csv`: subject-level observed, unobserved, and
  full-grid RMSE, relative RMSE, field correlation, and kinetic-energy and
  enstrophy errors;
- `reconstruction_summary.csv`: cohort summaries of those subject metrics;
- `recovery_by_time.csv` and `recovery_error_by_time.png`: dense-grid recovery
  diagnostics, including every unobserved time;
- `representative_physics_recovery.png`: latent and reconstructed kinetic
  energy and enstrophy for a median-recovery subject;
- `observation_schedule.csv` and `irregular_observation_schedule.png`: the
  exact sparse longitudinal design;
- `subject_parameters.csv`: the independently sampled PDE/initial-condition
  parameters, random-mode coefficients, and reproducibility seeds;
- `fluid_observed_data.rds`: the actual noisy irregular snapshot data as a
  `matfdata.irreg` object, before mean centering;
- `fluid_matfpca_model.rds`: fitted covariance, eigensystem, subject scores,
  mean, visit indices, and configuration;
- `run_summary.txt` and `session_info.txt`: headline recovery numbers, timing,
  and software versions.

The most direct numerical check is the `unobserved` row in
`reconstruction_summary.csv`. The per-subject PDF is intentionally retained as
the visual check that good averages are not hiding subjects whose vortices or
shear-layer deformation are missed.
