# MatrixFPCA

## Implementation

Reusable R implementation modules are kept in [`R/`](R/). Source
[`R/matfpca.R`](R/matfpca.R) from the project root to load the MatFPCA
implementation and its matrix-functional-data dependency.

## ADNI application

[`application/adni/adni_matfpca.R`](application/adni/adni_matfpca.R) applies MatFPCA to
the longitudinal axial ADNI slices, including training-only temporal mean
centering, a reproducible random image holdout, masked reconstruction
diagnostics, PCA/UMAP/t-SNE comparisons, a continuous morphology coordinate
with MRI decoding and age diagnostics, and leading-FPC plots. See
[`application/adni/README.md`](application/adni/README.md) for the design choices,
runtime controls, and generated outputs. The mathematical and computational
details are in [`application/adni/TECHNICAL_DETAILS.md`](application/adni/TECHNICAL_DETAILS.md).

## 2D fluid-flow application

[`application/fluid/fluid_flow_matfpca.R`](application/fluid/fluid_flow_matfpca.R)
simulates iid periodic 2D Navier--Stokes vorticity trajectories with random
viscosities and initial conditions, exposes only irregular noisy snapshots to
MatFPCA, and measures subject-level recovery on the complete latent time grid.
See [`application/fluid/README.md`](application/fluid/README.md) for
the PDE, numerical design, commands, and recovery outputs.

## Method-comparison simulation

[`simulation/simu-methods.R`](simulation/simu-methods.R) benchmarks MatFPCA against HAPP–PACE/FACE, mFACEs, mGSFPCA, Shi et al.'s longitudinal-image MD-FPCA, SPACO-minus, sparse Product FPCA, and a spline-PCA baseline. Package and research-prototype adapters live in [`simulation/competitors/`](simulation/competitors/). The statistical formulation, mathematical descriptions of every competitor, identifiability-safe metrics, preliminary local results, and commands are in [`simulation/README-method-comparison.md`](simulation/README-method-comparison.md). Repeated Compute Canada experiments are submitted through [`simulation/run_simu_methods.sh`](simulation/run_simu_methods.sh).

## Eigenfunction truncation

`trunc.matfeigen()` retains the full row and column eigensystems by default and applies its PVE threshold only to the nonnegative temporal eigenvalues:

```r
eig_trunc <- trunc.matfeigen(eig, pve = 0.95)
```

The earlier greedy rule, which can trim the row, column, and temporal eigensystems together, remains available with `pve_rule = "joint"`:

```r
eig_trunc <- trunc.matfeigen(eig, pve = 0.95, pve_rule = "joint")
```

Supplying all three ranks (`I`, `J`, and `K`) bypasses PVE selection.

## Score inversion

`get_matfpca_scores()` uses a rank-truncated generalized inverse when the estimated noise variance is exactly zero and the usual ridge BLUP when it is positive. This automatic behavior, the generalized-inverse cutoff, and the numerical ridge are configurable:

```r
scores <- get_matfpca_scores(
  data,
  eig,
  sig2 = covariance_fit$sig2,
  K_grid = covariance_fit$K,
  inverse_strategy = "auto", # or "ginv" / "ridge"
  ginv_tol = 1e-3,
  nugget = 1e-8
)
```

The generalized-inverse cutoff is relative to the largest eigenvalue of each subject's observed-time covariance. Setting `inverse_strategy = "ridge"` recovers the preceding score-inversion behavior.
