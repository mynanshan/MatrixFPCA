# MatrixFPCA

## Method-comparison simulation

[`experiment/simu-methods.R`](experiment/simu-methods.R) benchmarks MatFPCA against HAPP–PACE/FACE, mFACEs, mGSFPCA, Shi et al.'s longitudinal-image MD-FPCA, SPACO-minus, sparse Product FPCA, and a spline-PCA baseline. Package and research-prototype adapters live in [`experiment/competitors/`](experiment/competitors/). The statistical formulation, mathematical descriptions of every competitor, identifiability-safe metrics, preliminary local results, and commands are in [`experiment/README-method-comparison.md`](experiment/README-method-comparison.md). Repeated Compute Canada experiments are submitted through [`experiment/run_simu_methods.sh`](experiment/run_simu_methods.sh).

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
