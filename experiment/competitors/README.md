# Competitor implementations

`core.R` contains the common Hilbert-space metrics plus the maintained
HAPP--PACE and spline-PCA implementations. Every adapter returns the common
latent-covariance eigenfit used by
`experiment/simu-methods.R`: `psi[time,row,column,component]`, decreasing
eigenvalues `nu`, estimated noise `sig2`, and `tgrid`. The benchmark then
canonicalizes the functions in its trapezoidal-time/pixel-sum inner product.

## Primary competitors

- `shi-mdfpca.R`: independent, paper-faithful implementation of Shi, Jiang,
  Ma, Beg and Cao (2024). It directly targets sparsely observed longitudinal
  images with sequential alternating regression and a three-domain
  tensor-product B-spline. It adds a documented noise correction when turning
  fitted subject scores into a covariance estimate. No public author code for
  the longitudinal-image estimator was found; this is therefore not described
  as the authors' original implementation.
- `mfaces.R`: the CRAN `mfaces` implementation of Li, Xiao and Luo's smoothed
  multivariate auto/cross-covariance estimator. It is the strongest direct
  nonseparable covariance-smoothing baseline here, but fits `d(d+1)/2`
  surfaces for `d=p*q` pixels. A default safety limit prevents accidental
  large-image runs.
- `happ-face.R` and the existing `fit_happ_pace`: the Happ--Greven finite
  expansion, with FACE or PACE univariate FPCA at each pixel. Comparing these
  two diagnoses whether poor performance is caused by PACE rather than the
  multivariate construction.
- `mgsfpca.R`: the CRAN `mGSFPCA` Gaussian sparse multivariate FPCA pipeline.
  A single rank/basis value is the reproducible fast profile; comma-separated
  candidate lists enable the package's per-pixel AIC selection.

## Secondary structural competitors

- `spaco-minus.R`: a self-contained smooth probabilistic-PARAFAC-style
  comparator without auxiliary covariates. The paper authors' public
  `SPACOexperments` repository imports a `spaco` module that it does not ship,
  while the `spaco==1.0` PyPI wheel contains metadata but no Python module.
  This adapter is consequently labeled an independent SPACO- style
  implementation, not exact author software.
- `product-fpca.R`: CRAN `fdarep::SparseProductFPCA`. It treats a vectorized
  image as a one-dimensional spatial function and records that artificial
  ordering in diagnostics. It is useful as a structured sparse baseline but
  does not respect the two image axes.
- `spline_pca` in `core.R`: smooth-then-PCA, retained as a deliberately
  weak sparse-design baseline.

## Installation and complexity

Run from the repository root:

```bash
Rscript experiment/competitors/install-dependencies.R
```

The current verified package layer is `mfaces 0.1-4`, `mGSFPCA 0.2.2`, and
`fdarep 0.1.1`. Exact versions used in any run are retained inside fit
diagnostics. Package binaries and virtual environments are not copied into the
repository.

Approximate scaling in the number of pixel outcomes `d=p*q`:

| method | important fit-time/memory driver |
|---|---|
| HAPP--PACE/FACE | `d` independent univariate smoothers, then PCA of local scores |
| mFACEs | `d(d+1)/2` smoothed covariance surfaces |
| mGSFPCA | `d` univariate generalized sparse FPCAs, then score PCA |
| Shi MD-FPCA | alternating solves in `b_t*b_r*b_c` spline coefficients |
| SPACO- style | rank-`K` ALS; linear in observed tensor entries per sweep |
| SparseProductFPCA | spatial FPCA plus sparse temporal FPCA |

For dimension scaling, run separate `p=q=3,4,6,8,...` jobs and keep failures
or safety-limit exits in the results table; they are part of the computational
comparison rather than missing data to discard.
