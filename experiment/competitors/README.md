# Competitor implementations

This directory isolates package adapters and independent implementations from the simulation driver. The complete statistical formulation and mathematical description of every estimator are in [`../README-method-comparison.md`](../README-method-comparison.md).

Every adapter returns the same latent-covariance eigenfit: `psi[time,row,column,component]`, decreasing eigenvalues `nu`, estimated noise variance `sig2`, and `tgrid`. `core.R` canonicalizes these eigenfunctions in the trapezoidal-time/pixel-sum inner product before evaluation, so signs, numerical scaling, and arbitrary basis rotations do not change the metrics.

## File map

- `core.R` contains shared Hilbert-space metrics, HAPP–PACE, and the deliberately weak smooth-then-PCA spline baseline.
- `happ-face.R` applies the Happ–Greven finite-expansion construction with a FACE rather than PACE univariate decomposition at every pixel.
- `mfaces.R` calls `mfaces::mface.sparse` to estimate all $d(d+1)/2$ multivariate auto/cross-covariance surfaces for $d=pq$ pixels.
- `mgsfpca.R` calls the Gaussian path of `mGSFPCA::spMultFPCA`, preserving package-selected local ranks and basis sizes in diagnostics.
- `shi-mdfpca.R` independently implements Shi et al.'s sequential alternating regression in a time–row–column tensor-product spline basis, including Gram-orthogonality constraints and score-noise correction. It is paper-faithful code, not author-supplied software.
- `spaco-minus.R` independently implements a smooth CP/PARAFAC model for missing longitudinal tensor slabs without auxiliary covariates. It is labeled SPACO-minus because neither the authors' experiment repository nor the PyPI artifact provides an importable implementation of the original estimator.
- `product-fpca.R` calls `fdarep::SparseProductFPCA` after column-major vectorization of each image and records that artificial spatial ordering in diagnostics.
- `common.R` contains data reshaping, interpolation, truncation, and package-checking helpers shared by the adapters.
- `install-dependencies.R` installs the reproducible R package layer without copying package binaries into the repository.

## Installation

Run from the repository root:

```bash
Rscript experiment/competitors/install-dependencies.R
```

The locally verified package layer on 2026-08-06 is `fdapace 0.6.0`, `face 0.1-8`, `MFPCA 1.3-11`, `mfaces 0.1-4`, `mGSFPCA 0.2.2`, and `fdarep 0.1.1`. Exact versions used in a fit are retained in adapter diagnostics when the underlying package exposes them.

## Computational structure

Approximate scaling in the number of pixel outcomes $d=pq$ is:

| method | principal fit-time or memory driver |
|---|---|
| HAPP–PACE/FACE | $d$ independent univariate smoothers followed by PCA of concatenated local scores |
| mFACEs | $d(d+1)/2$ smoothed covariance surfaces followed by a block-kernel eigendecomposition |
| mGSFPCA | $d$ generalized sparse univariate FPCAs followed by multivariate score PCA |
| Shi MD-FPCA | alternating linear solves in $b_t b_r b_c$ tensor-spline coefficients for each component |
| SPACO-minus | rank-$K$ alternating least squares, approximately linear in observed tensor entries per sweep |
| sparse Product FPCA | marginal spatial FPCA followed by sparse temporal FPCA and product-component assembly |
| spline PCA | $nd$ subject/pixel spline fits followed by PCA of $b_td$ coefficients |

Dimension studies must retain failures, safety-limit exits, and timeouts in the result table. They are evidence about computational applicability rather than missing observations to discard.
