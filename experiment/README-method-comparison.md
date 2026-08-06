# Method-comparison simulation

`simu-methods.R` compares the proposed product-covariance estimator with
methods that treat a longitudinal image as a vector-valued function. It uses
the same generator as `simu-simple.R`, but defaults to `p = q = 8` because the
nonseparable competitors are intentionally much less scalable.

## Literature conclusion

There is one especially close published method:

- Shi, Jiang, Ma, Beg and Cao (2024), [Dynamic Survival Prediction Using
  Sparse Longitudinal Images via Multi-Dimensional Functional Principal
  Component Analysis](https://doi.org/10.1080/10618600.2024.2335182), studies
  images collected sparsely and irregularly over longitudinal time. Its
  multi-dimensional FPCA uses alternating regression with tensor-product
  B-splines to avoid a full covariance eigendecomposition. This is a direct
  methodological competitor and should be discussed in the paper. The article
  page exposes a methodological supplement, but no runnable public code was
  located as of 2026-08-06. It is therefore not silently approximated or
  attributed to the authors in the current benchmark.

The implemented competitor `happ_pace` follows Happ and Greven (2018),
[Multivariate Functional Principal Component Analysis for Data Observed on
Different (Dimensional) Domains](https://doi.org/10.1080/01621459.2016.1273115).
After vectorizing an image into `p*q` functional coordinates, it runs sparse
univariate PACE at each pixel and eigendecomposes the covariance of all local
scores. This is precisely the finite-expansion MFPCA construction. The current
implementation uses [`fdapace::FPCA`](https://cran.r-project.org/package=fdapace)
for the local fits; the maintained [`MFPCA`](https://cran.r-project.org/package=MFPCA)
package provides the same general framework and also documents a PACE
univariate expansion.

Other related methods are less direct:

- Chen and Jiang (2017), [Multi-dimensional Functional Principal Component
  Analysis](https://doi.org/10.1007/s11222-016-9679-5), handles general
  multidimensional domains with local-linear smoothing, FFT/GPU computation,
  and randomized eigendecomposition. It is relevant, but its generic scattered
  design and implementation are not naturally matched to complete images
  observed on sparse time slabs.
- Zipunnikov et al.'s
  [HD-LFPCA](https://www.biostat.jhsph.edu/~vzipunni/vadim_zipunnikov_hd-lfpca.pdf)
  is scalable for longitudinal high-dimensional images, but its base process is
  a random-intercept/random-slope plus visit-level decomposition, rather than
  the unrestricted smooth temporal covariance used here.
- Hu and Yao (2022), [Sparse Functional PCA in High
  Dimensions](https://doi.org/10.5705/ss.202020.0445), targets a large number of
  functional coordinates by assuming sparse eigenfunctions. It becomes a
  compelling comparator if a future simulation includes spatially sparse FPCs;
  it is not well aligned with the current dense AR(1)/Fourier truth.

## Implemented methods

- `matfpca_refined`: proposed MatFPCA with cross-fitted Gaussian-refined
  weights.
- `matfpca_identity`: the identity-weight ablation retained as a bridge to the
  existing simulations.
- `happ_pace`: nonseparable Happ--Greven MFPCA after vectorizing pixels, adapted
  from `/home/muye/Downloads/matfpca-happ.R` and made self-contained.
- `spline_pca`: a fast, transparent nonseparable baseline. It ridge-smooths
  every subject/pixel trajectory in a common B-spline basis and performs PCA in
  the correct basis-Gram inner product. It is a baseline, not an implementation
  of Shi et al. (2024).

## Metrics and identifiability

The primary covariance metric is relative squared Hilbert--Schmidt error,

`||C_hat - C||_HS^2 / ||C||_HS^2`,

computed for the latent covariance (the measurement-error diagonal is
excluded). Comparing the full operator makes the arbitrary factor rescaling
`a U`, `b V`, `(ab)^(-1) K` irrelevant.

Global FPCs are evaluated through principal angles between leading weighted
eigenspaces. `mean_sin2` and `max_sin2` are invariant to eigenfunction signs,
permutations, and rotations within a tied eigenspace. If the requested rank
cuts through a repeated true eigenvalue, the code expands it to include the
whole tie block. This matters particularly when `p = q` and the row and column
spectra coincide. `rel_eigenvalue_sse` compares the corresponding sorted global
eigenvalues. `runtime` measures fitting only, not metric calculation.

## Usage

A quick smoke run is:

```bash
Rscript experiment/simu-methods.R \
  -n 20 -p 3 -q 3 -nt 31 -s1 3 -s2 5 \
  --eval_ncomp 3 --happ_npc 2 --spline_nbasis 5 \
  --weight_folds 3 --save 0
```

The moderate-dimension default is:

```bash
Rscript experiment/simu-methods.R --seed 1
```

Select a subset with, for example,

```bash
Rscript experiment/simu-methods.R \
  --methods matfpca_refined,happ_pace,spline_pca \
  --ncores 4 --backend multicore --seed 1
```

Compact CSV files are written to `experiment/results-methods/`. Add
`--save_fits 1` only for diagnostic runs because the fitted eigenfunction
arrays can be large.

For a paper-quality study, vary at least sample size, `p*q`, visit sparsity,
noise, and temporal rank. The HAPP local truncation (`--happ_npc`) and spline
basis/penalty (`--spline_nbasis`, `--spline_smoothness`) should be fixed by a
predeclared rule or studied in a sensitivity grid. A second data-generating
family with nonseparable covariance is also important: the present
product-structured truth measures the efficiency gained when MatFPCA's model is
correct, not robustness to misspecification.
