# Method-comparison simulation

`simu-methods.R` compares the proposed product-covariance estimator with
methods for sparse multivariate functions, longitudinal image tensors, and
structured sparse tensors. It uses the same generator as `simu-simple.R` and
defaults to `p = q = 4`, where the full multivariate covariance smoothers are
still practical. Implementations are isolated in `experiment/competitors/`.

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
  located as of 2026-08-06. `shi_mdfpca` is an independent paper-faithful
  implementation of its alternating-regression tensor-product-spline
  estimator; it is not represented as author-supplied software.

The implemented competitors `happ_pace` and `happ_face` follow Happ and Greven (2018),
[Multivariate Functional Principal Component Analysis for Data Observed on
Different (Dimensional) Domains](https://doi.org/10.1080/01621459.2016.1273115).
After vectorizing an image into `p*q` functional coordinates, they run sparse
univariate PACE or FACE at each pixel and eigendecompose the covariance of all
local scores. This is precisely the finite-expansion MFPCA construction. The
PACE path uses [`fdapace::FPCA`](https://cran.r-project.org/package=fdapace);
the maintained [`MFPCA`](https://cran.r-project.org/package=MFPCA) package
provides the same general framework and also documents a PACE univariate
expansion. Running both local estimators directly tests whether FACE improves
the fits enough to explain HAPP--PACE's error.

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
- `happ_face`: the same second-stage Happ--Greven score PCA with `face.sparse`
  local expansions. This is the direct PACE-versus-FACE sensitivity check.
- `mfaces`: Li, Xiao and Luo's joint multivariate sparse covariance smoother.
  It estimates every pixel auto/cross-covariance and is therefore the primary
  fully nonseparable competitor at small-to-moderate `p*q`.
- `mgsfpca`: the `mGSFPCA` sparse Gaussian multivariate FPCA implementation.
- `shi_mdfpca`: paper-faithful alternating regression on a tensor-product
  B-spline over longitudinal time, image row, and image column.
- `spaco_minus`: an independently implemented SPACO- style smooth CP model for
  missing longitudinal tensor slabs. Exact public SPACO software was not
  installable; see `competitors/README.md` for the audit trail.
- `sparse_product_fpca`: `fdarep::SparseProductFPCA` on a column-major
  vectorization of the image. It imposes a one-dimensional artificial spatial
  order and is a secondary structured baseline.
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

For explicitly low-rank competitors, `oracle_rank_rel_hs` reports the
unavoidable error of the best true covariance truncated to the fitted rank,
and `excess_rel_hs_cov` subtracts that lower bound from the observed error.
This prevents a rank-five method from looking inaccurate merely because the
reported target covariance contains more than five nonzero components.

Global FPCs are evaluated through principal angles between leading weighted
eigenspaces. `mean_sin2` and `max_sin2` are invariant to eigenfunction signs,
permutations, and rotations within a tied eigenspace. If the requested rank
cuts through a repeated true eigenvalue, the code expands it to include the
whole tie block. This matters particularly when `p = q` and the row and column
spectra coincide. `rel_eigenvalue_sse` compares the corresponding sorted global
eigenvalues. `runtime` measures fitting only, not metric calculation.

## Usage

A quick package/adapter smoke run is:

```bash
Rscript experiment/simu-methods.R \
  -n 20 -p 3 -q 3 -nt 31 -s1 3 -s2 5 \
  --eval_ncomp 3 --happ_npc 2 --spline_nbasis 5 \
  --mgs_rank_candidates 2 --mgs_basis_candidates 5 \
  --shi_ncomp 3 --spaco_rank 3 \
  --weight_folds 3 --save 0
```

The moderate-dimension default is:

```bash
Rscript experiment/simu-methods.R --seed 1
```

Select a primary subset with, for example,

```bash
Rscript experiment/simu-methods.R \
  --methods matfpca_refined,happ_pace,happ_face,mfaces,mgsfpca,shi_mdfpca \
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

The dimension study must report runtime, errors, and method failures together.
In particular, mFACEs scales with `(p*q)^2`; the default
`--mfaces_max_outcomes 25` is an intentional safety guard. Override it with
`--allow_slow_competitors 1` only through a resource-controlled small-image
job such as the launcher below.

For Compute Canada multi-seed runs, `run_simu_methods.sh` creates separate
Slurm arrays for fast, parallel, slow, and mFACEs resource classes. Every task
contains only one method group, has a method-specific hard timeout, and writes
a scheduler status record. Install the package layer under the same R module,
then preview the manifests before submitting:

```bash
module load r/4.5.0
Rscript experiment/competitors/install-dependencies.R
experiment/run_simu_methods.sh 3 1 20 --profile pilot --dry-run
experiment/run_simu_methods.sh 100 1 40 --profile main --no-overwrite
Rscript experiment/summarize_simu_methods.R
```

The `main` profile varies image dimension (`4, 8, 16, 32`), sample size,
visit density, noise, spatial/temporal spectrum, heavy-tailed scores, and score
dependence. Eligibility is case-specific:

| method group | largest scheduled square image | additional restriction |
|---|---:|---|
| MatFPCA, SPACO-, sparse Product FPCA, Shi MD-FPCA | 32 x 32 | all main scenarios |
| spline PCA, HAPP--PACE, mGSFPCA | 16 x 16 | all eligible scenarios |
| HAPP--FACE | 8 x 8 | baseline, denser visits, and noise sensitivity |
| mFACEs | 8 x 8 | baseline dimensions and one denser-visit setting |

HAPP jobs use eight process workers with one BLAS thread each. Shi and mFACEs
use one R process with eight BLAS/OpenMP threads, avoiding nested
oversubscription. The command-level limits range from 50 minutes for fast
methods to 7.5 hours for mFACEs, below their corresponding Slurm limits. A
timed-out method records `timeout` in its `.status` file and can be resubmitted
with `--no-overwrite`, which skips only tasks whose status is `ok`.
`summarize_simu_methods.R` additionally combines the retained manifests and
status files into `scheduler-status.csv`, so timeouts and pending tasks remain
visible rather than silently disappearing from the accuracy table.
After the pilot, memory, Slurm time, and command timeouts can be adjusted
without editing the launcher through `MATFPCA_MEM_*`, `MATFPCA_TIME_*`, and
`MATFPCA_TIMEOUT_*` environment variables recognized by the launcher.
