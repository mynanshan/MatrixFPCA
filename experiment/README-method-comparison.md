# Method-comparison simulation

`simu-methods.R` compares the proposed product-covariance estimator with methods for sparse multivariate functions, longitudinal image tensors, and structured sparse tensors. It uses the generator from `simu-simple.R`; competitor implementations and package adapters are isolated in `experiment/competitors/`.

## Statistical problem and notation

For subject $i=1,\ldots,n$, an image is observed only at subject-specific times $t_{i\ell}$:

$$
Y_i(t_{i\ell},a,b)=X_i(t_{i\ell},a,b)+\varepsilon_i(t_{i\ell},a,b), \qquad a=1,\ldots,p,\quad b=1,\ldots,q.
$$

Write $r=(a,b)$ for one of the $d=pq$ pixels. The latent covariance is the operator with kernel

$$
C_{rr'}(s,t)=\operatorname{Cov}\{X_i(s,r),X_i(t,r')\}.
$$

The proposed MatFPCA model assumes a separable three-way kernel,

$$
C_{(a,b),(a',b')}(s,t)=K(s,t)U_{aa'}V_{bb'}=\{K\otimes U\otimes V\}_{(s,a,b),(t,a',b')},
$$

where the individual scales of $K,U,V$ are not identifiable because $c_1K,c_2U,(c_1c_2)^{-1}V$ define the same full covariance. `matfpca_refined` estimates the factors with cross-fitted Gaussian-refined weights; `matfpca_identity` is the identity-weight ablation. All external competitors below estimate a nonseparable or differently structured model and are evaluated through the induced full covariance, not through arbitrarily scaled factors.

## Competitor methodology

### Happ–Greven multivariate FPCA with PACE or FACE

Happ and Greven's finite-expansion construction first performs a separate univariate FPCA for every pixel. For pixel $r$,

$$
X_i^{(r)}(t)-\mu_r(t) \approx \sum_{m=1}^{M_r}\xi_{irm}\phi_{rm}(t).
$$

The local scores are concatenated into $\boldsymbol\xi_i=(\xi_{i11},\ldots,\xi_{idM_d})^\top$. If $c_k$ is the $k$-th eigenvector of the empirical covariance of these concatenated scores, the multivariate FPC and eigenvalue are

$$
\psi_k^{(r)}(t)=\sum_{m=1}^{M_r}c_{k,rm}\phi_{rm}(t), \qquad \nu_k=\operatorname{Var}(c_k^\top\boldsymbol\xi_i).
$$

Thus cross-pixel dependence is estimated only in the second-stage score covariance; the first-stage temporal bases are pixel-specific. `happ_pace` estimates each local mean and covariance from irregular observations using `fdapace::FPCA`, obtains conditional-expectation scores, retains `--happ_npc` local components, and then performs the score PCA above. `happ_face` uses the same second stage but replaces PACE by `face::face.sparse`, whose covariance smoothing uses a fast penalized-spline representation. Comparing the two isolates the effect of the univariate smoother while holding the multivariate construction fixed. The implementation is the Happ–Greven finite-expansion algorithm adapted to vectorized image pixels; no spatial product structure is imposed.

### mFACEs

mFACEs directly estimates every auto- and cross-covariance surface. After estimating coordinate means, its raw cross-products are

$$
G_{i,rr'}(t_{i\ell},t_{im})=\{Y_i(t_{i\ell},r)-\widehat\mu_r(t_{i\ell})\}\{Y_i(t_{im},r')-\widehat\mu_{r'}(t_{im})\}.
$$

For each $r\le r'$, a tensor-product spline surface $B(s)^\top\Theta_{rr'}B(t)$ is fitted to these cross-products with roughness regularization, schematically

$$
\widehat\Theta_{rr'}=\arg\min_\Theta \sum_{i,\ell,m}\left\{G_{i,rr'}(t_{i\ell},t_{im})-B(t_{i\ell})^\top\Theta B(t_{im})\right\}^2+\lambda_{rr'}\mathcal J(\Theta).
$$

The fitted surfaces are assembled into the block kernel $\widehat C_{rr'}(s,t)$, symmetrized, and jointly eigendecomposed. Measurement-error variance is estimated from the auto-covariance diagonals. This is the primary fully nonseparable covariance-smoothing baseline, but it requires $d(d+1)/2$ bivariate smoothers and consequently becomes expensive as $d=pq$ grows. The adapter calls `mfaces::mface.sparse` and intentionally restricts routine experiments to small images.

### mGSFPCA

mGSFPCA begins with a generalized sparse FPCA for every outcome. In the Gaussian case used here, the coordinate model reduces to

$$
X_i^{(r)}(t)=\mu_r(t)+\sum_{m=1}^{M_r}\xi_{irm}\phi_{rm}(t)+\varepsilon_i^{(r)}(t),
$$

with B-spline representations for $\mu_r$ and $\phi_{rm}$. The package estimates the local rank and basis size from supplied candidate sets, forms subject-level local scores, estimates their full cross-coordinate covariance, and eigendecomposes it to obtain global multivariate FPCs, just as in a finite basis representation. In non-Gaussian applications the same framework models an exponential-family conditional mean through a link function; the benchmark deliberately uses its Gaussian path so it targets the current continuous-image model. The adapter calls `mGSFPCA::spMultFPCA` and converts `fullEigFun` and `values` into the common covariance-eigenfit representation.

### Shi et al. multi-dimensional FPCA

Shi et al.'s method represents each longitudinal image component in a tensor-product spline basis. Let

$$
b(t,a,b)=B_t(t)\otimes B_r(a)\otimes B_c(b), \qquad \psi_k(t,a,b)=b(t,a,b)^\top\beta_k.
$$

For a subject-specific observation vector $y_i$ and its tensor-basis design $X_i$, the rank-$K$ model is

$$
y_i\approx X_i\sum_{k=1}^K\alpha_{ik}\beta_k.
$$

Components are fitted sequentially by alternating between least-squares score updates and a penalized coefficient update,

$$
\min_{\beta_k,\{\alpha_{ik}\}}\sum_i w_i\left\|r_{ik}-\alpha_{ik}X_i\beta_k\right\|_2^2+\rho\,\beta_k^\top P\beta_k,
$$

subject to $\beta_k^\top G\beta_k=1$ and $\beta_j^\top G\beta_k=0$ for $j<k$, where $G$ is the tensor-basis Gram matrix, $P$ combines second-derivative penalties over time, row, and column, and $w_i$ prevents subjects with more observations from dominating. After fitting, the score covariance is corrected for least-squares measurement-noise inflation and eigendecomposed to rotate the fitted span into uncorrelated FPCs. `shi_mdfpca` is an independent, paper-faithful implementation of this alternating-regression estimator; no runnable author implementation of the longitudinal-image method was located as of 2026-08-06, so it is not presented as author-supplied software.

### SPACO-minus smooth CP model

The implemented SPACO-style comparator uses a smooth rank-$K$ CP/PARAFAC representation without auxiliary covariates:

$$
X_i(t,r)=\mu(t,r)+\sum_{k=1}^K u_{ik}\phi_k(t)v_{rk}+\varepsilon_i(t,r).
$$

With $\Omega$ denoting observed subject-time slabs, it minimizes the masked penalized loss

$$
\sum_{(i,t)\in\Omega}\sum_{r=1}^d\left\{Y_i(t,r)-\mu(t,r)-\sum_{k=1}^Ku_{ik}\phi_k(t)v_{rk}\right\}^2+\lambda\sum_{k=1}^K\int\{\phi_k''(t)\}^2\,dt.
$$

Alternating least squares updates subject scores $u_{ik}$, pixel loadings $v_{rk}$, and B-spline coefficients for the temporal factors. Factors are normalized in the time-integral and pixel-sum inner product. Finally, the covariance of $u_i$ is eigendecomposed and the factors are rotated to produce orthogonal covariance eigenfunctions. The authors' public experiment repository imports a `spaco` module that it does not ship, while the `spaco==1.0` PyPI artifact contains metadata but no importable implementation; `spaco_minus` is therefore explicitly an independent SPACO-style baseline rather than exact author software.

### Sparse Product FPCA

`fdarep::SparseProductFPCA` constructs marginal spatial and temporal FPCAs and combines them into product components. After converting each image to the artificial one-dimensional coordinate $r=1,\ldots,pq$, the fitted eigenfunctions have the form

$$
\psi_{jk}(t,r)=\phi_k(t)\eta_j(r),
$$

with component variances estimated from the corresponding sparse product scores. This avoids a full $pq$-variate covariance smoother and is computationally attractive, but the spatial eigenfunctions $\eta_j(r)$ depend on the chosen column-major pixel ordering and do not distinguish the two image axes. It is therefore a structured secondary baseline, not an order-invariant matrix method.

### Smooth-then-PCA spline baseline

For each subject $i$ and pixel $r$, `spline_pca` first estimates a complete trajectory $B(t)^\top c_{ir}$ through penalized least squares,

$$
\widehat c_{ir}=\arg\min_c\sum_{\ell}\{Y_i(t_{i\ell},r)-B(t_{i\ell})^\top c\}^2+\lambda c^\top Rc,
$$

where $R$ is the integrated squared-second-derivative penalty. The coefficient vectors are whitened by the basis Gram matrix $G=\int B(t)B(t)^\top dt$, concatenated across pixels, and subjected to ordinary PCA. This yields the correct functional inner product but asks a very sparse subject-level smoother to reconstruct each complete trajectory before pooling information, so it is deliberately treated as a weak baseline.

## Relationship to other literature

Shi, Jiang, Ma, Beg and Cao (2024), [Dynamic Survival Prediction Using Sparse Longitudinal Images via Multi-Dimensional Functional Principal Component Analysis](https://doi.org/10.1080/10618600.2024.2335182), is the closest published formulation because it directly studies images observed sparsely and irregularly over longitudinal time.

The implemented `happ_pace` and `happ_face` follow Happ and Greven (2018), [Multivariate Functional Principal Component Analysis for Data Observed on Different (Dimensional) Domains](https://doi.org/10.1080/01621459.2016.1273115). The maintained [`MFPCA`](https://cran.r-project.org/package=MFPCA) package provides the same general finite-expansion framework; the benchmark uses explicit PACE and FACE adapters so their local smoothing behavior and failures remain observable.

Chen and Jiang (2017), [Multi-dimensional Functional Principal Component Analysis](https://doi.org/10.1007/s11222-016-9679-5), handles general multidimensional domains with local-linear smoothing, FFT/GPU computation, and randomized eigendecomposition. It is relevant, but its generic scattered design and implementation are not naturally matched to complete images observed on sparse time slabs.

Zipunnikov et al.'s [HD-LFPCA](https://www.biostat.jhsph.edu/~vzipunni/vadim_zipunnikov_hd-lfpca.pdf) is scalable for longitudinal high-dimensional images, but its base process is a random-intercept/random-slope plus visit-level decomposition rather than the unrestricted smooth temporal covariance used here. Hu and Yao (2022), [Sparse Functional PCA in High Dimensions](https://doi.org/10.5705/ss.202020.0445), assumes sparse eigenfunctions across a large number of functional coordinates and would be especially relevant in a future simulation with spatially sparse FPCs; that assumption is poorly aligned with the current dense AR(1)/Fourier truth.

## Metrics and identifiability

The primary metric is relative squared Hilbert–Schmidt covariance error,

$$
E_C=\frac{\|\widehat C-C\|_{HS}^2}{\|C\|_{HS}^2},
$$

computed for the latent covariance with the measurement-error diagonal excluded. Comparing the full operator makes factor scale and order irrelevant.

Explicitly low-rank competitors also receive an oracle truncation benchmark. If $C_K$ is the best rank-$K$ truncation of the true covariance, the table reports $\|C_K-C\|_{HS}^2/\|C\|_{HS}^2$ and the excess error $E_C-\|C_K-C\|_{HS}^2/\|C\|_{HS}^2$. This prevents a rank-eight fit from being penalized solely because the target contains additional nonzero components.

For leading eigenspaces, let $P$ and $\widehat P$ contain weighted-orthonormal true and estimated FPCs. The reported mean squared sine of the principal angles is

$$
E_{\mathrm{sub}}=\frac{1}{K}\sum_{k=1}^K\sin^2\theta_k=1-\frac{1}{K}\|P^\top W\widehat P\|_F^2.
$$

It is invariant to signs, permutations, and rotations inside tied eigenspaces. If the requested evaluation boundary cuts through a repeated true eigenvalue, the code expands the boundary to include the complete tie block. `rel_eigenvalue_sse` compares the corresponding sorted eigenvalues, and `runtime` measures fitting only.

## Temporary local illustration: three seeds

The following small batch was run locally on 2026-08-06 as an implementation check. It is deliberately stored here as a temporary illustration rather than treated as a simulation conclusion: there are only three seeds, one small dimension, and product-structured truth favoring MatFPCA.

The setting was $n=50$, $p=q=4$, 51 time-grid points, 10–20% of visits observed per subject, noise standard deviation 0.1, AR(1) spatial factors, temporal rank 5, independent Gaussian scores, competitor covariance rank 8, five evaluated leading components, three local HAPP components per pixel, four method workers, and one BLAS thread. Every one of the 30 method fits completed successfully. Values are means across seeds 1–3; parentheses contain Monte Carlo standard deviations for covariance error.

| method | relative HS error | excess relative HS error | mean $\sin^2$ angle | runtime (s) | successful seeds |
|---|---:|---:|---:|---:|---:|
| MatFPCA refined | 0.0457 (0.0525) | 0.0457 | 0.1718 | 3.100 | 3/3 |
| MatFPCA identity | 0.0632 (0.0639) | 0.0632 | 0.1743 | 0.612 | 3/3 |
| sparse Product FPCA | 0.1466 (0.0742) | 0.1189 | 0.1328 | 1.357 | 3/3 |
| HAPP–FACE | 0.1945 (0.0353) | 0.1667 | 0.1868 | 7.030 | 3/3 |
| mGSFPCA | 0.2189 (0.0577) | 0.1911 | 0.2104 | 5.028 | 3/3 |
| HAPP–PACE | 0.2249 (0.0479) | 0.1971 | 0.2185 | 1.638 | 3/3 |
| mFACEs | 0.2268 (0.0772) | 0.1991 | 0.2478 | 40.075 | 3/3 |
| SPACO-minus | 0.2423 (0.0859) | 0.2145 | 0.2260 | 4.342 | 3/3 |
| Shi MD-FPCA | 0.3000 (0.0758) | 0.2723 | 0.3076 | 13.524 | 3/3 |
| spline PCA | 0.4090 (0.2312) | 0.3812 | 0.3008 | 0.034 | 3/3 |

The exact command pattern was:

```bash
for seed in 1 2 3; do
  Rscript experiment/simu-methods.R \
    -n 50 -p 4 -q 4 -nt 51 -s1 2 -s2 4 -sig 0.1 \
    --rank_time 5 --eigset 2 --eigvalset 1 --eigensystem ar1 \
    --scores_distr normal --score_dependence independent \
    --eval_ncomp 5 --cov_ncomp 8 --happ_npc 3 \
    --mfaces_max_outcomes 100 --allow_slow_competitors 1 \
    --shi_ncomp 8 --spaco_rank 8 --ncores 4 --blas_threads 1 \
    --backend multicore --weight_folds 5 --seed "$seed" \
    --results_dir experiment/results-methods/local-illustration \
    --result_tag local_illustration_p4
done
```

The ignored local CSV files can be regenerated with that command. The numerical table above is the durable snapshot.

## General usage

A quick adapter smoke run is:

```bash
Rscript experiment/simu-methods.R \
  -n 20 -p 3 -q 3 -nt 31 -s1 3 -s2 5 \
  --eval_ncomp 3 --happ_npc 2 --spline_nbasis 5 \
  --mgs_rank_candidates 2 --mgs_basis_candidates 5 \
  --shi_ncomp 3 --spaco_rank 3 \
  --weight_folds 3 --save 0
```

Compact CSV files are written to `experiment/results-methods/`. Add `--save_fits 1` only for diagnostic runs because fitted eigenfunction arrays can be large.

For a paper-quality study, vary at least sample size, image dimension, visit sparsity, noise, temporal rank, and covariance family. HAPP local truncation and spline basis/penalty settings should be fixed by a predeclared rule or included in a sensitivity analysis. A nonseparable data-generating family is important because the present product truth measures efficiency when MatFPCA is correctly specified rather than robustness to misspecification. Runtime, accuracy, timeouts, and method-level failures must be reported together.

## Compute Canada runs

`run_simu_methods.sh` creates separate Slurm arrays by computational class and image-size tier. Each task contains one method group, has a method-specific hard timeout, and writes a scheduler-status record. Install packages under the same R module and preview manifests before submitting:

```bash
module load r/4.5.0
Rscript experiment/competitors/install-dependencies.R
experiment/run_simu_methods.sh 3 1 20 --profile pilot --dry-run
experiment/run_simu_methods.sh 100 1 40 --profile main --no-overwrite
Rscript experiment/summarize_simu_methods.R
```

The profiles are `pilot` for $4\times4$ and $8\times8$ calibration, `scaling` for baseline dimensions $4,8,16,32$, and `main` for dimension, sample-size, visit-density, noise, spectrum, heavy-tail, and score-dependence settings. Eligibility is method-specific:

| method group | largest scheduled square image | additional restriction |
|---|---:|---|
| MatFPCA, SPACO-minus, sparse Product FPCA, Shi MD-FPCA | 32 × 32 | all main scenarios |
| spline PCA, HAPP–PACE, mGSFPCA | 16 × 16 | all eligible scenarios |
| HAPP–FACE | 8 × 8 | baseline, denser visits, and noise sensitivity |
| mFACEs | 8 × 8 | baseline dimensions and one denser-visit setting |

The small-image allocations were calibrated from three successful pilot seeds at both $4\times4$ and $8\times8$. Peak observed memory was 0.39 GB for fast tasks, 0.45 GB for medium tasks, 1.89 GB for HAPP–PACE, 2.53 GB for the combined HAPP–FACE/Shi class, and 1.34 GB for mFACEs. The longest pilot wall time was 4 minutes 12 seconds for 8×8 mFACEs. The resulting defaults retain additional headroom:

| Slurm class | image tier | CPUs | memory | Slurm limit | command timeout |
|---|---:|---:|---:|---:|---:|
| `fast_small` | at most 8 × 8 | 4 | 2 GB | 15 min | 10 min |
| `medium_small` | at most 8 × 8 | 4 | 2 GB | 30 min | 20 min |
| `parallel_small` | at most 8 × 8 | 8 | 4 GB | 30 min | 20 min |
| `slow_small` | at most 8 × 8 | 8 | 6 GB | 30 min | 20 min |
| `mfaces` | at most 8 × 8 | 8 | 4 GB | 30 min | 20 min |

The original conservative `fast`, `medium`, `parallel`, and `slow` allocations remain in force for 16×16 and 32×32 jobs because those dimensions have not yet been resource-calibrated. CPU counts also remain unchanged: the available pilot accounting included `AllocCPUS` but not `TotalCPU`, so it does not establish whether fewer cores would preserve runtime. The launcher now prints a post-run `sacct` command containing both `TotalCPU` and `ExitCode` so the next scaling run can calibrate CPU efficiency as well as memory and time.

HAPP jobs use process workers with one BLAS thread each. Shi and mFACEs use one R process with threaded BLAS/OpenMP, avoiding nested oversubscription. Command-level timeouts are shorter than their corresponding Slurm limits. A timeout records `timeout` in its `.status` file and can be retried with `--no-overwrite`, which skips only scheduler tasks marked `ok`. Method-level errors are recorded in result CSV rows and must also be checked; a scheduler task can complete even when its individual method returns an error row.

`summarize_simu_methods.R` combines accuracy and runtime results and audits retained manifests/status files into `scheduler-status.csv`. Memory, Slurm time, and command timeouts can be adjusted without editing the launcher through its `MATFPCA_MEM_*`, `MATFPCA_TIME_*`, and `MATFPCA_TIMEOUT_*` environment variables.
