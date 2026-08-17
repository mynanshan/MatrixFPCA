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

subject to $\beta_k^\top G\beta_k=1$ and $\beta_j^\top G\beta_k=0$ for $j<k$, where $G$ is the tensor-basis Gram matrix, $P$ combines second-derivative penalties over time, row, and column, and $w_i$ prevents subjects with more observations from dominating. After fitting, the score covariance is corrected for measurement-noise inflation and eigendecomposed to rotate the fitted span into uncorrelated FPCs. The benchmark's `auto` score option chooses a dimensionless ridge multiplier by leave-one-visit-out prediction error, refits the subject scores with that ridge, and uses the corresponding ridge-sandwich noise correction. The selected multiplier and full CV path are saved in the fit diagnostics. `shi_mdfpca` is an independent, paper-faithful implementation of the alternating-regression estimator with this explicitly documented score-stabilization extension; no runnable author implementation of the longitudinal-image method was located as of 2026-08-06, so it is not presented as author-supplied software.

### SPACO-minus smooth CP model

The implemented SPACO-style comparator uses a smooth rank-$K$ CP/PARAFAC representation without auxiliary covariates:

$$
X_i(t,r)=\mu(t,r)+\sum_{k=1}^K u_{ik}\phi_k(t)v_{rk}+\varepsilon_i(t,r).
$$

With $\Omega$ denoting observed subject-time slabs, it minimizes the masked penalized loss

$$
\sum_{(i,t)\in\Omega}\sum_{r=1}^d\left\{Y_i(t,r)-\mu(t,r)-\sum_{k=1}^Ku_{ik}\phi_k(t)v_{rk}\right\}^2+\lambda\sum_{k=1}^K\int\{\phi_k''(t)\}^2\,dt.
$$

Alternating least squares updates subject scores $u_{ik}$, pixel loadings $v_{rk}$, and B-spline coefficients for the temporal factors. Factors are normalized in the time-integral and pixel-sum inner product. After this single factor fit, a dimensionless subject-score ridge is selected by leave-one-visit-out prediction error using the fixed fitted temporal and spatial factors and the standard one-standard-error rule. This tuning requires only rank-by-rank linear solves and does not repeat the expensive factor fit. The scores are then refitted using all visits, their covariance is eigendecomposed, and the factors are rotated to produce orthogonal covariance eigenfunctions. A separate $10^{-8}$ diagonal loading is retained only for numerical stability in the alternating solves. The selected ridge, minimum-CV ridge, candidate grid, CV errors, and subject-level standard errors are stored in every SPACO result row. The authors' public experiment repository imports a `spaco` module that it does not ship, while the `spaco==1.0` PyPI artifact contains metadata but no importable implementation; `spaco_minus` is therefore explicitly an independent SPACO-style baseline rather than exact author software.

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

The primary metric is the relative squared Hilbert--Schmidt error of a common-rank covariance reconstruction,

$$
E_C=\frac{\|\widehat C-C\|_{HS}^2}{\|C\|_{HS}^2},
$$

computed for the latent covariance with the measurement-error diagonal excluded. Here $\widehat C=\sum_{k=1}^K\widehat\lambda_k\widehat\psi_k\otimes\widehat\psi_k$, where the same retained rank $K$ is used for every method, including both MatFPCA variants. Thus the comparison concerns the induced covariance operator and is unaffected by component signs, labels, or the nonidentifiable scaling of the three MatFPCA covariance factors.

Rank 9 is the primary setting. To show that conclusions are not an artifact of that single truncation, the main study repeats the otherwise identical $8\times8$ baseline at the evenly spaced ranks $K\in\{6,9,12\}$: `rank6_p8` and `rank12_p8` provide the two sensitivity runs, while `dim_p8` supplies the $K=9$ run. The lower value is 6 rather than a smaller rank because, under the baseline truth, the requested leading-five eigenspace ends inside an eigenvalue-multiplicity block and the identifiable target therefore contains six FPCs. The chosen $K$ is stored as `cov_ncomp` in every result row. This common-rank design also fixes an earlier asymmetry in which competitors were truncated but MatFPCA was scored at full rank.

Every method also receives an oracle truncation benchmark. If $C_K$ is the best rank-$K$ truncation of the true covariance, the table reports $\|C_K-C\|_{HS}^2/\|C\|_{HS}^2$ and the excess error $E_C-\|C_K-C\|_{HS}^2/\|C\|_{HS}^2$. The oracle term is especially important as dimension grows: it separates unavoidable common-rank truncation error from estimation error.

The leading-eigenspace metric is defined using the discretized version of the functional inner product. Let $t_1,\ldots,t_m$ be the common time grid and let $w_1,\ldots,w_m$ be its trapezoidal integration weights. After vectorizing an FPC over time and the $p\times q$ pixels, define

$$
W=I_{pq}\otimes\operatorname{diag}(w_1,\ldots,w_m).
$$

If $\Phi_R=[\phi_1,\ldots,\phi_R]$ and $\widehat\Phi_R=[\widehat\phi_1,\ldots,\widehat\phi_R]$ are matrices of the true and estimated leading FPC evaluations, respectively, they are normalized so that

$$
\Phi_R^\top W\Phi_R=I_R,\qquad
\widehat\Phi_R^\top W\widehat\Phi_R=I_R.
$$

This is the quadrature approximation to the usual orthonormality condition
$\sum_{a=1}^p\sum_{b=1}^q\int \phi_j(t,a,b)\phi_k(t,a,b)\,dt=\delta_{jk}$; it is what the earlier phrase "weighted-orthonormal" meant. The singular values of $\Phi_R^\top W\widehat\Phi_R$ are the cosines of the principal angles $\theta_1,\ldots,\theta_R$ between the two subspaces. The reported loss is

$$
E_{\mathrm{sub}}=\frac{1}{R}\sum_{r=1}^R\sin^2\theta_r
=1-\frac{1}{R}\|\Phi_R^\top W\widehat\Phi_R\|_F^2.
$$

The scheduled study requests $R=5$, independently of the covariance-reconstruction rank $K$. This loss is invariant to component signs and to a change of orthonormal basis within either leading subspace; in particular, it remains well-defined when individual FPCs associated with a repeated eigenvalue are not identifiable. If the requested boundary $R$ cuts through a repeated true eigenvalue, the code expands $R$ to include the complete multiplicity block. `rel_eigenvalue_sse` compares the corresponding ordered eigenvalues, and `runtime` measures fitting only.

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

### Visit-density convention

`--s_obs_min` and `--s_obs_max` are integer units of 0.05. They define the approximate fraction of the 51-point common visit grid retained independently for each subject, not the fraction of image pixels observed. At each retained visit the complete $p\times q$ image is observed. With the current rounding rule, the two study settings are:

| launcher values | nominal grid fraction | possible visits per subject | interpretation |
|---|---:|---:|---|
| `2 4` | 10–20% | 6–10 | primary sparse-longitudinal setting |
| `4 8` | 20–40% | 11–20 | less-sparse sensitivity |

Thus the baseline is genuinely sparse in longitudinal time, while the 20–40% scenario provides a less-sparse sensitivity check. The realized minimum, mean, and maximum visit counts are stored in every result row.

### Predeclared `v3` competitor settings

The rerun records every tuning choice in its CSV rows. All methods use common covariance rank 9 in the primary scenarios, with the predeclared $K=6$ and $K=12$ baseline sensitivities described above. HAPP retains three local components per pixel; a three-seed baseline sensitivity increased this cap to five but worsened mean relative HS error from 0.423 to 0.667, so the original value is retained rather than making the competitor noisier. Likewise, expanding mGSFPCA from the fixed rank-3/basis-7 fit to a four-combination AIC search slightly worsened relative HS error (0.454 to 0.474) and nearly quadrupled runtime, so the fixed package inputs are retained.

Sparse Product FPCA now permits ten spatial and five temporal marginal components instead of five and three. In the same calibration this reduced mean relative HS error from 0.304 to 0.288, with a modest leading-subspace tradeoff. SPACO-minus selects its subject-score ridge from zero and a generic quarter-decade logarithmic grid from $10^{-4}$ through $1$ by leave-one-visit-out prediction after a single factor fit. The one-standard-error rule chooses the largest candidate within one subject-level standard error of the minimum; no truth-based error enters the selection. With the final 18-candidate grid, a serial timing pair added 0.061 seconds at $8\times8$ (0.7%), while the $16\times16$ pair showed no measurable overhead relative to refitting scores at a supplied value. The tuner operates only on rank-$K$ sufficient statistics and does not repeat ALS. Shi MD-FPCA likewise uses tensor-spline roughness 0.01 and visit-CV score ridge; across five baseline calibration seeds, mean relative HS error fell from 0.541 under the old unregularized fit to 0.304. These calibration seeds are diagnostics for fixing the rerun configuration, not results to include in the final Monte Carlo table.

For a paper-quality study, vary at least sample size, image dimension, visit sparsity, noise, temporal rank, and covariance family. HAPP local truncation and spline basis/penalty settings should be fixed by a predeclared rule or included in a sensitivity analysis. A nonseparable data-generating family is important because the present product truth measures efficiency when MatFPCA is correctly specified rather than robustness to misspecification. Runtime, accuracy, timeouts, and method-level failures must be reported together.

## Compute Canada runs

`run_simu_methods.sh` creates separate Slurm arrays by computational class and image-size tier. Each task contains one method group, uses the same one-hour per-fit cutoff, and writes a scheduler-status record. Install packages under the same R module and preview manifests before submitting:

```bash
module load r/4.5.0
Rscript experiment/competitors/install-dependencies.R
experiment/run_simu_methods.sh 3 1 20 --profile pilot --dry-run
experiment/run_simu_methods.sh 100 1 40 --profile main --no-overwrite
Rscript experiment/summarize_simu_methods.R
```

The profiles are `pilot` for $8\times8$ calibration, `scaling` for baseline dimensions $8,16,32$, and `main` for dimension, sample-size, visit-density, noise, spectrum, heavy-tail, and score-dependence settings. The obsolete $4\times4$ scenario has been removed. Eligibility is now common except for the structurally quadratic mFACEs estimator:

| method group | largest scheduled square image | additional restriction |
|---|---:|---|
| all methods except mFACEs | 32 × 32 | all main scenarios |
| mFACEs | 8 × 8 | all p=8 scenarios |

The mFACEs exception is based on its computational target rather than its observed accuracy: it fits 2,080 cross-covariance surfaces at $8\times8$, 32,896 at $16\times16$, and 524,800 at $32\times32$. All other former dimension and scenario allowlists were removed after the first 100-seed run showed that the retained methods were comfortably inside the runtime budget. In particular, every spline-PCA, HAPP–PACE, and mGSFPCA fit completed at $16\times16$.

The small-image allocations were calibrated from successful pilot runs. Peak observed memory was 0.39 GB for fast tasks, 0.45 GB for medium tasks, 1.89 GB for HAPP–PACE, 2.53 GB for the combined HAPP–FACE/Shi class, and 1.34 GB for mFACEs. The small-class Slurm limits now leave 15 minutes of scheduler headroom beyond the common command cutoff:

| Slurm class | image tier | CPUs | memory | Slurm limit | command timeout |
|---|---:|---:|---:|---:|---:|
| `fast_small` | at most 8 × 8 | 4 | 2 GB | 75 min | 60 min |
| `medium_small` | at most 8 × 8 | 4 | 2 GB | 75 min | 60 min |
| `parallel_small` | at most 8 × 8 | 8 | 4 GB | 75 min | 60 min |
| `slow_small` | at most 8 × 8 | 8 | 6 GB | 75 min | 60 min |
| `mfaces` | at most 8 × 8 | 8 | 4 GB | 75 min | 60 min |

The original conservative `fast`, `medium`, `parallel`, and `slow` allocations remain in force for 16×16 and 32×32 jobs because those dimensions have not yet been resource-calibrated. CPU counts also remain unchanged: the available pilot accounting included `AllocCPUS` but not `TotalCPU`, so it does not establish whether fewer cores would preserve runtime. The launcher now prints a post-run `sacct` command containing both `TotalCPU` and `ExitCode` so the next scaling run can calibrate CPU efficiency as well as memory and time.

HAPP jobs use process workers with one BLAS thread each. Shi and mFACEs use one R process with threaded BLAS/OpenMP, avoiding nested oversubscription. A timeout records `timeout` in its `.status` file and can be retried with `--no-overwrite`, which skips only scheduler tasks marked `ok`. Method-level errors are recorded in result CSV rows and must also be checked; a scheduler task can complete even when its individual method returns an error row. The default `v3` study tag prevents results from the earlier rank-8 configuration from being skipped or overwritten accidentally.

`summarize_simu_methods.R` combines accuracy and runtime results and audits retained manifests/status files into `scheduler-status.csv`. Memory and Slurm time can be adjusted without editing the launcher through its `MATFPCA_MEM_*` and `MATFPCA_TIME_*` environment variables; the common fit cutoff is controlled by `MATFPCA_METHOD_TIMEOUT`.

For the complete paper-style comparison and a strict expected-run audit, use:

```bash
Rscript experiment/analyze_simu_methods.R
```

By default this checks the 100-seed `main` design encoded by the launcher and
writes machine-readable audit files, a complete long-form metric table, a wide
relative-Hilbert--Schmidt-error table, and Markdown tables under
`experiment/summary-methods/`. Method-level failures are counted separately
from missing result files and scheduler failures; failed fits are shown in the
completion columns but excluded from metric means. Pass `--fail_on_incomplete`
in automated checks when any missing or failed expected fit should produce a
nonzero exit status. The `pilot` and `scaling` subsets can instead be selected
with `--profile pilot` or `--profile scaling`.
