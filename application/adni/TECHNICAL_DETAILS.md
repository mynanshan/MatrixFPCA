# Technical details of the ADNI MatFPCA analysis

This note describes the analysis implemented in
[`adni_matfpca.R`](adni_matfpca.R). It distinguishes three quantities that are
easy to conflate:

1. variance explained by the MatFPCA eigencomponents;
2. variation in the selected subject-score cloud explained by the morphology
   curve; and
3. held-out image reconstruction accuracy.

In particular, the reported 48.0% belongs to item 2. It is not the fraction of
total MRI intensity variation or temporal variation explained by one component.

## 1. Observations and preprocessing

For subject $i$, visit $v$, row $r$, and column $c$, let

\[
Y_{iv}(r,c) = Y_i(r,c,t_{iv})
\]

denote the axial MRI slice observed at time $t_{iv}$ years from baseline.
The current analysis has 393 subjects and 1,504 images through four years.

The application performs the following operations before covariance fitting:

- Images are block-averaged to $114\times96$ pixels.
- Screening visits are retained, baseline visits are assigned exactly $t=0$,
  and other times are snapped at one-day resolution.
- Subjects with fewer than two usable images are removed.
- A consensus brain mask is estimated using training images only. A pixel must
  be foreground in at least half of the training images; the largest connected
  component is retained and internal holes are filled.
- Twenty percent of subjects are selected reproducibly, and one visit from
  each selected subject is held out. The subject's remaining visits are still
  available for estimating that subject's scores. This is therefore a
  held-out-visit prediction experiment, not prediction for entirely new
  subjects.

No diagnosis or clinical outcome label enters sampling, model fitting,
continuum estimation, or plotting.

## 2. Time-varying mean

The model separates the population mean from subject-specific covariance:

\[
Y_i(r,c,t)=\mu(r,c,t)+X_i(r,c,t)+\varepsilon_i(r,c,t).
\]

For every pixel, $\mu(r,c,t)$ is represented by a natural-spline regression in
time with five total design columns. If subject $i$ contributes $n_i$
training visits, each of those visits receives weight $1/n_i$. Thus a subject
with many scans does not dominate the population mean. The fit uses training
images only. MatFPCA is applied to

\[
\widetilde Y_i(r,c,t)=Y_i(r,c,t)-\widehat\mu(r,c,t).
\]

The fitted mean is added back to every reconstruction.

## 3. Fully separable covariance model

The centered process is assumed to have a separable covariance:

\[
\operatorname{Cov}\{X_i(r,c,s),X_i(r',c',t)\}
=U(r,r')V(c,c')K(s,t).
\]

Let the marginal eigensystems be

\[
U\phi_a=\lambda_a^{(r)}\phi_a,\qquad
V\psi_b=\lambda_b^{(c)}\psi_b,
\]

and

\[
\int K(s,t)\theta_k(t)\,dt
=\lambda_k^{(t)}\theta_k(s).
\]

The corresponding product expansion is

\[
X_i(r,c,t)
=\sum_{a,b,k}\xi_{i,abk}\,
  \phi_a(r)\psi_b(c)\theta_k(t),
\]

with

\[
\operatorname{Var}(\xi_{i,abk})
=\lambda_a^{(r)}\lambda_b^{(c)}\lambda_k^{(t)}.
\]

Individual rescalings of $U,V,K$ are not identifiable, but the product
covariance, product eigenfunctions, product eigenvalues, and reconstructions
are identifiable. The plotted marginal eigenfunctions use the normalization
chosen by their eigendecompositions.

### Covariance estimation

The default `gaussian_refined` estimator uses three subject-level folds.
Identity-weight fits on the other folds provide pilot covariance and noise
estimates, from which Gaussian inverse-variance mode weights are constructed.
The held-out fold then contributes with weights learned without using its own
outcomes. The final temporal covariance uses the two-step estimator: smoothed
off-diagonal visit products estimate the latent temporal surface, while the
diagonal supplies the measurement-error correction. The fitted error variance
in the current run is $0.02779$.

The cross-fitting here concerns covariance weights. It is separate from the
one-visit reconstruction holdout described above.

## 4. Truncation, FVE, and subject scores

The row, column, and time ranks are the smallest marginal ranks reaching 95%
positive-eigenvalue cumulative variance. The current retained ranks are

\[
(I,J,K)=(34,24,1).
\]

Product components are ordered by

\[
\lambda_{abk}
=\lambda_a^{(r)}\lambda_b^{(c)}\lambda_k^{(t)}.
\]

Separate marginal fractions of variance explained are available because of
separability. For example,

\[
\operatorname{FVE}_{t,k}
=\frac{\lambda_k^{(t)}}{\sum_h\lambda_h^{(t)}},
\qquad
\operatorname{FVE}_{s,ab}
=\frac{\lambda_a^{(r)}\lambda_b^{(c)}}
{\sum_{p,q}\lambda_p^{(r)}\lambda_q^{(c)}}.
\]

The temporal diagnostic uses trapezoidal quadrature on the irregular compact
time grid. It finds 99.017% temporal FVE for the nearly constant first mode and
0.983% for the strongly non-constant second mode. The automatic 95% rule keeps
only the first mode. Consequently, the fitted scores used below describe
predominantly stable morphology; the morphology coordinate is not a coordinate
of temporal change.

For a subject observed at irregular times, the retained score vector is the
Gaussian conditional-mean or BLUP score. In generic matrix notation, if
$x_i$ stacks that subject's centered observed pixels, $B_i$ evaluates the
retained product eigenfunctions at the observed visits, and
$\Lambda=\operatorname{diag}(\lambda_{abk})$, then

\[
\widehat\xi_i
=\Lambda B_i^\top
\left(B_i\Lambda B_i^\top+\widehat\sigma^2 I\right)^{-1}x_i.
\]

The implementation exploits the separable structure rather than constructing
this full matrix literally. Scores are then reordered by decreasing product
eigenvalue.

## 5. The continuous morphology coordinate

### 5.1 Input score cloud and the role of $D$

Let

\[
z_i=(\widehat\xi_{i1},\ldots,\widehat\xi_{iD})^\top,
\qquad z_i\in\mathbb R^D,
\]

contain the leading non-degenerate product-FPC scores for subject $i$. The
dimension is controlled by `--score-pcs=D`; the default is $D=30$. Every
continuum-dependent output has a suffix such as `_pcs30`, so fits using several
values of $D$ can coexist. In the default fit, these 30 components account for
18.86% of the variance among all retained product components.

The scores are not standardized before fitting the curve. Hence
Euclidean distance gives more influence to high-variance components, which is
appropriate when the goal is to represent variation in the reconstructed MRI
process. Standardizing would give weak FPCs the same influence as dominant
ones.

### 5.2 What the principal curve represents

The morphology curve is a smooth map

\[
g:[0,1]\longrightarrow\mathbb R^D.
\]

Its fitting algorithm is implemented in
[`adni_continuum_utils.R`](adni_continuum_utils.R).

It helps to separate three objects:

- $z_i\in\mathbb R^D$ is the observed estimated score vector for subject
  $i$, represented by one point in the score cloud.
- $g(u)\in\mathbb R^D$ is a point on a one-dimensional smooth curve. Although
  $g$ has $D$ output coordinates, it has only one input coordinate $u$.
- $p_i=\Pi_g(z_i)$ is the point on the fitted curve closest to $z_i$.
  Subject $i$'s morphology coordinate is the normalized distance travelled
  along the curve from its starting endpoint to $p_i$.

Thus the method does not fit 30 unrelated coordinates. It fits 30 coordinate
functions $g_1(u),\ldots,g_D(u)$ that share the same scalar argument $u$;
together they trace one path through the $D$-dimensional cloud. Nearby values
of $u$ correspond to nearby fitted combinations of all selected FPC scores.

This is a nonlinear analogue of the first PCA axis. PCA projects each point
onto one straight line. The principal curve replaces that line by a smooth
path, while retaining a single scalar coordinate. Unlike UMAP or t-SNE, this
curve is fitted in the original score space, supplies an explicit projection
and residual for every subject, and can be decoded back into an MRI.

The implementation uses an alternating algorithm because neither the curve nor
the subject coordinates are known initially:

1. guess a one-dimensional ordering of subjects;
2. smooth every score as a function of that ordering to obtain a curve;
3. project every subject onto the new curve and replace their coordinate by
   normalized curve arc length; and
4. repeat steps 2--3 until the projection loss stabilizes.

The following sections give the exact calculations.

### 5.3 Initialization

The score matrix is column-centered, and ordinary PCA supplies an initial
one-dimensional ordering. The first-PC scores are linearly rescaled to
$[0,1]$. This linear initialization by itself explains 46.08% of empirical
score variation.

### 5.4 Alternating spline and projection steps

Given current coordinates $u_i\in[0,1]$, every score dimension is smoothed
separately:

\[
g_d(u)=\operatorname{SmoothSpline}
\{(u_i,z_{id}):i=1,\ldots,n\},
\qquad d=1,\ldots,D.
\]

Each smoother has fixed effective degrees of freedom 5. The $D$ smoothers are
evaluated on a dense common grid; with 393 subjects the grid has
$\max(800,4n)=1{,}572$ points. Consecutive fitted points define a polyline
approximation to $g$.

For subject $i$, the algorithm finds the closest point on every polyline
segment. For a segment from $g_m$ to $g_{m+1}$, write
$d_m=g_{m+1}-g_m$. Its candidate projection is

\[
p_{im}=g_m+\alpha_{im}d_m,
\qquad
\alpha_{im}=\left[
\frac{(z_i-g_m)^\top d_m}{d_m^\top d_m}
\right]_{[0,1]},
\]

where the brackets truncate the scalar to $[0,1]$. The segment giving the
smallest $\lVert z_i-p_{im}\rVert^2$ supplies the subject's projection
$p_i$. The new coordinate is cumulative arc length to $p_i$, divided by the
total curve length.

Spline smoothing and projection alternate until the relative change in total
squared projection loss is below $10^{-4}$, or 20 iterations have been
performed.

In pseudocode, [`fit_score_principal_curve()`](adni_continuum_utils.R) performs:

```text
Z = matrix of n subjects by D raw scores
u = rescaled first-PCA score of each subject

repeat:
    for d = 1,...,D:
        fit a df-5 smoothing spline to (u_i, Z_id)
        evaluate it on a common dense u-grid
    join the grid evaluations into a polyline g
    for each subject i:
        project Z_i onto every segment of g
        retain the closest projection p_i
        u_i = arc length from the curve start to p_i / total arc length
    loss = sum_i ||Z_i - p_i||^2
until relative loss change < 10^-4 or 20 iterations
```

The current $D=30$ fit converges in seven iterations. The helper
`project_to_polyline()` implements the segment-by-segment nearest-point
calculation. `orient_principal_curve()` optionally reverses the final direction,
and `interpolate_curve_scores()` evaluates the fitted score curve at requested
landmark coordinates. None of these latter two operations refits the curve.

### 5.5 Why the reported value is 48.0%

Define the centered total sum of squares of the input score cloud as

\[
\operatorname{TSS}
=\sum_{i=1}^n\lVert z_i-\bar z\rVert^2,
\]

and the curve residual sum of squares as

\[
\operatorname{RSS}_{\text{curve}}
=\sum_{i=1}^n\lVert z_i-p_i\rVert^2.
\]

The reported quantity is

\[
R^2_{\text{curve}}
=1-\frac{\operatorname{RSS}_{\text{curve}}}{\operatorname{TSS}}
=0.4804.
\]

Thus the fitted one-dimensional curve accounts for 48.04% of the empirical
squared Euclidean variation in the selected 30-dimensional raw score cloud.
The first linear PCA axis accounts for 46.08%, so allowing curvature adds about
1.96 percentage points. The fitted curve coordinate and the linear first-PC
coordinate have Spearman correlation 0.9967; the morphology continuum is
therefore almost linear in this dataset.

This 48.04% must not be described as any of the following:

- 48.04% of all image intensity variation;
- 48.04% of temporal variation;
- 48.04% of the complete retained MatFPCA variance; or
- a predictive clinical $R^2$.

It is a geometric goodness-of-fit measure for a one-dimensional curve in the
chosen raw score space.

### 5.6 Direction and age

Arc length has an arbitrary direction. After the curve is fitted without age,
the coordinate is replaced by $1-u_i$ if necessary so its Spearman
correlation with baseline age is nonnegative. This changes only which endpoint
is called 0 or 1; it does not change the curve, projections, subject ordering,
or 48.04% value. The resulting coordinate-age Spearman correlation is 0.4323.
Age therefore orients the coordinate but does not define it.

## 6. MRI displays along the continuum

Seven target positions are empirical percentiles from 5% through 95%. Three
displays are constructed:

1. **Curve-decoded MRI.** The score curve $g(u)$ is interpolated at the target
   position, omitted scores are set to zero, and the MatFPCA expansion is
   evaluated at baseline:

   \[
   \widehat Y_u(r,c,0)=\widehat\mu(r,c,0)
   +\sum_{\ell=1}^{D}g_\ell(u)e_\ell(r,c,0).
   \]

2. **Locally averaged observed MRI.** Each subject contributes their observed
   image closest to baseline. Gaussian kernel weights in morphology-coordinate
   distance form a local mean. The bandwidth reaches approximately 15% of the
   cohort, subject to a minimum of 0.03 coordinate units.
3. **Closest actual subject.** An unused subject closest to each target
   coordinate is selected, with curve residual and subject ID as tie-breakers.

These panels distinguish what the fitted model decodes, what nearby observed
images average to, and what a real individual scan looks like.

## 7. Reconstruction evaluation

For each retained training or held-out visit, subject scores are estimated from
the visits available to the fit, the centered image is reconstructed, and the
estimated temporal mean is restored. RMSE, relative RMSE, correlation, and PSNR
are computed only inside the training-derived brain mask. The current mean
held-out RMSE is 0.0893 and mean held-out correlation is 0.9506.

These metrics assess visit reconstruction and interpolation under the fitted
model. They do not validate the morphology coordinate against a clinical
endpoint.

## 8. Reproducible numerical outputs

The quantities discussed above are saved in:

- [`results/latent_coordinate_diagnostics_pcs30.csv`](results/latent_coordinate_diagnostics_pcs30.csv),
  including the curve and linear explained-variation values;
- [`results/latent_principal_curve_scores_pcs30.csv`](results/latent_principal_curve_scores_pcs30.csv),
  containing the fitted score-space curve;
- [`results/fpc_scores_embeddings_pcs30.csv`](results/fpc_scores_embeddings_pcs30.csv),
  containing subject coordinates, residual distances, ordered scores, and
  display embeddings;
- [`results/latent_coordinate_landmarks_pcs30.csv`](results/latent_coordinate_landmarks_pcs30.csv),
  documenting the selected continuum landmarks; and
- [`results/adni_matfpca_model_pcs30.rds`](results/adni_matfpca_model_pcs30.rds), containing
  the compact covariance fit, eigensystem, scores-derived continuum, and
  preprocessing configuration.

The concise run-level values are also recorded in
[`results/run_summary_pcs30.txt`](results/run_summary_pcs30.txt).

## 9. Sensitivity to the score dimension

The dimension is an explicit runtime argument. For example, the following fits
five choices while reusing the same image and covariance caches:

```bash
for d in 10 20 30 50 100; do
  Rscript application/adni/adni_matfpca.R --score-pcs=$d
done
```

Every continuum-derived file receives `_pcsD`, so none of these fits overwrites
another. Then run:

```bash
Rscript application/adni/summarize_continuum_sensitivity.R \
  --reference-pcs=30
```

The resulting
[`results/continuum_dimension_sensitivity.csv`](results/continuum_dimension_sensitivity.csv)
reports, for each $D$, curve and linear explained score variation, association
with age, rank and linear correlation with the $D=30$ coordinate, and absolute
coordinate differences. Coordinate correlation is the main stability measure;
the explained-variation percentages are conditional on different-dimensional
denominators and should not be compared as though their target were identical.

The current sensitivity run gives:

| $D$ | Retained product variance | Curve variation explained | Spearman correlation with $D=30$ | Mean absolute coordinate difference |
|---:|---:|---:|---:|---:|
| 10  | 8.37%  | 64.03% | 0.9893 | 0.0347 |
| 20  | 14.15% | 54.04% | 0.9983 | 0.0138 |
| 30  | 18.86% | 48.04% | 1.0000 | 0.0000 |
| 50  | 26.65% | 42.28% | 0.9981 | 0.0312 |
| 100 | 41.12% | 34.96% | 0.9960 | 0.0548 |

The ranking is exceptionally stable from 20 through 100 dimensions and remains
very similar even at 10. The declining explained percentage is expected: each
larger $D$ adds secondary directions to the total score variation in the
denominator, while the fitted object remains one-dimensional. These results
support 30 as a compact working choice, although they do not identify it as a
uniquely optimal cutoff.
