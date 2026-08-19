# ADNI MatFPCA application

`adni_matfpca.R` runs the complete longitudinal analysis on the axial ADNI
slices and the newest `data/ADNI/filtered_manifest_*.csv` file.

The pipeline applies these design rules before fitting:

- it keeps observations at or before 48 months and retains screening visits;
- it measures time in years, snaps dates at one-day resolution, assigns every
  visit coded `bl` to exactly zero, and starts the grid below zero so screening
  images are represented. Only observed grid locations plus 0 and 4 years are
  retained, avoiding a wasteful dense daily covariance matrix;
- it removes subjects with fewer than two available images;
- it chooses a reproducible random 20% of subjects and holds out one randomly
  selected image from each, without using a clinical outcome label;
- it estimates a subject-weighted spline mean image using training images only,
  subtracts the mean at each visit before MatFPCA, and adds it back to every
  reconstruction;
- it infers a consensus brain mask from training-image foreground boundaries,
  requiring support in at least half of training images, keeps its largest
  connected component, fills internal holes, and applies the mask before fitting
  and after reconstruction;
- it block-averages the images to 114 x 96 by default. This retains substantially
  more anatomical detail while keeping three-fold covariance cross-fitting
  practical;
- it evaluates the full temporal eigensystem even when marginal-PVE truncation
  retains only one time mode. A trapezoidal quadrature correction accounts for
  the irregular compact time grid. For each numerically supported mode, it
  reports marginal PVE and the fraction of time-domain energy left after
  subtracting the best constant approximation. A diagnostic screen calls a mode
  materially dynamic when these quantities are at least 0.1% and 50%,
  respectively;
- it fits a fixed-complexity one-dimensional principal curve through 30 leading
  raw product-FPC scores. The raw scale makes dominant, high-variance modes count
  more than weak modes. Arc length along that curve gives every subject a continuous
  morphology coordinate in [0,1]. Age is not used to learn the coordinate; it
  only fixes the otherwise arbitrary direction so the age correlation is
  nonnegative;
- it displays six empirical-percentile landmarks in three complementary ways:
  smooth curve-decoded baseline MRIs, kernel-smoothed means of observed
  near-baseline MRIs, and the closest actual subject MRI at each landmark;
- it compares PCA, UMAP, and t-SNE in two dimensions. UMAP is the primary
  local-structure display, PCA checks global geometry, and t-SNE is a local
  neighborhood sensitivity analysis.

Within this script, clinical outcome groups are not used for subject sampling,
holdout selection, continuum fitting, embedding, plotting, or diagnostics. The
input cohort still comes from the supplied filtered manifest.

See [`TECHNICAL_DETAILS.md`](TECHNICAL_DETAILS.md) for the mathematical model,
estimation steps, score construction, and a derivation of the continuous
morphology coordinate and its 48.0% explained score variation.

Run the full analysis from the repository root:

```bash
Rscript application/adni/adni_matfpca.R
```

The default estimator uses three-fold cross-fitted Gaussian-refined weights.
For a quick smoke test:

```bash
Rscript application/adni/adni_matfpca.R \
  --max-subjects=30 \
  --pool=4 \
  --weight-scheme=identity \
  --output-dir=/tmp/adni_matfpca_smoke
```

Use `--help` for the main runtime controls. Image and model caches are stored
below the selected result directory and are reused only when the analysis
signature matches. The default result directory is `application/adni/results`.

## Outputs

The result directory contains:

- `analysis_manifest.csv`: every retained image, snapped year, and train/test
  assignment;
- `reconstruction_metrics.csv` and `reconstruction_summary.csv`: training and
  held-out RMSE, relative RMSE, correlation, and PSNR, evaluated inside the
  inferred brain mask;
- `inferred_brain_mask.png`: the training-only consensus mask used by the fit
  and reconstructions;
- `top_fpc_spatiotemporal_panels.png`: a publication-oriented 2 x 5 display.
  Each row begins with one of the two variance-bearing temporal eigenfunctions
  and then shows the four leading spatial product eigenfunctions. Captions give
  separate marginal FVE, and all panels of a given type share a common scale;
- `temporal_eigencomponent_diagnostics.csv` and
  `temporal_eigencomponent_diagnostics.png`: full-eigensystem evidence for
  distinguishing variance-bearing temporal dynamics from non-constant but
  numerically negligible modes;
- `fpc_scores_embeddings_pcs30.csv`: subject age, continuous morphology coordinate,
  curve residual, all ordered FPC scores, and PCA/UMAP/t-SNE coordinates;
- `latent_coordinate_diagnostics_pcs30.csv`,
  `latent_principal_curve_scores_pcs30.csv`, and
  `latent_coordinate_landmarks_pcs30.csv`: the fitted one-dimensional curve, its
  linear-axis comparison, age association, and selected MRI landmarks;
- `fpc_score_umap_by_continuum_pcs30.png` and
  `fpc_score_embedding_comparison_by_continuum_pcs30.png`: continuous-color
  views of the subjects;
- `smoothed_observed_mri_continuum_pcs30.png`,
  `observed_mri_continuum_pcs30.png`, and
  `reconstructed_mri_continuum_pcs30.png`: population-smoothed, actual-subject, and
  MatFPCA-decoded MRI changes from the low to high end of the coordinate;
- `latent_coordinate_vs_age_pcs30.png`: the continuous, imperfect relationship with
  baseline age;
- `embedding_quality_metrics_pcs30.csv`: PCA, UMAP, and t-SNE trustworthiness,
  neighbor overlap, and global distance correlation;
- `continuum_dimension_sensitivity.csv`: stability of coordinates and explained
  score variation across the fitted score dimensions;
- `adni_matfpca_model_pcs30.rds`: the compact fitted covariance, retained
  eigensystem, mean function, grid, morphology continuum, embeddings, and
  preprocessing configuration; and
- `run_summary_pcs30.txt` and `session_info_pcs30.txt`: the principal numerical results and
  software versions.

The `_pcs30` suffix records the actual number of score dimensions used. Change
`--score-pcs` to fit another dimension without overwriting the default
continuum. After fitting multiple dimensions, run
`summarize_continuum_sensitivity.R` to create
`continuum_dimension_sensitivity.csv`.

The nonlinear two-dimensional plots are visualization aids. The morphology
coordinate itself is fitted in the original selected-dimensional raw score
space; 30 dimensions is the default.
