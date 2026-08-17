# ADNI MatFPCA application

`adni_matfpca.R` runs the complete longitudinal analysis on the axial ADNI
slices and the newest `data/ADNI/filtered_manifest_*.csv` file.

The pipeline applies the requested design rules before fitting:

- it keeps observations at or before 48 months and retains screening visits;
- it measures time in years, snaps dates at one-day resolution, assigns every
  visit coded `bl` to exactly zero, and starts the grid below zero so screening images are
  represented. Only observed grid locations plus 0 and 4 years are retained,
  avoiding a wasteful dense daily covariance matrix;
- it removes subjects with fewer than two available images;
- it chooses a reproducible, outcome-stratified 20% of subjects and holds out
  one randomly selected image from each;
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
- it embeds 30 standardized product-FPC scores into two dimensions with UMAP.
  Patient labels are used only afterward for color and representation checks.

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
- `reconstruction_examples_by_group.png`: a representative baseline-AD,
  converting-MCI, and non-converting-MCI subject, with observed/reconstructed
  training images and observed/predicted held-out images;
- `top_fpc_spatial_panels.png` and `top_fpc_temporal_panels.png`: the leading
  product FPC effects in one-row panels;
- `fpc_scores_and_umap.csv`: ordered FPC scores and the two UMAP coordinates;
- `fpc_score_umap_by_patient_group.png`: the label-free UMAP embedding colored
  afterward by clinical outcome group;
- `umap_representation_diagnostics.csv`, `umap_neighbor_diagnostics.csv`, and
  `knn_probe_*.csv`: known-group silhouettes, neighborhood enrichment, and a
  stratified cross-validated k-nearest-neighbor probe in the original FPC space.
  The probe is reported both for all three groups and for converter versus
  non-converter MCI subjects, since the baseline-AD group is much smaller;
- `adni_matfpca_model.rds`: the compact fitted covariance, retained
  eigensystem, mean function, grid, and preprocessing configuration;
- `run_summary.txt` and `session_info.txt`: the principal numerical results and
  software versions.

The clinical comparison uses three manifest-defined outcome groups: prevalent
AD at baseline, incident MCI-to-AD conversion, and right-censored non-converting
MCI. These labels are not supplied to MatFPCA or UMAP. UMAP preserves local
neighborhoods rather than global distances and can create visually compelling
islands even without clinically meaningful separation. The FPC-space silhouette,
neighbor enrichment, and cross-validated probe should therefore be considered
alongside the plot; none turns this exploratory analysis into a diagnostic model.
