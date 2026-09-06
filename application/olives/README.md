# OLIVES retinal OCT feature experiment

This directory documents the reusable OLIVES data layout and the feature-based
longitudinal analysis. The raw images and labels are kept under `data/OLIVES`
and `data/OLIVES_Dataset_Labels`; generated features and fitted models belong in
the ignored `application/olives/output/` directory.

The local release contains 87 patients and 96 study eyes. The TREX-DME
collection has 47 patients represented by GILA, TREX, and Monthly arms. Visits
are recorded as `V1`, `V2`, … and are not interchangeable with elapsed weeks.
The image-level clinical workbook maps full image paths to `Patient_ID`,
`Eye_ID`, BCVA, and CST. Join on the full relative path, never on a basename.
Within a complete acquisition there are 49 ordered OCT B-scans; a visit is a
longitudinal observation and the B-scans are spatial positions within that
observation.

`DATA_STRUCTURE.md` contains the detailed inventory, path conventions, label
tables, and image examples. `subject_inventory.csv` and
`patient_inventory.csv` are filename inventories and are not independent
observations.

## Frozen RETFound representation

RETFound is a retinal foundation model trained with masked autoencoding: a
vision transformer learns to reconstruct masked retinal-image patches from the
visible context. We use the released OCT-specific ViT-Large checkpoint, frozen
for inference. Its encoder has 24 transformer blocks, 16 × 16 image patches,
and a 1,024-dimensional hidden state. The model is a representation learner;
we do not use its decoder, a disease-specific classifier, or any fine-tuning on
the OLIVES clinical labels.

`prepare_retfound.py` builds an outcome-blind manifest using one nominal central
B-scan (index 24) per eligible eye-visit, one study eye per patient, and the
documented visit index. It does not align, crop, downsample, or otherwise
modify the image.

`extract_retfound.py` runs the official OCT RETFound MAE ViT-L/16 encoder with
deterministic evaluation preprocessing: RGB conversion, bicubic resize of the
short side to 256, center crop to 224, and ImageNet normalization. It saves the
final encoder patch tokens, excluding CLS, as a `196 × 1024` matrix for every
eye-visit. Rows are the 14 × 14 raster patch positions and columns are learned
feature coordinates. The checkpoint is loaded frozen; the MAE decoder and any
task head are discarded. The script records checkpoint and manifest hashes and
checks repeatability against the official encoder sequence.

The resulting entries are latent activation coordinates, not segmented retinal
layers or hand-designed biomarkers. A row corresponds to one spatial 16 × 16
patch in the standardized 224 × 224 input; a column is one learned channel in
the final transformer representation. Neighboring rows preserve image-patch
order, while neighboring columns do not have a direct clinical interpretation.
The separate `pooled.npy` file averages the 196 patch vectors and is used only
for a low-dimensional exploratory screen. MatFPCA uses the full patch-by-
feature matrix, so it can model covariance across patch locations and learned
feature coordinates jointly.

The model checkpoint is the OCT weight linked by the
[RETFound repository](https://github.com/openmedlab/RETFound_MAE). The local
Compute Canada environment is used through
`~/projects/def-jiguocao/muye/env/neuro/bin/python`; the supplied Slurm script
requests one H100 GPU. The saved matrices are the only inputs to the subsequent
MatFPCA fit; raw pixel matrices are not used.

The encoder dependencies are listed in `requirements-retfound.txt`. The
checkpoint itself is intentionally kept outside Git.

```bash
module load scipy-stack/2025a
python application/olives/prepare_retfound.py
sbatch application/olives/run_retfound.sbatch
```

After extraction, `matfpca_retfound.R` fits the existing matrix-valued FPCA
implementation to the patch-by-feature trajectories and writes the retained
scores and diagnostics. `associate_retfound_scores.R` tests the leading
matrix-FPCA product scores against the clinical endpoints. `continuum_retfound.R`
forms the joint score geometry and fits an ADNI-style principal curve, including
age, BCVA, CST, and available baseline biomarker-coordinate associations.
Its primary joint test uses all nonzero score-geometry coordinates through a
distance/kernel permutation statistic; the first-ten-coordinate regression is
reported only as a sensitivity analysis because downstream truncation can lose
distributed signal.
`pca_retfound_scores.R` applies ordinary PCA to that complete score geometry
and tests every nonzero PCA mode against the same broad clinical panel. It
reports both a global correction over all mode/outcome pairs and a correction
within each outcome, so low-variance modes are not silently removed before
clinical screening. Its output also flags constant or extremely sparse binary
labels; those labels remain in the audit table but are excluded from the
screenable sensitivity family because a one-positive association can be driven
by a single patient.
`prepare_retfound_biomarkers.py` joins the sparse image-level biomarker export
at the visit level. `check_retfound_clinical.py` performs an exploratory
patient-level association screen using frozen pooled RETFound features and
clinical baseline/change variables; this is a representation check, not a
confirmatory clinical analysis.
