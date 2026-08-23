# ADNI BrainIAC + MatrixFPCA

This application applies MatrixFPCA to frozen BrainIAC representations from
295 ADNI subjects and 1,308 longitudinal MRI visits through month 48. Only the
full-cohort analysis is retained.

## Encoder tensors

BrainIAC uses a 12-block MONAI ViT-B/16³ encoder. Its hidden width of 768 is an
architectural checkpoint choice: the learned patch-projection layer maps each
flattened `16 × 16 × 16 = 4,096` voxel patch to 768 coordinates. It is not
estimated by MatrixFPCA.

The tensor flow for one MRI is:

```text
registered MRI: 193 × 229 × 193
  -> fixed common foreground crop: 161 × 201 × 177
  -> trilinear resize and nonzero-voxel normalization: 96 × 96 × 96
  -> 6 × 6 × 6 grid of 16³ patches
  -> learned patch projection: 216 × 768 initial token matrix
  -> Transformer Blocks 1,...,12: 216 × 768 after every block
  -> final LayerNorm: 216 × 768 final-normalized matrix
```

Thus an initial `216 × 768` patch embedding does exist before Block 1. This
analysis does not treat it as a fourth representation. It separately analyzes
the output after Block 3, after Block 6, and after Block 12 plus final LayerNorm.
Every selected layer has the same shape but different values and semantics.
There is no CLS token because the encoder uses `classification=False`; all 216
rows are spatial patch tokens.

## Common foreground crop

All source MRIs are affine-registered to the same `193 × 229 × 193` grid. The
extractor applies the same zero-based crop `[16:177, 16:217, 0:177]` to every
scan before resizing to `96³`. The resulting crop is `161 × 201 × 177` and
raises the representative nonzero fraction of the model input from about 25%
to 38%.

A subject-specific tight crop is deliberately not used: it would move the same
anatomy between token indices across subjects and visits, invalidating the
spatial MatrixFPCA mode. Every input is checked against the registered source
shape before cropping.

## Dimension retention

Feature PCA is fitted only on the 236 population-training subjects and now
selects the smallest rank retaining at least 98% variance. There is no fixed
`768 -> 64` reduction.

MatrixFPCA is not intrinsically a dimension-reduction method. In this analysis,
PVE-based truncation of the token and feature modes is disabled. Every
numerically positive estimated token/feature covariance mode is retained; only
the temporal mode uses a 98% PVE target with a cap of four. Zero or negative
estimated covariance eigenvalues cannot define valid FPCs and are excluded.

The revised full-cohort fits are:

| Layer | Encoder matrix | Feature-PCA matrix | Feature-PCA PVE | Positive token × feature × temporal ranks | Product modes |
|---|---:|---:|---:|---:|---:|
| Block 3 | `216 × 768` | `216 × 83` | 98.01% | `144 × 73 × 2` | 21,024 |
| Block 6 | `216 × 768` | `216 × 96` | 98.01% | `216 × 96 × 3` | 62,208 |
| Final normalized | `216 × 768` | `216 × 42` | 98.01% | `7 × 9 × 2` | 126 |

The final layer is still much lower-rank than the earlier layers, but it is no
longer artificially summarized by only three scores. Its 126 modes reflect the
positive rank of the estimated covariance under the revised pipeline.

For manageable tabular output, `subject_fpc_scores.csv` exports the leading 100
ordered product scores. Each FPC is one scalar per subject, so each exported
score matrix is `295 × 100`. The fitted model retains the complete positive-mode
score tensor.

## Current scientific summary

The leading morphology-associated coordinates are:

| Layer | Selected product FPC | Ventricular Spearman rho | Adjusted partial R² |
|---|---:|---:|---:|
| Block 3 | FPC7 = R2 × C2 × T1 | 0.854 after orientation | 0.673 |
| Block 6 | FPC1 = R1 × C1 × T1 | 0.700 after orientation | 0.417 |
| Final normalized | FPC1 = R1 × C1 × T1 | 0.623 after orientation | 0.313 |

Block 3 FPC7 is strongly spatially localized to ventricular tokens
(ROI-to-remote absolute-loading ratio 44.3). Label-free principal-curve
continua from the leading 20 FPCs also track ventricular morphology in all
three layers (rho 0.763, 0.612, and 0.454 for Block 3, Block 6, and final).

See [`BRAINIAC_MATFPCA_INTERPRETATION.md`](BRAINIAC_MATFPCA_INTERPRETATION.md)
for interpretation details and limitations.

## Main files

- `extract_brainiac_latents.py`: deterministic common-crop preprocessing and
  frozen Block 3, Block 6, or final-normalized token extraction.
- `adni_latent_matfpca.R`: training/validation split, 98%-PVE streaming feature
  PCA, temporal mean, untruncated positive spatial modes, scores, and diagnostics.
- `interpret_brainiac_matfpca.R`: phenotype association, crop-aware atlas/token
  localization, score stability, unsupervised continua, and observed-MRI panels.
- `compare_brainiac_layers.R`: full-cohort layer comparison.
- `adni_latent_utils.R` and `latent_io.R`: disk-backed utilities and NumPy input.
- `test_extract_brainiac_latents.py`, `test-interpret-brainiac-matfpca.R`, and
  `../../test/test-adni-latent-utils.R`: regression and source-contract tests.

## Reproduction

Create an isolated environment and install the pinned dependencies:

```bash
python -m venv application/adni-latent/.venv-brainiac
application/adni-latent/.venv-brainiac/bin/pip install \
  -r application/adni-latent/requirements-brainiac.txt
```

The local BrainIAC checkpoint is `data/models/BrainIAC.ckpt`, with SHA-256
`f22bdbcae26823a9d9e8aee883c6f24386ba4617339c12269848b6666cc62693`.
Its use must comply with the model's non-commercial academic-research terms.

Run the extractor three times, changing `LAYER` and `OUTPUT` as follows:

```text
3      application/adni-latent/latents/brainiac_patch_full_layer03
6      application/adni-latent/latents/brainiac_patch_full_layer06
final  application/adni-latent/latents/brainiac_patch_full
```

```bash
application/adni-latent/.venv-brainiac/bin/python \
  application/adni-latent/extract_brainiac_latents.py \
  --manifest data/ADNI/filtered_manifest_20260819.csv \
  --nifti-dir /home/muye/Projects/SurvADNI/ADNI_raw_data/downloaded_images/preprocessed/affine \
  --checkpoint data/models/BrainIAC.ckpt \
  --output-dir OUTPUT --max-months 48 --min-visits 3 \
  --max-subjects all --seed 20260819 --layer LAYER --device cpu --batch-size 2
```

Then run:

```bash
Rscript application/adni-latent/adni_latent_matfpca.R \
  --metadata=application/adni-latent/latents/brainiac_patch_full_layer03/metadata.csv \
  --output-dir=application/adni-latent/results/brainiac_patch_full_layer03
Rscript application/adni-latent/adni_latent_matfpca.R \
  --metadata=application/adni-latent/latents/brainiac_patch_full_layer06/metadata.csv \
  --output-dir=application/adni-latent/results/brainiac_patch_full_layer06
Rscript application/adni-latent/adni_latent_matfpca.R \
  --metadata=application/adni-latent/latents/brainiac_patch_full/metadata.csv \
  --output-dir=application/adni-latent/results/brainiac_patch_full
Rscript application/adni-latent/interpret_brainiac_matfpca.R
Rscript application/adni-latent/compare_brainiac_layers.R
```

Generated full-cohort latents and results are ignored by Git and remain under
`application/adni-latent/latents` and `application/adni-latent/results`.
