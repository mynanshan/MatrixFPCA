# ADNI BrainIAC + MatrixFPCA

This directory studies whether frozen BrainIAC features preserve useful global
and local brain morphology, and whether MatrixFPCA can aggregate those
longitudinal matrix features into interpretable, fixed-shape subject scores.

Start with:

- [`BRAINIAC_FEATURE_EXTRACTION.md`](BRAINIAC_FEATURE_EXTRACTION.md) for an
  intuitive description of BrainIAC and its `216 x 768` feature matrix;
- [`BRAINIAC_CONCEPT_AUDIT.md`](BRAINIAC_CONCEPT_AUDIT.md) for the independent
  no-MatrixFPCA global/local information probes;
- [`BRAINIAC_MATFPCA_INTERPRETATION.md`](BRAINIAC_MATFPCA_INTERPRETATION.md) for
  the current full-cohort scientific results.

## Current analysis: all eligible subjects

The reported analysis uses all 295 eligible ADNI subjects and 1,308 MRI visits
through month 48. The earlier 30-subject run was only a computational pilot;
its correlations and nonlinear plots are not used in the current conclusions.

- 236 population-training subjects and 1,033 visits fit feature PCA, the mean,
  and MatrixFPCA covariance.
- 59 subjects are held out from population fitting; 58 contribute an interior
  held-out visit for model diagnostics.
- Block 3, Block 6, and final features use exactly the same subject-visit keys.
- Every visit is represented by all 216 spatial tokens and 768 BrainIAC
  features before training-only feature PCA.
- Image reconstruction is disabled. Anatomical panels show original observed
  MRIs for subjects selected along learned coordinates.

The full fitted models are:

| Layer | Token x feature x time ranks | Product modes | Feature-PCA PVE | Separability correlation |
|---|---:|---:|---:|---:|
| Block 3 | 24 x 24 x 1 | 576 | 96.91% | 0.250 |
| Block 6 | 3 x 7 x 2 | 42 | 96.95% | 0.950 |
| Final normalized | 1 x 3 x 1 | 3 | 99.22% | 0.995 |

The full cohort materially changes the interpretation. Block 3 contains a
strongly localized ventricular mode; Block 3 and final-layer score continua
track ventricular enlargement and atrophy in actual observed MRIs. Block 6 has
anatomy-associated individual modes, but its unsupervised principal curve is
not a clear physical continuum. Several regional longitudinal-slope
associations now survive strict layer-wide multiplicity correction.

## Files

- `extract_brainiac_latents.py`: deterministic, resumable frozen-encoder
  extraction, exact checkpoint loading, atomic NumPy writes, hashes, and
  provenance. Its default is now all eligible subjects.
- `adni_latent_matfpca.R`: training/validation split, streaming feature PCA,
  temporal mean, MatrixFPCA fitting, scores, and diagnostics.
- `interpret_brainiac_matfpca.R`: phenotype association, eigenmode
  localization, score stability, unsupervised continua, embeddings, and
  observed-MRI landmarks. It requires at least 250 subjects by default.
- `compare_brainiac_layers.R`: full-cohort layer summary and comparison plot.
- `adni_latent_utils.R` and `latent_io.R`: disk-backed utilities and a
  dependency-free reader for C-order NumPy matrices.
- `evaluate_brainiac_reconstruction.py`: legacy no-MatrixFPCA linear
  decodability audit; this is not a BrainIAC decoder.
- `evaluate_brainiac_concepts.py`: subject-held-out global/local anatomical
  probes with matched distant-token controls; it does not use MatrixFPCA.
- `test_extract_brainiac_latents.py` and
  `test-interpret-brainiac-matfpca.R`: extraction and interpretation safety
  tests.

The core disk-backed MatrixFPCA path also required a fix in
[`../../R/matfdata.R`](../../R/matfdata.R), with a regression test in
[`../../test/test-disk-backed-matfdata.R`](../../test/test-disk-backed-matfdata.R).

## Reproduce the full extraction

The ignored environment in `.venv` contains the validated local dependencies.
For a clean environment, use a separate path:

```bash
python -m venv application/adni-latent/.venv-brainiac
application/adni-latent/.venv-brainiac/bin/pip install \
  -r application/adni-latent/requirements-brainiac.txt
```

The official checkpoint is locally stored, and ignored by Git, at
`data/models/BrainIAC.ckpt`. Its SHA-256 is:

```text
f22bdbcae26823a9d9e8aee883c6f24386ba4617339c12269848b6666cc62693
```

Confirm that use and redistribution comply with BrainIAC's non-commercial
academic-research terms. Never commit the checkpoint or an access token.

Run one extraction for each layer. `--max-subjects all` is explicit below for
clarity, although all eligible subjects is now the default.

```bash
application/adni-latent/.venv/bin/python \
  application/adni-latent/extract_brainiac_latents.py \
  --manifest data/ADNI/filtered_manifest_20260819.csv \
  --nifti-dir /home/muye/Projects/SurvADNI/ADNI_raw_data/downloaded_images/preprocessed/affine \
  --checkpoint data/models/BrainIAC.ckpt \
  --output-dir application/adni-latent/latents/brainiac_patch_full_layer03 \
  --max-months 48 --min-visits 3 --max-subjects all \
  --seed 20260819 --layer 3 --device cpu --batch-size 2

application/adni-latent/.venv/bin/python \
  application/adni-latent/extract_brainiac_latents.py \
  --manifest data/ADNI/filtered_manifest_20260819.csv \
  --nifti-dir /home/muye/Projects/SurvADNI/ADNI_raw_data/downloaded_images/preprocessed/affine \
  --checkpoint data/models/BrainIAC.ckpt \
  --output-dir application/adni-latent/latents/brainiac_patch_full_layer06 \
  --max-months 48 --min-visits 3 --max-subjects all \
  --seed 20260819 --layer 6 --device cpu --batch-size 2

application/adni-latent/.venv/bin/python \
  application/adni-latent/extract_brainiac_latents.py \
  --manifest data/ADNI/filtered_manifest_20260819.csv \
  --nifti-dir /home/muye/Projects/SurvADNI/ADNI_raw_data/downloaded_images/preprocessed/affine \
  --checkpoint data/models/BrainIAC.ckpt \
  --output-dir application/adni-latent/latents/brainiac_patch_full \
  --max-months 48 --min-visits 3 --max-subjects all \
  --seed 20260819 --layer final --device cpu --batch-size 2
```

## Reproduce MatrixFPCA and interpretation

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

Generated features and results are ignored by Git. Current outputs are:

- `latents/brainiac_patch_full_layer03`,
  `latents/brainiac_patch_full_layer06`, and `latents/brainiac_patch_full`;
- `results/brainiac_patch_full_layer03`,
  `results/brainiac_patch_full_layer06`, and `results/brainiac_patch_full`;
- `results/brainiac_matfpca_interpretation_full`;
- `results/brainiac_layer_comparison_full`.

Each fitted-layer directory contains its model, manifest, scores,
eigenspectra, temporal functions, token maps, and diagnostics. The
interpretation directory contains figures and machine-readable association,
localization, continuum, embedding, and observed-MRI landmark tables.

## What the MRI panels mean

BrainIAC is an encoder without a compatible image decoder, so its feature
matrices are not inverted into generated MRIs. This does not prevent anatomical
validation: each feature vector retains its source subject and visit, allowing
the analysis to display that subject's original MRI. The continuum panels
select subjects at score-coordinate landmarks and show their actual baseline
images. They test whether a learned ordering corresponds to visible anatomy;
they are not reconstruction claims.

The earlier linear probe found that coarse population structure is decodable
but fine image detail is not reliably reconstructed. The concept audit asks the
more relevant question—whether global and local anatomical quantities remain
predictable from frozen features. Both were exploratory 30-subject,
pre-MatrixFPCA audits and are retained as provenance, not as the evidence for
the current full-cohort MatrixFPCA conclusions.

## Legacy pilot outputs

Directories ending in `_pilot`, the old `brainiac_layer_comparison` directory,
and `brainiac_matfpca_interpretation` are retained only for provenance. Default
scripts now point to the `_full` data/results, and the interpretation guard
prevents accidental reporting of the 30-subject pilot.

## Coronal orientation utility

The old PNG exporter requires only a top-to-bottom flip. It must not be rotated
180 degrees or horizontally flipped. `fix_coronal_orientation.py` writes
corrected copies without modifying source images.
