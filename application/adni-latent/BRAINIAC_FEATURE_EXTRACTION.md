# BrainIAC feature extraction: what it preserves and what it cannot reconstruct

## Short answer

BrainIAC is a **3D MRI encoder**, not an autoencoder. It was pretrained with a
contrastive objective to make useful MRI representations, not to reproduce its
input. The released checkpoint has no image decoder, so BrainIAC features alone
cannot be converted back to the original MRI.

There are actually two irreversible targets:

1. The native NIfTI cannot be recovered after resampling to `96 x 96 x 96` and
   per-scan intensity normalization.
2. Even the normalized `96 x 96 x 96` model input cannot be directly recovered,
   because the encoder is many-to-one and no inverse was trained.

An auxiliary decoder can be trained as a **probe** of retained information. Our
subject-held-out probe recovers only a blurry population-like brain and improves
relative RMSE by just 1--2.5% over predicting the training-mean MRI. This is not
enough evidence to treat BrainIAC latents as a lossless image representation.
It does not imply that all local anatomy is absent: the independent
[`BRAINIAC_CONCEPT_AUDIT.md`](BRAINIAC_CONCEPT_AUDIT.md) finds strongly
decodable ventricular deformation and weaker, layer-dependent hippocampal and
amygdala information.

## What BrainIAC is

[BrainIAC](https://github.com/AIM-KannLab/BrainIAC) (Brain Imaging Adaptive
Core) is a foundation encoder for structural brain MRI. Its selected backbone
was pretrained with SimCLR: augmented views of the same anatomy are pulled
together in representation space and different samples are pushed apart. This
encourages features useful for classification, regression, and segmentation;
it does **not** reward voxel reconstruction. See the
[BrainIAC paper](https://doi.org/10.1038/s41593-026-02202-6).

The model used here is:

| Part | This application |
|---|---|
| Input | one affine-registered 3D T1 MRI |
| Model input | nonzero-voxel normalized `96 x 96 x 96` volume |
| Patch embedding | `16 x 16 x 16` voxels to 768 numbers |
| Spatial grid | `6 x 6 x 6 = 216` patch tokens |
| Backbone | 12 transformer blocks, 12 attention heads, 3,072-wide MLP |
| Saved feature | all `216 x 768` tokens from block 3, block 6, or final layer |

The feature path is:

```text
native NIfTI
    -> resize to 96^3 and normalize intensities
    -> split into 216 patches
    -> project each 4096-voxel patch to 768 values
    -> 12 transformer blocks mix information across patches
    -> 216 x 768 feature matrix
```

The official wrapper returns token zero and calls it a “CLS” token. The MONAI
ViT is constructed with `classification=False`, however, so it has no CLS token;
token zero is the first spatial patch. This project keeps all 216 tokens rather
than discarding 215 of them.

## How feature extraction is used here

[`extract_brainiac_latents.py`](extract_brainiac_latents.py) freezes the
checkpoint and applies the official validation preprocessing to every visit.
It saves one finite `float32` matrix of shape `216 x 768` per MRI. Nothing is
trained and MatrixFPCA is not involved.

For downstream work, treat a row as a spatial patch representation and a column
as a learned feature. The rows are spatially ordered, but a late transformer row
is not simply a compressed copy of its corresponding image patch: self-attention
has mixed information from the whole brain into every token.

Use these features for tasks that can be validated in their own target space,
such as diagnosis, age, or outcome prediction. Do not assume that good latent
distance, latent PVE, or latent reconstruction implies good anatomy recovery.

## Why there is no direct MRI reconstruction

The local checkpoint provides 137 backbone tensors plus 12 SimCLR projection-
head tensors. It contains zero decoder, reconstruction, or UNETR tensors. The
projection head is used by the contrastive loss; it maps features to another
feature space and is not an image decoder.

Information is also discarded before and inside the encoder:

- Resizing loses native sampling information.
- Nonzero-voxel z-normalization loses native intensity location and scale unless
  those values are stored separately; the current extraction does not store them.
- Each raw patch has 4,096 voxels but its first embedding has only 768 values, so
  that linear step alone has rank at most 768.
- Later attention, nonlinearities, and normalization were never constrained to
  remain invertible.

The segmentation code in the official repository adds and trains a UNETR task
head. That predicts a segmentation mask; it is not a pretrained MRI decoder.

## Legacy information-retention audit (no MatrixFPCA)

This 30-subject probe was an early test of invertibility, not the sample used
for the current MatrixFPCA results. The full longitudinal analysis uses all 295
eligible subjects.

[`evaluate_brainiac_reconstruction.py`](evaluate_brainiac_reconstruction.py)
fits a separate linear kernel-ridge probe from the complete token matrix to the
normalized `96^3` model input. The split is by subject: 19 subjects train the
probe, 5 select its ridge penalty, and 6 subjects (25 visits) remain untouched
until final evaluation. A training-mean MRI is the required baseline. The same
split and target images are used for all layers.

| Feature layer | Mean RMSE | Mean relative RMSE | Correlation | Mean relative-RMSE reduction vs mean MRI |
|---|---:|---:|---:|---:|
| Block 3 | 0.346 | 0.688 | 0.727 | 2.48% |
| Block 6 | 0.351 | 0.698 | 0.716 | 1.15% |
| Final | 0.350 | 0.696 | 0.717 | 1.35% |
| Training-mean MRI | 0.355 | 0.706 | 0.707 | 0% |

Block 3 retains the most linearly decodable detail, but the gain over the mean
MRI is small. The examples show smooth, population-like reconstructions; fine
cortical boundaries and subject-specific ventricular detail remain inaccurate.

- [Block 3 examples](reconstruction_probe/reconstruction_examples_block_3.png)
- [Block 6 examples](reconstruction_probe/reconstruction_examples_block_6.png)
- [Final-layer examples](reconstruction_probe/reconstruction_examples_final.png)
- [Aggregate metrics](reconstruction_probe/reconstruction_summary.csv)

This probe is a lower-bound test of linear decodability, not an official or
generative BrainIAC decoder. BrainIAC pretraining also included ADNI, so the
subject split isolates probe fitting but cannot undo possible encoder pretraining
overlap. Most importantly, a probe trained on observed BrainIAC features is not
automatically valid for an arbitrary MatrixFPCA reconstruction that may lie off
the observed feature manifold.

## Consequence for MatrixFPCA

MatrixFPCA models a trajectory in **BrainIAC feature space**. Its latent errors
would answer whether that feature trajectory is reproduced, not whether the
original anatomy is reproduced. The current analysis therefore disables latent
reconstruction. Because every latent is linked to its source scan, the
full-cohort continuum panels instead show each selected subject's actual
observed baseline MRI. These panels validate anatomical ordering but are not
decoded reconstructions.

The full-cohort interpretation provides the validation bridge: rather than
decoding every voxel, it tests learned coordinates against independently
computed regional morphology levels and slopes and displays actual observed
MRIs along those coordinates.

Therefore:

- Keep BrainIAC + MatrixFPCA only if the scientific estimand is explicitly a
  trajectory of predictive BrainIAC representations and validate it with a
  downstream clinical task.
- If image-space reconstruction is essential, switch the primary representation
  to a model with a compatible pretrained image decoder (for example, an MRI
  autoencoder or masked autoencoder), and score held-out reconstructions in the
  normalized image space before applying MatrixFPCA.
- If a new decoder is trained for BrainIAC, validate it on held-out subjects and
  then test decoded MatrixFPCA latents for off-manifold failures. Do not describe
  the result as native-space recovery unless resampling and intensity transforms
  are also invertible and their per-scan parameters are retained.

## Reproduce the audit

With the three existing extraction caches:

```bash
application/adni-latent/.venv/bin/python \
  application/adni-latent/evaluate_brainiac_reconstruction.py
```

The script writes anonymous visit-level metrics, tuning results, a provenance
manifest, and three example contact sheets to `reconstruction_probe/`.

For the complementary global/local anatomical audit, run:

```bash
application/adni-latent/.venv/bin/python \
  application/adni-latent/evaluate_brainiac_concepts.py
```
