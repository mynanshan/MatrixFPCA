# Guide: Latent-Matrix Sparse Functional PCA for Longitudinal Medical Imaging

## 1. Goal

The purpose of these experiments is to test the matrix-valued sparse functional-data method on **learned longitudinal representations** rather than directly on raw pixels.

For subject \(i\) at irregular visit time \(t_{ij}\), start with a medical image or volume

\[
X_i(t_{ij}),
\]

apply a **fixed pretrained image encoder** \(E\), and retain a structured matrix representation

\[
Z_i(t_{ij}) = E\{X_i(t_{ij})\}\in\mathbb{R}^{p\times q}.
\]

Then apply the proposed sparse matrix-valued FPCA/covariance model to

\[
\{(t_{ij},Z_i(t_{ij})):j=1,\ldots,m_i\}.
\]

The key methodological division of labor is:

- the deep model performs **representation learning / nuisance suppression**;
- the statistical model performs **population-level longitudinal covariance estimation and FPCA under sparse irregular sampling**.

The encoder should initially be treated as **fixed**. Do not jointly train the encoder and FPCA model in the first analysis.

---

## 2. Shared principles for both experiments

### 2.1 Freeze the encoder

Use the pretrained encoder in deterministic inference mode.

Do not:

- fine-tune on diagnosis;
- fine-tune using time;
- jointly optimize the encoder with the FPCA objective;
- use random test-time augmentation when extracting the final representation.

This keeps the statistical problem close to the current theory: after feature extraction, \(Z_i(t)\) is treated as the observed matrix-valued functional datum.

### 2.2 Prefer a matrix with meaningful modes

Avoid taking a single pooled vector unless it is only a baseline.

Preferred forms are:

\[
Z_i(t)\in\mathbb{R}^{\text{spatial unit}\times\text{feature}},
\]

where the first mode retains anatomical/spatial organization and the second mode contains learned feature coordinates.

Two particularly useful constructions are:

1. **slice/scan-position × latent feature**, obtained by stacking one embedding from each ordered slice;
2. **patch/token × latent feature**, obtained by retaining all spatial patch tokens from one 2D image.

These have a natural Kronecker interpretation:

\[
\operatorname{Cov}\{\operatorname{vec} Z_i(s),\operatorname{vec} Z_i(t)\}
\approx
K(s,t)\,V_{\mathrm{feature}}\otimes U_{\mathrm{spatial}}.
\]

### 2.3 Do not assume the latent representation is automatically well behaved

Before the full analysis, explicitly check:

- temporal signal strength;
- within-subject versus between-subject variation;
- approximate Kronecker separability;
- eigenvalue/effective-rank behavior of both matrix modes;
- number of FPCs required to explain a reasonable fraction of latent variance.

The representation is useful only if it makes the longitudinal statistical problem cleaner than the raw-image representation.

### 2.4 Use subject-level splits whenever a prediction/validation task is performed

All visits from one subject must stay in the same split.

For repeated eyes from the same patient, either:

- retain only one eye per patient for the primary analysis, or
- account for patient-level clustering in a sensitivity analysis.

### 2.5 Store representations before fitting the statistical model

Feature extraction should be a separate preprocessing stage.

A useful data contract is:

```text
metadata.parquet
    subject_id
    visit_id
    time
    image_path
    diagnosis / group
    auxiliary clinical variables
    encoder_name
    encoder_checkpoint
    representation_type

latents/
    <subject_id>/
        <visit_id>.npy
```

Each `.npy` file should contain one `float32` matrix of identical shape `(p, q)` across all visits.

Also save an extraction manifest containing:

- encoder commit/checkpoint hash;
- preprocessing choices;
- input resolution;
- layer used;
- token type used;
- whether CLS/storage tokens were removed;
- matrix dimensions;
- feature centering/scaling choices.

---

# 3. Experiment A: ADNI longitudinal T1 MRI + frozen BrainDINO

## 3.1 Scientific question

Can a brain-MRI-specific frozen representation expose longitudinal neuroanatomical variation more efficiently than raw image pixels, so that sparse matrix-valued FPCA identifies meaningful population and subject-specific trajectories?

The main application should emphasize longitudinal structural progression rather than diagnostic classification.

Possible downstream clinical quantities include:

- diagnosis / conversion status;
- MMSE;
- ADAS-Cog;
- CDR or CDR-SB;
- hippocampal or ventricular measurements if available.

These are for interpretation/validation, not for training the encoder.

---

## 3.2 Data choice

Use longitudinal **T1-weighted structural MRI** from ADNI.

For the first pilot, do not use the whole ADNI database. Construct a manageable cohort such as:

- subjects with at least 3 usable T1 visits;
- at least approximately 2 years between first and last usable visit where possible;
- preferential enrichment for baseline MCI, AD, or clinically progressing subjects;
- optionally keep a smaller cognitively-normal/stable group as a comparison.

Define functional time as

\[
t_{ij}
=
\frac{\text{scan date}_{ij}-\text{baseline scan date}_i}
     {365.25}
\]

in years.

A good initial goal is roughly 50–150 subjects rather than maximizing sample size.

---

## 3.3 MRI preprocessing

The first implementation should be deliberately conservative.

Recommended pipeline:

1. select the desired ADNI T1 acquisition per visit;
2. apply ADNI quality-control information;
3. skull-strip;
4. orient volumes consistently;
5. crop to brain foreground;
6. resize/resample to a common volume grid;
7. normalize intensity;
8. optionally apply rigid or affine template alignment;
9. avoid deformable registration in the first version.

A reasonable BrainDINO-compatible preprocessing target is a standardized \(256^3\) volume followed by axial slice extraction.

Important principle:

> Use alignment to remove gross scanner/pose nuisance, but do not aggressively deform every brain into the same anatomy if that risks removing the biological longitudinal variation of interest.

---

## 3.4 Encoder

Primary intended encoder:

**BrainDINO**, frozen.

Architecture-level facts relevant to this project:

- slice-wise brain MRI foundation model;
- ViT-B/16 backbone;
- 16×16 patches;
- trained on axial 2D brain MRI slices;
- learns both CLS-level global representations and patch-level representations;
- suitable for frozen-feature extraction.

### Availability gate

Before implementing the full experiment, check whether the pretrained BrainDINO checkpoint has become publicly downloadable.

At the time this guide was prepared, the repository contained downstream code but stated that pretrained weights were not yet publicly released.

If the checkpoint is unavailable:

1. request the checkpoint from the authors if appropriate;
2. prototype the complete pipeline using a frozen alternative such as DINOv3 or an available brain-MRI foundation model;
3. keep the representation-extraction interface model-agnostic so BrainDINO can later be swapped in without changing the FPCA code.

Do not spend substantial effort reproducing BrainDINO pretraining.

---

## 3.5 Preferred latent representation: axial-slice × latent-feature matrix

This is the preferred ADNI representation.

For each visit, extract a fixed ordered set of axial slices from the standardized volume:

\[
x_{ij1},\ldots,x_{ijS}.
\]

For each slice, extract one frozen slice-level global representation, e.g. the final-layer CLS token:

\[
c_{ijs}\in\mathbb R^D.
\]

Stack them in anatomical order:

\[
Z_i(t_{ij})
=
\begin{bmatrix}
c_{ij1}^{\top}\\
\vdots\\
c_{ijS}^{\top}
\end{bmatrix}
\in\mathbb R^{S\times D}.
\]

Recommended first choice:

```text
row mode    = standardized axial slice position
column mode = latent feature coordinate
```

A natural value is \(S=128\), matching the slice-wise representation strategy used by BrainDINO, but the code should make this configurable.

### Why this representation is attractive

It:

- uses the entire 3D MRI rather than arbitrarily selecting one 2D slice;
- remains a matrix, so no tensor extension is needed;
- retains coarse superior–inferior anatomical organization;
- removes enormous pixel-level dimensionality;
- gives a scientifically interpretable row covariance \(U\) across axial positions;
- gives a learned feature covariance \(V\) across representation coordinates.

BrainDINO's own subject-level analyses average slice embeddings into a vector. Here, **do not average them**; preserving the slice-by-feature matrix is the central idea.

---

## 3.6 Alternative ADNI representation: patch-token × feature matrix from one slice

For a simpler 2D version, choose a fixed anatomically corresponding axial slice per visit.

Pass it through the frozen encoder and retain the patch tokens from one selected transformer block:

\[
Z_i(t_{ij})
\in
\mathbb R^{P\times D},
\]

where

```text
P = number of spatial patch tokens
D = token embedding dimension
```

Exclude:

- CLS token;
- storage/register tokens;
- any task-specific head.

For a ViT-B/16 encoder, \(P\) is determined by the input resolution and patch size. Infer the final tensor shape from the implementation rather than hard-coding it.

Interpretation:

```text
row mode    = spatial patch location
column mode = latent feature coordinate
```

Use this version if the slice-stack representation has poor cross-subject anatomical correspondence.

---

## 3.7 Which transformer layer to use?

Start with the final encoder block.

As a sensitivity analysis, compare:

- an early block;
- a middle block;
- the final block.

The desired layer should balance:

- spatial/anatomical information;
- disease sensitivity;
- longitudinal smoothness;
- avoidance of overly task-like or overly invariant representations.

Do not select the layer based solely on downstream diagnosis accuracy. Prefer the layer that gives the cleanest longitudinal structure and acceptable separability.

---

## 3.8 Latent normalization before FPCA

Do not normalize every visit independently to unit Frobenius norm by default; that could erase scientifically meaningful amplitude variation.

Recommended first pass:

1. extract raw frozen features;
2. estimate a grand mean matrix from the analysis cohort;
3. center the latent matrices;
4. optionally standardize feature columns if their scales are extremely heterogeneous;
5. keep a no-column-scaling version as the primary analysis if the encoder dimensions are already comparable.

Any centering/scaling constants used for held-out prediction must be estimated from the training subjects only.

---

## 3.9 FPCA analysis

Fit the proposed sparse matrix-valued covariance model to

\[
\{Z_i(t_{ij})\}.
\]

Estimate:

\[
U,\qquad V,\qquad K(s,t),\qquad \sigma^2
\]

or their exact counterparts in the current implementation.

For the slice-stack representation:

- \(U\): covariance across axial positions;
- \(V\): covariance across latent feature coordinates;
- \(K\): longitudinal covariance over years.

Report:

- row/feature covariance summaries;
- eigenvalue spectra;
- effective ranks;
- temporal covariance surface;
- temporal eigenfunctions;
- leading separable FPC structures;
- subject-level scores/trajectories.

---

## 3.10 ADNI diagnostics and validation

### Diagnostic A: Is there visible temporal signal in latent space?

Compute within-subject latent distances

\[
d_{ijk}
=
\frac{\|Z_i(t_{ij})-Z_i(t_{ik})\|_F}
{\sqrt{SD}}
\]

and examine them versus

\[
|t_{ij}-t_{ik}|.
\]

Compare with between-subject distances.

A useful representation should show nontrivial, reasonably smooth dependence on time lag.

### Diagnostic B: Does latent space improve over raw pixels?

On a matched subset, compare:

1. downsampled raw-image matrix FPCA;
2. pooled latent vector FPCA;
3. proposed latent-matrix FPCA.

Useful comparisons:

- effective rank;
- number of components required for 80/90/95% latent-space PVE;
- stability of estimated temporal covariance;
- association of FPC scores with clinical variables;
- held-out-visit prediction in the representation space.

Do not interpret lower-dimensional latent PVE as proof of better image reconstruction; the spaces have different meanings.

### Diagnostic C: Kronecker approximation

Compute a diagnostic measure of how well the empirical latent covariance is approximated by the fitted separable covariance.

The goal is not perfect separability. The key question is whether the approximation is plausible enough that the structure is useful rather than grossly misspecified.

### Diagnostic D: Clinical interpretability

Relate leading subject scores or reconstructed trajectories to:

- MCI/AD status;
- conversion;
- MMSE;
- ADAS-Cog;
- CDR-SB.

These analyses are primarily external validation of the learned functional variation.

---

# 4. Experiment B: OLIVES longitudinal OCT + frozen RETFound

## 4.1 Scientific question

Can retina-specific learned representations turn longitudinal OCT into a cleaner matrix-valued functional process, with FPC trajectories reflecting treatment- and disease-related structural changes?

This application is attractive because it has:

- genuinely repeated longitudinal observations;
- irregular clinical/treatment timing;
- a relatively small cohort;
- retina-specific clinical and biomarker annotations;
- a directly available retinal foundation model.

---

## 4.2 Data

Use the **OLIVES** dataset.

Primary modality:

**OCT**, not fundus photography.

Reasons:

- OCT anatomy has a more standardized spatial organization than retinal vessel topology in fundus images;
- retinal-layer/fluid morphology is directly relevant to disease and treatment response;
- RETFound provides a dedicated OCT encoder;
- the fundus images can remain available for optional qualitative comparison.

Define subject as an eye unless patient-level pairing requires a stricter independence treatment.

Define time in weeks (or days) from each eye's baseline visit.

Preserve:

- diagnosis;
- treatment/injection timing;
- available biomarker labels;
- clinical outcome variables.

Do not initially attempt causal treatment modeling. Treatment timing should be retained as metadata for interpretation.

---

## 4.3 Encoder

Use the **RETFound OCT encoder**, frozen.

Recommended checkpoint:

```text
RETFound OCT ViT-Large
```

Relevant architecture:

```text
patch size        = 16 × 16
encoder            = ViT-Large
embedding dimension = 1024
```

Use the official pretrained OCT weights.

Do not fine-tune RETFound for DR/DME labels in the primary analysis. The goal is a generic retinal representation, not a disease-supervised embedding.

---

## 4.4 Preferred OLIVES representation: central-B-scan patch × feature matrix

This should be the first OLIVES implementation because it stays closest to RETFound's documented OCT preprocessing.

For each longitudinal OCT volume:

1. identify the central/middle B-scan using the dataset's acquisition ordering;
2. resize/preprocess it exactly as expected by RETFound;
3. perform deterministic inference;
4. extract the final encoder patch tokens;
5. discard the CLS token if present;
6. stack the patch tokens as rows.

With a 224×224 input and 16×16 patches, the natural matrix is

\[
Z_i(t_{ij})
\in
\mathbb R^{196\times1024}.
\]

Interpretation:

```text
row mode    = ordered retinal image patches
column mode = learned retinal feature coordinates
```

Do not reshape the 196 patch locations into a vector-valued global embedding. The preservation of the token-by-feature matrix is the point of the analysis.

---

## 4.5 Alternative OLIVES representation: B-scan-position × feature matrix

A second, more ambitious representation uses the whole OCT volume.

For each visit, let the ordered OCT B-scans be

\[
x_{ij1},\ldots,x_{ijB}.
\]

Extract one global frozen RETFound feature per B-scan:

\[
c_{ijb}\in\mathbb R^{1024},
\]

then stack:

\[
Z_i(t_{ij})
\in
\mathbb R^{B\times1024}.
\]

If all visits contain at least 49 B-scans, a natural first target is:

```text
B = 49
```

If more scans are available, resample/select a common ordered set of 49 scan positions.

Interpretation:

```text
row mode    = OCT B-scan position across the volume
column mode = latent retinal feature
```

### Important validation before using this version

Confirm that:

- B-scan ordering is consistent across visits;
- central and peripheral scan positions correspond approximately across eyes;
- RETFound features remain sensible away from the central scan;
- the row-mode covariance is not dominated by acquisition-order artifacts.

If these checks fail, use the central-B-scan patch-token representation instead.

---

## 4.6 OCT preprocessing

Follow the official RETFound preprocessing as closely as possible.

The deterministic inference pipeline should approximately include:

1. extract selected OCT B-scan;
2. resize using the model's expected preprocessing;
3. normalize consistently with RETFound;
4. avoid random crop/flip at final feature extraction;
5. run frozen encoder;
6. save tokens before the classification head.

Do not perform complicated deformable registration unless a pilot clearly shows it is required.

The latent representation is partly motivated by avoiding an excessive raw-image registration burden.

---

## 4.7 FPCA analysis

Apply the same sparse matrix-valued functional model:

\[
\operatorname{Cov}\{\operatorname{vec}Z_i(s),
                    \operatorname{vec}Z_i(t)\}
\approx
K(s,t)\,V\otimes U.
\]

For the central-B-scan representation:

- \(U\): covariance across spatial retinal patches;
- \(V\): covariance across latent retinal features;
- \(K\): longitudinal/treatment-period covariance.

For the whole-volume representation:

- \(U\): covariance across B-scan positions;
- \(V\): covariance across latent features;
- \(K\): longitudinal covariance.

---

## 4.8 OLIVES diagnostics and validation

### Diagnostic A: longitudinal change

Plot normalized latent distance against visit lag.

Also stratify or annotate by:

- treatment intervals;
- biomarker change;
- diagnosis;
- clinical response.

### Diagnostic B: FPC-score interpretation

Test whether leading FPC scores/trajectories are associated with known retinal biomarkers or clinical measurements.

The objective is not necessarily prediction. The important question is whether the statistical modes correspond to recognizable longitudinal disease/treatment variation.

### Diagnostic C: latent matrix versus raw OCT

Compare:

1. raw central B-scan FPCA;
2. global RETFound embedding;
3. RETFound patch-token matrix FPCA.

Evaluate:

- effective ranks;
- temporal covariance stability;
- number of components required for target PVE;
- held-out-visit latent prediction;
- association with biomarkers.

### Diagnostic D: treatment-time sensitivity

Because OLIVES subjects undergo repeated treatment, inspect whether major FPC-score changes concentrate around treatment windows.

Do not claim a causal effect without a treatment model. Use this only as interpretive evidence that the extracted trajectories contain clinically meaningful dynamics.

---

# 5. Recommended primary representation for each application

| Application | Raw observation | Encoder | Recommended latent matrix | Meaning of rows | Meaning of columns |
|---|---|---|---|---|---|
| ADNI | 3D T1 MRI visit | BrainDINO | \(S\times D\) | ordered axial slices | frozen MRI features |
| ADNI alternative | one axial 2D slice | BrainDINO | \(P\times D\) | spatial patch tokens | frozen MRI features |
| OLIVES | central OCT B-scan | RETFound-OCT | \(196\times1024\) | spatial patch tokens | frozen retinal features |
| OLIVES alternative | full OCT volume | RETFound-OCT | \(49\times1024\) | B-scan positions | frozen retinal features |

The first ADNI representation and first OLIVES representation are the recommended starting points.

---

# 6. Minimum viable pilot before full analysis

Do not begin by processing the complete datasets.

## Phase 0: environment and checkpoint verification

Confirm:

- data access works;
- encoder checkpoint loads;
- deterministic inference works;
- exact token tensor dimensions are understood;
- metadata contains usable longitudinal dates.

## Phase 1: tiny extraction pilot

Use approximately:

```text
ADNI:   10 subjects × several visits
OLIVES: 10 eyes × several visits
```

Verify visually and numerically that latent matrices have identical shapes and finite values.

## Phase 2: representation diagnostics

Expand to approximately 30–50 subjects/eyes.

Compute:

1. within-subject distance versus time lag;
2. between-subject distance;
3. row/column covariance spectra;
4. effective ranks;
5. approximate Kronecker separability;
6. basic latent PVE curves.

### Go/no-go criteria

Proceed if:

- there is measurable temporal signal;
- latent trajectories are not dominated by obvious acquisition artifacts;
- separability is at least approximately plausible;
- effective ranks are not pathological;
- a modest number of FPCs explains interpretable variation, or the scores show useful clinical associations.

If these fail, change representation layer/token type before changing the FPCA model.

## Phase 3: full FPCA

Fit the proposed estimator on the selected cohort.

Produce all covariance/FPC/score outputs.

## Phase 4: external validation

Use clinical variables and/or held-out visits.

## Phase 5: comparison baselines

At minimum compare:

- raw image representation;
- pooled latent vector;
- latent matrix.

---

# 7. Suggested analysis outputs for the paper

A compact real-data section could contain:

### Figure 1: pipeline

```text
longitudinal images
        ↓
frozen medical foundation encoder
        ↓
longitudinal latent matrices
        ↓
sparse matrix-valued covariance / FPCA
        ↓
population modes + subject trajectories
```

### Figure 2: evidence that the latent process is functional

Plot:

- latent distance versus time lag;
- representative subject latent-FPC trajectories.

### Figure 3: fitted covariance structure

Show:

- eigenvalue spectra of \(U\) and \(V\);
- effective ranks;
- estimated temporal covariance \(K(s,t)\);
- leading temporal eigenfunctions.

### Figure 4: clinical interpretation

Show leading subject scores against:

- disease/progression variables for ADNI;
- biomarkers/treatment/clinical response for OLIVES.

### Optional Figure 5: raw versus latent representation

Compare:

- PVE curves;
- effective ranks;
- held-out-visit latent prediction;
- stability under subsampling.

---

# 8. Important statistical cautions

## 8.1 Latent coordinates are not raw anatomical coordinates

The FPCs describe variation in the chosen encoder representation.

Interpret them as:

> dominant modes of longitudinal variation in learned medical-image features,

not automatically as physical deformation fields or raw image FPCs.

## 8.2 Encoder uncertainty is ignored in the first analysis

The current theory conditions on the frozen representation map \(E\).

This is acceptable for the initial application, but it should be stated explicitly.

## 8.3 Kronecker separability is an empirical approximation

A neural encoder does not guarantee separability.

The fitted model should therefore be accompanied by diagnostics and, if possible, sensitivity analysis to mild nonseparability.

## 8.4 More latent dimensions do not necessarily mean more effective information

Report effective ranks, not only nominal dimensions.

## 8.5 Non-Gaussian/dependent latent scores are expected

Transformer features will not generally produce independent or Gaussian FPC scores.

This is not necessarily a problem for the current project and may actually motivate the general fourth-order dependence theory already developed for the matrix-valued FPCA estimator.

---

# 9. Recommended order of execution

If implementing both applications, use this order:

1. **OLIVES + RETFound central-B-scan patch matrix**
   - checkpoint already available;
   - dataset relatively small;
   - easiest test of the full idea.

2. **ADNI + provisional frozen encoder**
   - build the MRI preprocessing and latent-matrix pipeline;
   - use an available encoder if BrainDINO weights are still inaccessible.

3. **ADNI + BrainDINO**
   - swap in BrainDINO once checkpoint access is available;
   - use the axial-slice × feature matrix as the main representation.

4. Only after these work, explore:
   - different encoder layers;
   - whole-OCT-volume stacking;
   - alternative brain encoders;
   - joint/supervised representation learning.

---

# 10. Final intended story

The intended paper-level story is not:

> deep learning is used to improve FPCA.

It is more specifically:

> Raw longitudinal medical images are often too high-dimensional, nuisance-sensitive, and poorly compressed by separable image-space FPCs. A frozen medical foundation model first maps each scan into a structured latent matrix that preserves spatial/anatomical organization while suppressing pixel-level detail. The proposed sparse matrix-valued functional-data method then models population-level longitudinal covariance and subject-specific variation in this representation space, retaining irregular observation times and allowing non-Gaussian/dependent FPC scores.

This keeps the contribution of the statistical method distinct from the representation-learning model while making the real-data application substantially more realistic.
