# EchoNet-Dynamic MatFPCA application: plan

Status (2026-11-23): preprocessing, Phase-0 checks, RDS conversion and the
pilot fits (56 × 56, n = 500, sparse and dense) are done; all go/no-go checks
pass (see `README.md`). Pilot lesson: use PCA of the MatFPCA-predicted video
(full score tensor) as the EF representation, not the top product-FPC
scores, which are individually tiny (rank-one eigenimages). Main fits run
at 112 × 112 with n = 500 TRAIN subjects (decision 2026-11-23: MatFPCA should
need few subjects when the images are rich; the sample is stratified by EF
quintile × image-appearance cluster and checked for representativeness, and a
56 × 56 sample-size study, n = 250…2,000, checks that 500 is enough). The whole pipeline runs on the cluster (Compute
Canada, R 4.5.0, same conventions as `application/olives/run_retfound.sbatch`).

## 0. Why this dataset

Our ADNI fit (`application/adni/results/temporal_eigencomponent_diagnostics.csv`)
measured time as years since baseline (0–4 y). Its leading temporal mode is
essentially constant (99.0% of the temporal variance), and the only material
dynamic mode carries 0.98%. That invites the question: *why model time at all?*

EchoNet-Dynamic reverses this. Each subject is a cardiac ultrasound video, and the
main signal is the heart contracting over time. The clinical label, ejection
fraction (EF), is by definition a contrast between two time points:

    EF = (EDV − ESV) / EDV,

where EDV is the left-ventricular (LV) volume at end-diastole (ED, the heart
fully filled) and ESV the volume at end-systole (ES, the heart fully
contracted). A static image summary cannot represent EF; a functional
representation should. This is the argument the application has to make.

## 1. What the dataset is

**Source.** EchoNet-Dynamic (Ouyang et al., *Nature* 2020, Stanford).
De-identified echocardiograms from routine clinical care, released for research
use. Each video is the **apical four-chamber (A4C) view**, with the probe at the
heart's apex. It shows all four chambers; the LV appears on the right side of
the image and the apex at the top of the fan-shaped ultrasound sector. The
authors already removed text, ECG and respirometer overlays outside the sector
and resized every frame to 112 × 112. Most of the usual image preprocessing has
therefore been done for us.

**Per video the release provides:**
- `Videos/<id>.avi`: the full clip, usually several heartbeats.
- `FileList.csv`: EF, ESV, EDV, frame size, frame rate (FPS), frame count, and
  split (TRAIN/VAL/TEST).
- `VolumeTracings.csv`: expert LV tracings on exactly two frames, one at ED and
  one at ES. Each tracing is a long-axis line (first row) followed by
  perpendicular chords (later rows). The ED frame is the one with the larger
  traced area.

**Checked locally (`data/EchoNet-Dynamic`, 7.5 GB):**

| Item | Value |
|---|---|
| Videos in `FileList.csv` / on disk | 10,030 / 10,030 |
| Frame size | 10,024 at 112×112; 6 at 768×1024 or 768×1040 (exclude) |
| Split | TRAIN 7,465 · VAL 1,288 · TEST 1,277 |
| FPS | 18 – 138 (median 50) |
| Frames per video | 28 – 1,002 (median 171) |
| EF (%) | 6.9 – 97.0 (median 59.2) |
| Videos with tracings | 10,025, each with exactly 2 traced frames; 1 tracing filename not in `FileList.csv` |
| Order of the traced frames | ED before ES in 9,921; ES before ED in 103 |
| ED→ES gap (seconds) | 5%: 0.22 · median: 0.32 · 95%: 0.44 · max: 4.4 (outliers exist) |
| Decoded pixels | RGB uint8 0–255; channels are **not** identical (a colour overlay is present in some or all frames; check before converting to grey) |

**Visual check** of video `0X100009310A3BD7FC`: a typical A4C fan on a black
background. White tick marks (a depth scale) sit on the right edge of the fan,
and small bright artefacts appear near the valves.

## 2. Statistical formulation

- Subject *i*: one video (one patient).
- Matrix observation Y_i(t) ∈ R^{p×q}: one greyscale frame (p = q = 112, or 56
  for pilots).
- Time t: **cardiac phase**, not clock time. Heart rates and frame rates differ,
  so after registration each subject's frames fall on a subject-specific,
  irregular set of phases. This gives a natural irregular design without
  artificial setup.
- Model: Y_i(t) = μ(t) + Σ ξ_ijkl φ_j ψ_kᵀ θ_l(t) + ε, fitted with the existing
  `fit_matfcov` → `eigen.matfcov` → `trunc.matfeigen` → `get_matfpca_scores`
  chain, exactly as in `application/adni/adni_matfpca.R`.

### 2.1 Time registration (the key preprocessing step)

**Primary: systole only.** Use the two labelled landmarks,

    t = (frame − ED) / (ES − ED) ∈ [0, 1],

keeping frames ED … ES. Use only videos with ED < ES. Drop videos whose ED→ES
gap is under 0.10 s or over 0.80 s; these are probably tracing or labelling
problems (document the count). At 50 fps this keeps about 16 frames per
subject.

Snap each phase to a common grid (step 0.02, 51 points). If two frames of one
subject snap to the same grid point, keep the one closest to it. Record every
original frame index and phase in the manifest.

**Sensitivity: full cycle.** Estimate the next ED as the frame after ES, within
[ES + 0.5·(ES−ED), ED + 2.0 s], that maximises the within-sector correlation
with the ED frame. Then apply piecewise-linear landmark registration:
ED → 0, ES → τ, next ED → 1, where τ is the median systolic fraction across
subjects. Report how often the next-ED estimate is unreliable (for example,
peak correlation < 0.8) and exclude those videos from this variant.

### 2.2 Spatial preprocessing

1. **Grey conversion:** grey = mean of R, G and B, rescaled to [0, 1]. First
   measure colour content per video: the fraction of in-sector pixels where
   max(RGB) − min(RGB) > 20. Report the distribution, and exclude or flag
   videos with substantial colour Doppler.
   *Finding (2026-11-23):* there is no colour Doppler. Colour pixels are
   the static overlays (apex marker at rows 1–5, depth-scale ticks on the
   right fan edge) plus a few videos with a sepia tint over the whole
   frame. The overlays are removed by a per-video fan mask (largest bright
   component, convex hull); the tint is harmless after the RGB mean.
2. **Sector mask:** *(implemented as: pixel inside the per-video fan, after
   any alignment, for ≥ 90% of TRAIN videos.)* Originally: pixels that are
   non-zero in at least 50% of training videos' mean frames (same idea as the ADNI consensus brain mask). Outside the
   mask, set pixels to 0 and leave them out of all metrics. Check whether the
   depth-scale ticks fall inside the mask. They are static, so they can only
   feed constant modes, but this is worth stating.
3. **Intensity normalisation:** ultrasound gain varies between exams.
   Primary: divide each video by its mean in-mask intensity over its retained
   frames (a scale factor only, so it does not distort the dynamics).
   Sensitivity: no normalisation.
4. **Resolution:** 112 × 112 for the main fit. A 2 × 2 block average to
   56 × 56 for pilots and for the dense/sparse comparison grid.
5. **No deformable registration.** *(Update: measured misalignment is
   larger than wall motion: LV centre spreads ~15 px horizontally and ~24 px
   vertically, axis tilt −2° to 27°, axis length 45–75 px (5–95%). The
   alignment is therefore a preprocessing switch, `--align none | rigid |
   similarity`, compared in README.)* The A4C acquisition and the authors'
   cropping already standardise the view.
   Optional sensitivity, LV-centred alignment: apply one similarity transform
   (rotation, scale, translation) per video, taken from the **ED** long-axis
   endpoints, so the LV long axis lands in a canonical position. Using one
   transform per video does not remove wall motion. This should sharpen μ(t)
   and the eigenimages, because between-patient misalignment of the LV is the
   main nuisance.

### 2.3 Sampling design

- **Dense fit:** every retained systolic frame (about 16 per subject).
- **Sparse fit:** for each training subject, sample m_i ~ Uniform{3,…,8}
  frames from that subject's retained frames (fixed seed). This is the regime
  MatFPCA is built for.
- Compare the dense and sparse eigenfunctions (subspace angles), score
  correlations, and EF performance. Stable results show the sparse method
  recovers what the dense data contain.

### 2.4 Subjects

- **Training:** 500 eligible TRAIN videos (was 2,000). They are the first 500
  in a systematic stratified random ranking (EF quintile × 8 appearance
  clusters), so every kind of heart and image is represented in proportion;
  any prefix of the ranking is itself a stratified sample, which the
  sample-size study uses.
- **Test:** all eligible TEST videos (about 1,250), used for held-out score
  prediction, reconstruction and EF evaluation. Leave VAL unused, or use it
  for tuning (for example, the ridge penalty in the EF prediction).
- **Mean function:** a natural-spline mean over phase (df 5–7), fitted on
  training images only with subject-balanced weights. The ADNI
  `fit_mean_function` can be reused as is.

## 3. Analyses and figures

Every figure should be understandable by someone who has not seen this
dataset; README text explains each one.

1. **Dataset explainer:** ED and ES frames from 3 patients (low, normal and
   high EF) with the LV tracing overlaid. Below them, the traced-area or
   in-LV intensity curve over the clip, with ED and ES marked. A registration
   schematic shows clock time on the x-axis mapped to phase.
2. **Mean dynamics μ(t):** frames at t = 0, 0.25, 0.5, 0.75 and 1, plus the
   difference μ(1) − μ(0), which shows average LV contraction.
3. **Temporal diagnostics:** reuse the ADNI constant-versus-dynamic screen and
   show this dataset next to ADNI's 99.0% / 0.98% split. **This is the
   headline.**
4. **Leading product FPCs:** spatial eigenimages φ_j ψ_kᵀ with their θ_l(t),
   and modes of variation μ(t) ± 2√λ · eigencomponent shown as frame strips
   (and a GIF for slides).
5. **EF interpretation:**
   - Fit a principal curve through the scores (reuse
     `application/adni/adni_continuum_utils.R`; `princurve` is not required)
     and colour it by EF; report Spearman correlations with EF, ESV and EDV.
   - Predict EF (and ESV, EDV) with ridge regression; fit on training scores,
     evaluate on TEST. Report R² and MAE, and AUC for EF < 40%
     (reduced-ejection-fraction heart failure).
   - Baselines on the same subjects:
     1. PCA of the ED frame only (static image).
     2. PCA of the average of each subject's observed frames (time
        discarded).
     3. PCA of the stacked dense systolic video (vectorised; time kept but
        not modelled functionally).
     4. MatFPCA dense versus sparse.

     The expected pattern is that static baselines do poorly on EF and
     MatFPCA does well. If the static ED-frame PCA predicts EF almost as
     well, EF is mostly carried by heart size (EDV) and the story must be
     reframed. Report that honestly.
6. **Reconstruction:** for held-out TEST subjects, predict the full systolic
   sequence from 3–8 observed frames. Compute in-mask RMSE and correlation
   against: mean only, mean plus static PCA, and dense versus sparse MatFPCA.
   Show qualitative strips of truth versus reconstruction for 3 subjects.

## 4. Go / no-go checks (cheap, run first)

- **Phase 0 (minutes):** for about 300 videos, decompose in-mask pixel
  variance on the registered systolic grid into between-subject variance
  (subject means) and within-subject temporal variance. If the within-subject
  share is not clearly larger than in ADNI, stop and rethink registration or
  representation before fitting.
- **After the pilot fit (56 × 56, n = 500, sparse):**
  - dynamic temporal modes carry at least 10% of the temporal variance;
  - EF test R² from MatFPCA scores is clearly above that of the static
    ED-frame PCA;
  - the eigenimages are not dominated by fan edges, scale ticks or other
    acquisition artefacts.

## 5. Files to create

```
application/echonet/
├── PLAN.md                    # this file
├── README.md                  # dataset explainer, design, commands, results (after the runs)
├── requirements.txt           # av, numpy, pillow, scipy (Python preprocessing)
├── prepare_echonet.py         # decode, QC, register, mask, normalise, write cache + manifest
├── echonet_matfpca.R          # mean, MatFPCA fit (dense/sparse), scores, reconstruction
├── analyze_echonet.R          # EF prediction, baselines, principal curve, figures
└── run_echonet.sbatch         # cluster pipeline
```

Outputs go to `application/echonet/output/`; add this to `.gitignore` next to
`application/olives/output/`.

**Preprocessing output:**
- `manifest.csv`: one row per retained frame, with video id, split, EF, ESV,
  EDV, FPS, ED, ES, original frame index, phase, grid index, colour fraction,
  normalisation scale and exclusion reason.
- `qc_summary.csv`: counts excluded at each step.
- One uint8 `.npz` per video holding the registered frames.
- An R step that converts the frames into per-subject `rTensor` RDS files.
  `create_matrix_funcdata_irreg()` already accepts a vector of RDS paths
  (disk-backed), which keeps memory low.

Decode videos with **PyAV** (its wheel bundles ffmpeg; no system ffmpeg or
OpenCV needed). Decoding 10k short 112×112 clips is CPU-bound and embarrassingly
parallel: use `multiprocessing` over videos.

## 6. Cluster resources (first guesses; adjust after the pilot)

| Job | Resources | Notes |
|---|---|---|
| Phase 0 variance check | 4 CPU, 8 GB, 20 min | about 300 videos |
| Preprocess all eligible videos | 16 CPU, 16 GB, 1 h | PyAV decoding, writes .npz |
| Pilot fit (56², n = 500, sparse) | 8 CPU, 16 GB, 1 h | go/no-go checks |
| Main fit (112², n = 2,000, dense) | 16 CPU, 64 GB, 6–12 h | dense array ≈ 2000 × 16 × 12544 × 8 B ≈ 3.2 GB, plus covariance work |
| Main fit, sparse | 16 CPU, 32 GB, 4–8 h | |
| Analysis and figures | 8 CPU, 32 GB, 2 h | |

No GPU is needed unless the optional LV-segmentation variant (Section 7) is
added.

**Environment:**
- R 4.5.0 (`module load r/4.5.0`) with rTensor, foreach, mgcv, pracma,
  splines, data.table, png, uwot.
- A Python venv with `application/echonet/requirements.txt`.

**Data transfer:** `rsync -a data/EchoNet-Dynamic/` to the cluster
project/scratch space (7.5 GB). Only `Videos/`, `FileList.csv` and
`VolumeTracings.csv` are needed.

## 7. Optional extensions (only if the main analysis works)

- **LV masks:** replace raw pixels with LV segmentation masks or
  signed-distance maps from the public EchoNet-Dynamic segmentation model.
  This gives a cleaner LV shape signal at the cost of a learned preprocessing
  step (GPU).
- **Full-cycle registration:** the Section 2.1 sensitivity variant, which
  covers both systole and diastole.
- **ADNI re-run with age as the time axis:** a separate, cheap check of the
  AABC idea (accelerated design over ages 36–90) with data already on disk.

## 8. Known risks

- **Spatial misalignment:** between-patient differences in LV position and
  scale may dominate the constant modes. The LV-centred alignment sensitivity
  addresses this.
- **EF and heart size:** EF may correlate with static anatomy (dilated LV ↔
  low EF). The static baselines are there to measure this.
- **Landmark quality:** the two tracing frames are expert-selected but may
  not be exactly the true ED/ES; 103 videos label ES before ED, and gap
  outliers exist.
- **Colour overlays and static artefacts:** colour content, depth ticks and
  valve speckle could produce artefact eigenimages. Check the colour content
  first.
- **Weak temporal identifiability:** with roughly 16 frames per subject over
  systole, the temporal covariance is well covered on [0, 1]. This is not the
  "functional snippet" problem that AABC would have.
