# Frozen RETFound feature analysis

Completed 2026-09-05 on the local TREX-DME release. The cohort contains 47
patients, one study eye per patient, and 785 eligible eye-visits through V27.
Each visit is represented by the nominal central OCT B-scan (slice 24).

The official OCT RETFound MAE ViT-L/16 encoder was loaded frozen from the local
OCT checkpoint. RETFound is a masked-autoencoder vision transformer: its
representation was learned by reconstructing masked retinal-image patches. The
ViT-Large encoder has 24 transformer blocks, 16 × 16 patches, and a 1,024-
dimensional hidden state. Evaluation preprocessing was RGB conversion, bicubic
resize of the short side to 256 pixels, center crop to 224 × 224, and ImageNet
normalization. The final normalized encoder tokens excluding CLS form a
196 × 1024 patch-by-feature matrix. No pixel alignment, downsampling, outcome
fitting, or random augmentation was used. Repeat inference differed by zero
at the saved precision, and the implementation matched the official encoder
sequence exactly.

Each row of a saved matrix is one of the 14 × 14 spatial image patches. Each
column is a learned latent channel; individual channels are not expected to be
human-readable biomarkers. The separate pooled vector averages the patch
tokens. MatFPCA receives the full matrix, preserving patch order and feature
coordinates, rather than receiving raw pixels or only the pooled vector.

The saved extraction metadata are in
`output/retfound/features/extraction.json`; the manifest is
`output/retfound/manifest.csv`. The matrix-valued fit is in
`output/retfound/matfpca/model.rds`. It retained all 196 patch and 1,024
feature directions and three temporal components under the 95% temporal-PVE
rule. The covariance fit took 162 seconds locally.

As a representation screen, PCA of the pooled patch tokens showed observable
clinical signal. Component 3 had arm-adjusted r = −0.45 with baseline BCVA
(BH q = 0.016) and r = 0.51 with baseline CST (BH q = 0.0096). For the 31
patients with V13, component 4 had arm-adjusted r = 0.57 with V13−V1 CST
change (BH q = 0.015); component 3 showed a similar association. These tests
are descriptive and use the full longitudinal feature extraction, so they are
not prospective prediction results.

The pooled RETFound screen and the individual MatFPCA score screen answer
different questions. The pooled screen first summarizes each baseline feature
vector with ordinary PCA, then tests those unsupervised coordinates. The
individual MatFPCA screen tests one separable product score at a time. A
feature representation can contain clinical information spread across many
weak product scores even when no single score is convincing.

In this fit, the leading individual product scores did not show a surviving
arm-adjusted association after correction; the smallest reported q-values were
above 0.58. The joint score analysis is therefore the more appropriate next
test.

## Joint scores and morphology coordinate

The complete retained MatFPCA score tensor was converted to an exact Euclidean
coordinate using its patient-by-patient Gram matrix. This preserves all
pairwise distances in the full product-score vector while reducing the numeric
dimension to at most 46, the maximum rank available for 47 patients. Joint
clinical tests use a full-geometry distance statistic with 4,999 arm-blocked
residual permutations; a first-ten-coordinate linear regression is retained as
a sensitivity analysis.

The joint score vector was screened against enrollment age, baseline BCVA,
baseline CST, V13−V1 BCVA/CST changes, and all 16 available baseline biomarker
fields, plus available numeric enrollment and treatment quantities from the
full DME workbook. The biomarker export supplies baseline labels for all 47
selected patients but does not provide V13 labels for this cohort, so biomarker
tests are baseline-only. The primary full-geometry distance test uses all 46
nonzero score coordinates and found no BH-significant result across the
37-outcome family. Its nominal strongest signals were Fluid (IRF) (p = 0.0254;
q = 0.612), height (p = 0.0844; q = 0.612), and preretinal
tissue/hemorrhage (p = 0.0524; q = 0.612). The secondary first-ten-coordinate
linear test gives similar nominal signals but is not the primary evidence.
The principal curve through the complete score geometry explained 54.6% of
score variation, compared with 35.9% for a straight axis. Its coordinate had no
corrected association with any tested quantity. Nominal coordinate associations
included Fluid (IRF) (arm-adjusted p = 0.0146; q = 0.423) and DRT/ME (p =
0.0396; q = 0.423) after correcting over the expanded 37-outcome
coordinate-screen family. Age was only one of the tested quantities
(Spearman rho = 0.03; arm-adjusted p = 0.742). Thus the curve is a useful
geometric summary to inspect, but it is not yet an interpretable clinical
morphology axis in this cohort.

The curve and the primary distance test use all 46 nonzero Gram coordinates.
The first-ten-coordinate regression is retained only as a rank sensitivity
analysis because 47 patients cannot support an unrestricted 46-coordinate
linear regression with arm adjustment. The MatFPCA fit itself retains all
196 × 1,024 spatial/product directions and three temporal components.
Downstream truncation can discard distributed signal and should be treated as a
sensitivity choice, not as a harmless preprocessing detail.

## PCA modes of the complete MatFPCA score representation

As a complementary linear analysis, ordinary PCA was applied to the complete
patient-by-patient MatFPCA score representation after the exact Gram-coordinate
embedding. This produces 46 nonzero score modes; no MatFPCA product score or
PCA mode was discarded before the clinical screen. The 46 modes were tested
against the same 37-outcome panel, giving 1,702 mode/outcome tests. The output
contains both a global BH correction over all 1,702 tests and a within-outcome
correction over the 46 modes.

The unfiltered audit table contains one apparently strong result: score-PCA
mode 19 explains only 0.899% of score variance and has an arm-adjusted
correlation of −0.509 with baseline OCT **Atrophy/thinning** (p = 0.000256;
within-outcome BH q = 0.0118; global q = 0.342). This is not evidence of a
meaningful clinical mode in this cohort: the label is one positive patient and
46 negative patients. The correlation is therefore vulnerable to single-point
leverage and is retained only to make the sparse-label failure visible.

After excluding constant labels, labels with fewer than five observations in
one level, and outcomes without arm variation, 25 outcomes remain (1,150
mode/outcome tests). No mode has global BH q < 0.05 in this screenable family.
The strongest
within-outcome result is mode 2 or 4 versus baseline preretinal
tissue/hemorrhage (q = 0.084 after correcting over modes); mode 39 versus
V13−V1 CST change has r = 0.517, p = 0.00288, q = 0.132. Neither is corrected
across the full clinical panel. Thus the extra PCA does not currently deliver
a robust, multiplicity-safe clinical mode, although the audit still shows why
variance truncation must be treated as a sensitivity choice.

## A promising representation and framing direction

The saved RETFound CLS token gives a useful diagnostic contrast with the patch
matrix. I fitted an outcome-blind PCA to the 47 baseline CLS vectors only,
then tested the resulting modes against the same 37-outcome panel, including
baseline BCVA/CST and V13−V1 changes. CLS mode 2 explains 11.49% of baseline
embedding variance and has arm-adjusted r = −0.576 with baseline CST
(p = 2.25 × 10⁻⁵; global BH q = 0.0196 across all 1,110 tests; q = 0.0169
after excluding sparse/constant labels). Its rank correlation is also strong
(rho = −0.560, p = 4.21 × 10⁻⁵). This is the strongest association in the
broader screen that remains significant after testing a prespecified clinical
panel and has a continuous, clinically interpretable endpoint.

The same baseline CLS mode is associated with V13−V1 CST change (r = 0.631,
p = 1.40 × 10⁻⁴; within-outcome q = 0.0042), but its broad-panel q = 0.0611
is just above 0.05. That result is therefore a promising prognostic hypothesis,
not a confirmed treatment-response marker. Baseline BCVA also has a nominal
association (r = 0.482, p = 0.00061) that does not survive the broad-panel
correction.

This contrast suggests a concrete reframing. The current patch-token
longitudinal MatFPCA is asking one covariance model to represent both static
baseline morphology and temporal change, while its subject BLUP scores are
estimated from every visit. A global CLS direction can therefore be diluted or
rotated away in the full score tensor even when it is clinically real. The next
analysis should separate: (1) a baseline matrix-FPCA/global-embedding axis for
structural morphology and baseline CST, and (2) a change-specific functional
analysis for treatment response. The CLS axis should be retained as a
prespecified covariate when testing longitudinal outcomes, rather than selected
after looking at the outcome. A pilot that routed the CLS trajectories through
the current full longitudinal MatFPCA still did not reproduce the baseline CST
association, so simply swapping the token type is not sufficient; the baseline
versus change decomposition is part of the proposed fix.

The baseline CLS screen is in
`output/retfound/embedding_screen/cls_baseline_pca_clinical_associations.csv`
with its spectrum in `cls_baseline_pca_spectrum.csv` and patient-level scores
in `cls_baseline_pca_scores.csv`. These results do not
replace the patch-matrix MatFPCA result; they identify a representation and
problem-framing direction that should be validated before adding more flexible
supervised models.

Joint tests are in `output/retfound/matfpca/joint_score_omnibus_tests.csv`;
curve coordinates and clinical variables are in
`output/retfound/matfpca/morphology_coordinates_clinical.csv`; curve
diagnostics are in `output/retfound/matfpca/curve_diagnostics.csv`; PCA mode
associations are in `output/retfound/matfpca/pca_mode_clinical_associations.csv`
and the mode spectrum is in `output/retfound/matfpca/score_pca_spectrum.csv`.

## Multiple-testing correction

BH means the Benjamini–Hochberg procedure. Given many hypothesis-test p-values,
it ranks them and adjusts each according to its rank and the number of tests.
The resulting q-value estimates the expected false-discovery proportion when
calling tests significant at a chosen threshold. Here, q-values are used to
avoid treating an isolated small p-value among many score/outcome tests as
strong evidence. A q-value below 0.05 means the result passes that analysis
family's false-discovery threshold; it does not make the result confirmatory or
remove dependence on the cohort, representation, or analysis choices.

Files produced by the reusable workflow:

- `prepare_retfound.py`: cohort and provenance manifest;
- `prepare_retfound_demographics.py`: reusable enrollment and treatment quantities;
- `prepare_retfound_biomarkers.py`: visit-level join to sparse OCT biomarkers;
- `extract_retfound.py`: frozen encoder and patch matrices;
- `check_retfound_clinical.py`: pooled-feature clinical screen;
- `screen_retfound_embeddings.py`: baseline CLS/pooled embedding screen against
  the broad clinical panel;
- `matfpca_retfound.R`: feature-matrix MatFPCA fit;
- `associate_retfound_scores.R`: leading product-score association screen;
- `continuum_retfound.R`: joint-score permutations and principal curve;
- `pca_retfound_scores.R`: PCA of the complete MatFPCA score geometry and
  clinical association screen;
- `run_retfound.sbatch`: GPU workflow for the supplied environment.
