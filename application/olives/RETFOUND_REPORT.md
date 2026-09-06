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
mode/outcome tests). No
mode has global BH q < 0.05 in this screenable family. The strongest
within-outcome result is mode 2 or 4 versus baseline preretinal
tissue/hemorrhage (q = 0.084 after correcting over modes); mode 39 versus
V13−V1 CST change has r = 0.517, p = 0.00288, q = 0.132. Neither is corrected
across the full clinical panel. Thus the extra PCA does not currently deliver
a robust, multiplicity-safe clinical mode, although the audit still shows why
variance truncation must be treated as a sensitivity choice.

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
- `matfpca_retfound.R`: feature-matrix MatFPCA fit;
- `associate_retfound_scores.R`: leading product-score association screen;
- `continuum_retfound.R`: joint-score permutations and principal curve;
- `pca_retfound_scores.R`: PCA of the complete MatFPCA score geometry and
  clinical association screen;
- `run_retfound.sbatch`: GPU workflow for the supplied environment.
