# Full-cohort BrainIAC + MatrixFPCA interpretation

## Bottom line

The full analysis uses all 295 eligible subjects and 1,308 MRI visits through
month 48. It replaces the earlier 30-subject pilot, whose temporal ranks,
correlations, and nonlinear plots were too unstable for scientific use.

The full results are much clearer:

- Block 3 contains a strongly localized ventricular eigenmode and the strongest
  single morphology coordinate.
- Final-layer features provide a compact global ventricular/atrophy continuum.
- Block 6 contains useful anatomy-linked individual modes, but its unsupervised
  principal curve is not aligned with a clear physical continuum.
- Multiple regional longitudinal-slope associations now survive strict
  layer-wide multiplicity correction.
- Actual observed MRIs ordered along the Block 3 and final continua show
  visible ventricular enlargement and global atrophy, broadly resembling the
  raw-image MatrixFPCA result.

No image reconstruction is used. Observed-image panels retrieve the actual
baseline MRI of subjects nearest six continuum percentiles.

## Cohort and validation

- 295 subjects and 1,308 visits, with at least three retained visits per subject.
- 236 population-training subjects and 1,033 visits fit feature PCA, the mean,
  and covariance.
- 59 entirely separate validation subjects; 58 have one interior visit held
  out while their scores are estimated from their remaining visits.
- The same subject-image keys are used for Block 3, Block 6, and final features.
- Feature PCA, temporal mean fitting, and subject splitting remain
  population-training-only.

The full retained ranks are:

| Layer | Token x feature x time ranks | Product modes | Feature-PCA PVE | Separability correlation |
|---|---:|---:|---:|---:|
| Block 3 | 24 x 24 x 1 | 576 | 96.91% | 0.250 |
| Block 6 | 3 x 7 x 2 | 42 | 96.95% | 0.950 |
| Final | 1 x 3 x 1 | 3 | 99.22% | 0.995 |

The pilot had temporal rank two for every layer. In the full cohort, Block 3
and final retain only one temporal mode, demonstrating that the pilot's second
temporal mode was unstable.

## 1. Individual continuous morphology coordinates

Each displayed FPC is learned without phenotype labels; the choice of which FPC
to display is post hoc and based on its physical association. Signs are oriented
so larger values correspond to larger ventricular level.

| Layer | Selected product mode | Mode variance | Ventricular rho | Bootstrap 95% interval | Partial rank R2 | Strict adjusted q |
|---|---|---:|---:|---:|---:|---:|
| Block 3 | FPC15 = R2 x C3 x T1 | 1.00% | 0.811 | 0.760--0.853 | 0.582 | 4.9e-54 |
| Block 6 | FPC1 = R1 x C1 x T1 | 52.83% | 0.421 | 0.311--0.526 | 0.133 | 6.2e-8 |
| Final | FPC1 = R1 x C1 x T1 | 77.91% | 0.687 | 0.614--0.750 | 0.391 | 3.8e-31 |

These associations remain after rank-based adjustment for age, sex, visit
count, and observed time span, and after BH adjustment across every tested FPC
and phenotype in the layer.

See
[`continuous_morphology_coordinates.png`](results/brainiac_matfpca_interpretation_full/continuous_morphology_coordinates.png)
and
[`selected_morphology_coordinates.csv`](results/brainiac_matfpca_interpretation_full/selected_morphology_coordinates.csv).

## 2. Spatially localized eigenmodes

Block 3 FPC15 is not merely globally correlated with ventricular size. Its
token loading has a ventricular ROI-to-matched-remote ratio of 26.8. It is thus
the clearest candidate for a genuinely local anatomical MatrixFPCA mode.
Block 3 FPC13 uses the same localized row mode and also strongly tracks
ventricular level (rho = -0.712; strict adjusted q = 4.3e-35).

Block 6 has localized row modes, but their strongest loading regions are
amygdala or thalamus while their best scalar association is ventricular level.
Those modes should not be described as region-matched local discoveries.
The final fit has token rank one and cannot express a localized token contrast.

See
[`localized_eigenmode_candidates.png`](results/brainiac_matfpca_interpretation_full/localized_eigenmode_candidates.png)
and
[`localized_eigenmode_candidates.csv`](results/brainiac_matfpca_interpretation_full/localized_eigenmode_candidates.csv).

## 3. Unsupervised morphology continua

Principal curves were fit without anatomy, age, diagnosis, or clinical labels.

| Layer | Score dimensions | Score variation explained | Ventricular-level rho | Adjusted partial rank R2 | Adjusted q |
|---|---:|---:|---:|---:|---:|
| Block 3 | 20 | 50.6% | 0.573 | 0.206 | 7.1e-15 |
| Block 6 | 20 | 61.7% | not selected | -- | -- |
| Final | 3 | 77.3% | 0.669 | 0.370 | 1.8e-29 |

Block 3 and final therefore provide statistically and visually interpretable
global morphology continua. Their coordinates have rho = 0.680 across subjects.
Block 6's continuum agrees weakly with either layer and its strongest regional
association does not survive within-layer FDR, even though individual Block 6
FPCs retain significant anatomy.

Observed-image landmark rows make this distinction visible. The Block 3 and
final rows progress toward larger ventricles and greater atrophy; the Block 6
row is less consistently ordered. Every panel is an original MRI, not a decoded
or nearest-neighbor-generated image.

See
[`raw_vs_brainiac_observed_mri_continua.png`](results/brainiac_matfpca_interpretation_full/raw_vs_brainiac_observed_mri_continua.png),
[`latent_continuum_interpretation.png`](results/brainiac_matfpca_interpretation_full/latent_continuum_interpretation.png),
and
[`latent_continuum_diagnostics.csv`](results/brainiac_matfpca_interpretation_full/latent_continuum_diagnostics.csv).

## 4. Longitudinal morphology

The full cohort supports associations that were invisible in the pilot:

- Block 3 FPC15 versus hippocampal annual slope: rho = 0.422, partial rank
  R2 = 0.167, strict adjusted q = 5.2e-11.
- Block 6 FPC4 versus ventricular annual slope: rho = 0.249, partial rank
  R2 = 0.087, strict adjusted q = 2.4e-5.
- Final FPC1 versus hippocampal annual slope: rho = -0.406, partial rank
  R2 = 0.145, strict adjusted q = 5.4e-10.
- Final FPC3 versus ventricular annual slope: rho = 0.253, partial rank
  R2 = 0.082, strict adjusted q = 1.4e-5.

These show that longitudinal BrainIAC + MatrixFPCA scores correspond to
independently measured regional deformation rates. A temporal component index
must still not be read naively as "level" or "slope"; the interpretation comes
from the external deformation association.

## 5. Age, cognition, and conversion

The Block 3 continuum correlates with age (rho = 0.349, adjusted q = 3.2e-7),
and the final continuum has rho = 0.263 (adjusted q = 3.5e-4). These are weaker
than their ventricular associations, so age does not define either continuum.
The raw-image age correlation is 0.432 in 393 subjects.

The full continua and individual modes also associate with MMSE, CDR-SB, FAQ,
sex, and incident conversion after adjustment. For example, continuum-to-
conversion rho is 0.216 for Block 3 and 0.250 for final, with adjusted q below
0.01. APOE4 is not significant after strict adjustment.

These are associations, not prospective prediction results. The conversion
labels are not held out as an external test outcome, and BrainIAC may have seen
ADNI during pretraining. Clinical claims require a separate cohort-level
validation design.

See
[`raw_vs_brainiac_latent_coordinate_vs_age.png`](results/brainiac_matfpca_interpretation_full/raw_vs_brainiac_latent_coordinate_vs_age.png)
and
[`eigenmode_phenotype_associations.csv`](results/brainiac_matfpca_interpretation_full/eigenmode_phenotype_associations.csv).

## Comparison with raw-image MatrixFPCA

The full BrainIAC result now has the same three forms of evidence as the
raw-image analysis:

1. a continuous score-space ordering;
2. PCA/UMAP/t-SNE embeddings colored by that ordering;
3. actual observed MRI landmarks along the coordinate.

Block 3 and final recover a recognizable global ventricular/atrophy sequence.
Block 3 additionally supplies localized ventricular eigenmodes that raw-image
variance components did not emphasize as directly. Final is much more compact
but sacrifices local token contrast. This global/local tradeoff is the main
scientific reason to retain both representations in the next comparison.

## Reproduce and inspect

Run the fitted-model interpretation with:

```bash
Rscript application/adni-latent/interpret_brainiac_matfpca.R
```

The script requires at least 250 subjects by default, so the 30-subject pilot
cannot silently become the reported analysis. Outputs are written to
`application/adni-latent/results/brainiac_matfpca_interpretation_full`.

Important machine-readable files include:

- [`eigenmode_phenotype_associations.csv`](results/brainiac_matfpca_interpretation_full/eigenmode_phenotype_associations.csv);
- [`latent_continuum_associations.csv`](results/brainiac_matfpca_interpretation_full/latent_continuum_associations.csv);
- [`latent_continuum_observed_mri_landmarks.csv`](results/brainiac_matfpca_interpretation_full/latent_continuum_observed_mri_landmarks.csv);
- [`latent_score_embeddings.csv`](results/brainiac_matfpca_interpretation_full/latent_score_embeddings.csv);
- [`token_roi_localization.csv`](results/brainiac_matfpca_interpretation_full/token_roi_localization.csv);
- [`cross_layer_score_concordance.csv`](results/brainiac_matfpca_interpretation_full/cross_layer_score_concordance.csv).

The legacy 30-subject directories are retained only for provenance and are not
used by the default full-cohort scripts or current conclusions.
