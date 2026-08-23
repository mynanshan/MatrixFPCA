# Full-cohort BrainIAC + MatrixFPCA interpretation

## Scope

This is the main 295-subject, 1,308-visit ADNI analysis through month 48. It
models longitudinal trajectories in frozen BrainIAC feature space. It does not
decode BrainIAC features or claim MRI reconstruction.

The same subject-visit keys are used for Block 3, Block 6, and final-normalized
features. Population feature PCA, temporal means, and MatrixFPCA covariance are
fitted using 236 subjects and 1,033 visits. Fifty-nine subjects are excluded
from population fitting; 58 contribute one held-out interior visit while their
scores are estimated from their remaining visits.

## Revised preprocessing and retention

Every registered `193 × 229 × 193` MRI receives the same crop
`[16:177, 16:217, 0:177]`, followed by resize to `96³`, nonzero-voxel
normalization, and frozen BrainIAC encoding. A common crop improves brain
coverage without destroying token-to-anatomy correspondence.

All three encoder outputs have shape `216 × 768`. Training-only feature PCA
retains the smallest rank reaching 98% variance. MatrixFPCA then retains every
numerically positive token and feature covariance mode; only temporal modes are
truncated at 98% PVE, with a maximum rank of four.

| Layer | Feature-PCA rank and PVE | Positive token × feature × time ranks | Product modes | Separability correlation |
|---|---:|---:|---:|---:|
| Block 3 | 83; 98.01% | `144 × 73 × 2` | 21,024 | 0.138 |
| Block 6 | 96; 98.01% | `216 × 96 × 3` | 62,208 | 0.901 |
| Final normalized | 42; 98.01% | `7 × 9 × 2` | 126 | 0.995 |

The final layer remains intrinsically concentrated relative to the earlier
layers, but it is no longer reduced to the previous three product scores. The
126 retained modes are all modes supported by positive estimated marginal
covariance eigenvalues under the revised fit.

`subject_fpc_scores.csv` contains the leading 100 ordered product scores for
each layer, giving a `295 × 100` score matrix before metadata columns. The full
positive-mode tensor remains in the fitted model.

## Individual morphology coordinates

Each product FPC is learned without phenotype labels. Selection of one mode per
layer is post hoc and based on regional morphology association among the leading
20 reliable scores. Signs are oriented so positive values correspond to larger
ventricular morphology levels.

| Layer | Selected mode | Mode variance | Ventricular rho | Bootstrap 95% interval | Partial rank R² | Global adjusted q |
|---|---|---:|---:|---:|---:|---:|
| Block 3 | FPC7 = R2 × C2 × T1 | 0.51% | 0.854 | 0.815–0.885 | 0.673 | 1.4e-69 |
| Block 6 | FPC1 = R1 × C1 × T1 | 24.80% | 0.700 | 0.620–0.766 | 0.417 | 4.2e-33 |
| Final normalized | FPC1 = R1 × C1 × T1 | 70.07% | 0.623 | 0.533–0.701 | 0.313 | 1.1e-22 |

Partial rank models adjust for baseline age, sex, visit count, and observed time
span. Global q-values adjust across all tested FPC-phenotype associations in all
three layers.

Block 3 FPC7 is also spatially localized: its ventricular-token to matched
remote-token mean absolute-loading ratio is 44.3. Block 6 has a weaker
ventricular-matched localized mode at FPC10 (ratio 2.66). Some final-layer row
modes are localized, but the strongest regions do not match their scalar
phenotype associations and should not be presented as region-specific findings.

Primary files:

- [`selected_morphology_coordinates.csv`](results/brainiac_matfpca_interpretation_full/selected_morphology_coordinates.csv)
- [`continuous_morphology_coordinates.png`](results/brainiac_matfpca_interpretation_full/continuous_morphology_coordinates.png)
- [`localized_eigenmode_candidates.csv`](results/brainiac_matfpca_interpretation_full/localized_eigenmode_candidates.csv)
- [`selected_eigenmode_token_time_maps.png`](results/brainiac_matfpca_interpretation_full/selected_eigenmode_token_time_maps.png)

## Label-free score continua

Principal curves are fitted from the leading 20 product scores without anatomy,
age, diagnosis, or clinical labels. Their strongest regional associations are:

| Layer | Regional target | Spearman rho | Partial rank R² |
|---|---|---:|---:|
| Block 3 | Lateral-ventricle level | 0.763 | 0.499 |
| Block 6 | Lateral-ventricle level | 0.612 | 0.311 |
| Final normalized | Lateral-ventricle level | 0.454 | 0.147 |

The principal curve explains 56.6%, 46.0%, and 80.2% of the selected score-space
variation for Block 3, Block 6, and final, respectively. Its improvement over a
single linear axis is modest, so these continua are close to one-dimensional
linear orderings rather than strongly curved manifolds.

Observed-MRI panels show actual baseline scans from subjects nearest continuum
percentiles. They validate anatomical ordering but are not generated images or
feature inversions.

Primary files:

- [`latent_continuum_top_morphology.csv`](results/brainiac_matfpca_interpretation_full/latent_continuum_top_morphology.csv)
- [`latent_continuum_diagnostics.csv`](results/brainiac_matfpca_interpretation_full/latent_continuum_diagnostics.csv)
- [`latent_continuum_interpretation.png`](results/brainiac_matfpca_interpretation_full/latent_continuum_interpretation.png)
- [`raw_vs_brainiac_observed_mri_continua.png`](results/brainiac_matfpca_interpretation_full/raw_vs_brainiac_observed_mri_continua.png)

## Interpretation limits

- BrainIAC was not trained as an image autoencoder and has no compatible
  pretrained decoder. Feature-space fit does not establish voxel reconstruction.
- The common crop changes the encoder input distribution relative to the
  checkpoint's uncropped validation preprocessing. The resulting representations
  are empirically validated here but should not be described as an official
  BrainIAC preprocessing recipe.
- The very high final-layer separability and low positive marginal rank indicate
  strong global concentration. This can be useful for dominant morphology but
  may suppress lower-variance personalized structure.
- Block 3 has much richer spatial rank and the strongest localized ventricular
  result, but its separability correlation is only 0.138. Product-mode inference
  for that layer therefore relies on a substantially weaker Kronecker
  approximation than Block 6 or final.
- Phenotype-based mode selection is exploratory even though the FPCs themselves
  are label-free.

## Reproduction

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
