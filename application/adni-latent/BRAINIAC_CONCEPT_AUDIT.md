# Does BrainIAC retain global and local brain structure?

## Legacy pre-MatrixFPCA pilot result

This audit used 30 subjects as an early representation check. It is retained
for provenance and is not part of the current BrainIAC + MatrixFPCA evidence,
which uses all 295 eligible subjects. See
[`BRAINIAC_MATFPCA_INTERPRETATION.md`](BRAINIAC_MATFPCA_INTERPRETATION.md) for
the full-cohort analysis.

Yes, but unevenly. BrainIAC retains a mixture of global and local anatomical
information even though it cannot reconstruct a sharp MRI. Ventricular shape
is exceptionally strong at every tested layer. Hippocampal and amygdala signal
is most useful at Block 3 and becomes weaker in later layers. The present audit
does not support a thalamic signal.

This refines the image-reconstruction result: poor pixel reconstruction does
not imply that every local anatomical measurement was discarded. A model can
retain a task-relevant regional statistic without retaining all voxels needed
to draw the region accurately.

No MatrixFPCA is used anywhere in this audit.

## What was tested

For each of the same 133 visits, an independently computed ANTs nonlinear
registration provides a log-Jacobian map on the MNI152NLin2009cAsym 1 mm grid.
A fixed TemplateFlow `aseg` atlas defines seven physical targets:

- broad: whole-brain, cortex, and white-matter mean log-Jacobian;
- local: bilateral hippocampus, lateral ventricles, amygdala, and thalamus.

These are regional deformation measurements, not voxel-reconstruction scores
and not automatically measures of atrophy. The sign depends on registration
direction; this audit only asks whether between-scan variation is encoded.

For each BrainIAC layer, linear kernel-ridge probes use:

- all 216 tokens;
- the mean of all tokens;
- for each local target, up to eight tokens with the greatest atlas overlap;
- an equal number of distant brain tokens as a spatial control.

Subjects—not visits—are split: 19 tune-train, 5 tune-validation, and 6 final
test subjects. Visits are weighted so every subject contributes equally. Ridge
penalties are selected using only the validation subjects. Skill is
`1 - probe MSE / development-mean MSE`; positive skill beats a target-mean
prediction. The fixed-test intervals bootstrap test subjects. Ten additional
nested subject splits check whether the result depends on the six test people.

## Physical results

The table reports median skill over the ten repeated splits. “Positive” is the
fraction of splits in which skill exceeded zero.

| Target and feature scope | Block 3 median (positive) | Block 6 median (positive) | Final median (positive) |
|---|---:|---:|---:|
| Whole brain, all tokens | 0.332 (9/10) | 0.014 (7/10) | 0.020 (7/10) |
| White matter, all tokens | 0.525 (10/10) | 0.307 (9/10) | 0.321 (9/10) |
| Lateral ventricles, all tokens | 0.909 (10/10) | 0.927 (10/10) | 0.911 (10/10) |
| Hippocampus, all tokens | 0.238 (8/10) | 0.096 (6/10) | 0.155 (6/10) |
| Hippocampus, local tokens | 0.421 (9/10) | 0.147 (7/10) | 0.265 (6/10) |
| Amygdala, local tokens | 0.288 (10/10) | 0.184 (7/10) | 0.064 (6/10) |
| Thalamus, local tokens | 0.002 (6/10) | 0.013 (8/10) | 0.011 (6/10) |

The cortex-wide target is omitted from the compact table because median skill
is near zero at every layer. This does not mean cortex is absent; it means this
particular regional mean deformation is not predicted reliably.

Spatial specificity is strongest and most stable at Block 3:

| Local target | Block 3 median ROI-minus-remote skill | ROI better, Block 3 | Block 6 | Final |
|---|---:|---:|---:|---:|
| Hippocampus | 0.362 | 10/10 | 10/10 | 8/10 |
| Lateral ventricles | 0.342 | 10/10 | 10/10 | 10/10 |
| Amygdala | 0.294 | 10/10 | 10/10 | 8/10 |
| Thalamus | -0.015 | 5/10 | 8/10 | 5/10 |

The absolute and relative criteria must be read together. For example, the
final-layer amygdala tokens often beat distant tokens, but their median absolute
skill is only 0.064. That is weak evidence, not a successful local model.

Useful figures and complete tables:

- [atlas regions and selected tokens](concept_audit/anatomical_targets_and_tokens.png)
- [broad physical probes](concept_audit/physical_probe_performance.png)
- [absolute local probe performance](concept_audit/local_probe_absolute_performance.png)
- [ROI-versus-remote spatial specificity](concept_audit/local_spatial_specificity.png)
- [fixed-split metrics](concept_audit/physical_probe_metrics.csv)
- [repeated-split summary](concept_audit/split_stability_summary.csv)
- [repeated spatial-specificity summary](concept_audit/local_specificity_stability.csv)

## Exploratory clinical results

The earliest scan per subject was used, leaving only six final-test subjects.
Age is robustly encoded (all-token skill 0.596, 0.599, and 0.655 for Block 3,
Block 6, and final). On this one split, CDR-SB and FAQ have modest signal at the
middle/final layers, while MMSE is essentially no better than its development
mean. Conversion is not supported: every representation has negative Brier
skill and test AUCs range from 0.44 to 0.67.

These clinical numbers are descriptive only. Six test subjects cannot establish
a clinical association, and BrainIAC pretraining included ADNI. The physical
probes are also measurements derived from the input MRI, so they establish
information retention—not disease causality.

- [clinical plot](concept_audit/clinical_probe_performance.png)
- [clinical metrics](concept_audit/clinical_probe_metrics.csv)

## Consequence for the MatrixFPCA experiment

There is now enough evidence to try longitudinal aggregation, but layer choice
must be treated as part of the scientific design:

1. Use **Block 3 as the local candidate**, because hippocampal and amygdala
   information is strongest and most spatially specific there.
2. Keep the **final layer as the compact/global comparator**. Its near-rank-one
   token structure is convenient for separable MatrixFPCA, but it is less
   convincing for weak local targets.
3. Score MatrixFPCA reconstructions in concept space: fit the same held-out
   regional probes to observed and reconstructed latent trajectories and compare
   regional trajectory error. Latent Frobenius error alone is insufficient.
4. Compare with the existing raw-image ADNI analysis using the same physical
   and clinical targets. A latent method is better only if it improves a target
   of interest without merely reproducing ventricular/global variation.

A two-stream analysis is also scientifically natural: a compact final-layer
global stream plus Block-3 ROI-token streams. It should be evaluated after the
single-layer baselines so extra model flexibility does not receive an unfair
comparison.

## Reproduce

```bash
application/adni-latent/.venv/bin/python \
  application/adni-latent/evaluate_brainiac_concepts.py
```

Use `--stability-repeats=0` for a fast fixed-split run. Outputs include exact
atlas labels, hashes, split counts, ridge tuning, target coverage, and an
explicit `matrixfpca_used: false` provenance field in
[`concept_audit/run_manifest.json`](concept_audit/run_manifest.json).
