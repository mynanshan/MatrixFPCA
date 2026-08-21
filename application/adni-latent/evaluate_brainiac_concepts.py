#!/usr/bin/env python3
"""Audit global and local anatomy in frozen BrainIAC tokens, before MatrixFPCA.

The targets are regional mean log-Jacobians from an independently computed
ANTs deformation field on the same MNI template grid.  A fixed aseg atlas
defines the regions.  Every probe is tuned and evaluated with subjects—not
visits—held out.  Local probes compare tokens overlapping an ROI with the same
number of distant brain tokens, which distinguishes spatially specific signal
from a global subject/disease signal.

This is a representation audit, not a causal or clinical validation study.
The default pilot has only 30 subjects; clinical results are exploratory.
MatrixFPCA is never imported by this script.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any, Iterable, Sequence

import numpy as np

import evaluate_brainiac_reconstruction as reconstruction
import extract_brainiac_latents as extractor


SCRIPT_VERSION = "1.0.0"
DEFAULT_SEED = 20260819
DEFAULT_OUTPUT_DIR = "application/adni-latent/concept_audit"
DEFAULT_ATLAS = (
    "application/adni-latent/resources/templateflow/"
    "tpl-MNI152NLin2009cAsym_res-01_seg-aseg_dseg.nii.gz"
)
DEFAULT_NONLINEAR_DIR = (
    "/home/muye/Projects/SurvADNI/ADNI_raw_data/downloaded_images/"
    "preprocessed/nonlinear"
)
DEFAULT_ALPHA_RATIOS = reconstruction.DEFAULT_ALPHA_RATIOS
BOOTSTRAP_REPLICATES = 2000
MAX_LOCAL_TOKENS = 8
DEFAULT_STABILITY_REPEATS = 10

# FreeSurfer aseg integer labels in the TemplateFlow segmentation.
ROI_LABELS: dict[str, tuple[int, ...]] = {
    "whole_brain": tuple(range(1, 300)),
    "cortex": (3, 42),
    "white_matter": (2, 41),
    "hippocampus": (17, 53),
    "lateral_ventricles": (4, 5, 43, 44),
    "amygdala": (18, 54),
    "thalamus": (10, 49),
}
LOCAL_TARGETS = ("hippocampus", "lateral_ventricles", "amygdala", "thalamus")
CLINICAL_TARGETS = {
    "baseline_age": ("continuous", "baseline_age"),
    "baseline_mmse": ("continuous", "baseline_mmse"),
    "baseline_cdrsb": ("continuous", "baseline_cdrsb"),
    "baseline_faq": ("continuous", "baseline_faq"),
    "incident_ad_conversion": ("binary", "incident_ad_conversion"),
}


class ConceptAuditError(RuntimeError):
    """A user-actionable concept-audit failure."""


def log(message: str) -> None:
    print(time.strftime("%Y-%m-%d %H:%M:%S"), "|", message, flush=True)


def finite_float(value: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan
    return result if math.isfinite(result) else math.nan


def parse_binary(value: str) -> float:
    normalized = value.strip().lower()
    if normalized in extractor.TRUE_VALUES:
        return 1.0
    if normalized in extractor.FALSE_VALUES:
        return 0.0 if normalized not in {"", "na", "nan", "null"} else math.nan
    return math.nan


def subject_weights(rows: Sequence[dict[str, str]], indices: np.ndarray) -> np.ndarray:
    counts: dict[str, int] = {}
    for index in indices:
        subject = rows[int(index)]["subject_id"]
        counts[subject] = counts.get(subject, 0) + 1
    weights = np.asarray(
        [1.0 / counts[rows[int(index)]["subject_id"]] for index in indices],
        dtype=np.float64,
    )
    # Mean one is convenient; the relative weighting is what matters.
    return weights * (len(weights) / weights.sum())


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.sum(values * weights[:, None], axis=0) / weights.sum()


def prepare_kernel(
    features: np.ndarray,
    fit_indices: np.ndarray,
    predict_indices: np.ndarray,
    fit_weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    x_fit = np.asarray(features[fit_indices], dtype=np.float32)
    x_predict = np.asarray(features[predict_indices], dtype=np.float32)
    mean = weighted_mean(x_fit, fit_weights).astype(np.float32)
    weighted_centered = (x_fit - mean) * np.sqrt(fit_weights)[:, None].astype(np.float32)
    scale = float(
        np.einsum("ij,ij->", weighted_centered, weighted_centered, dtype=np.float64)
        / fit_weights.sum()
    )
    if not math.isfinite(scale) or scale <= 0:
        raise ConceptAuditError("A probe feature representation has zero variance.")
    kernel = (weighted_centered @ weighted_centered.T).astype(np.float64) / scale
    cross = ((x_predict - mean) @ weighted_centered.T).astype(np.float64) / scale
    return kernel, cross


def fit_from_kernel(
    kernel: np.ndarray,
    cross: np.ndarray,
    fit_targets: np.ndarray,
    fit_weights: np.ndarray,
    alpha: float,
) -> tuple[np.ndarray, np.ndarray]:
    target_mean = weighted_mean(fit_targets, fit_weights)
    right = (fit_targets - target_mean) * np.sqrt(fit_weights)[:, None]
    regularized = kernel.copy()
    regularized.flat[:: regularized.shape[0] + 1] += alpha
    try:
        coefficients = np.linalg.solve(regularized, right)
    except np.linalg.LinAlgError as error:
        raise ConceptAuditError(f"Kernel-ridge solve failed at alpha={alpha:g}: {error}") from error
    return target_mean + cross @ coefficients, target_mean


def weighted_mse(truth: np.ndarray, estimate: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(weights * np.square(estimate - truth)) / weights.sum())


def weighted_correlation(truth: np.ndarray, estimate: np.ndarray, weights: np.ndarray) -> float:
    truth_mean = float(np.sum(weights * truth) / weights.sum())
    estimate_mean = float(np.sum(weights * estimate) / weights.sum())
    truth_centered = truth - truth_mean
    estimate_centered = estimate - estimate_mean
    numerator = float(np.sum(weights * truth_centered * estimate_centered))
    denominator = math.sqrt(
        float(np.sum(weights * truth_centered**2))
        * float(np.sum(weights * estimate_centered**2))
    )
    return numerator / denominator if denominator > 0 else math.nan


def auc(truth: np.ndarray, score: np.ndarray) -> float:
    positive = score[truth == 1]
    negative = score[truth == 0]
    if len(positive) == 0 or len(negative) == 0:
        return math.nan
    comparisons = positive[:, None] - negative[None, :]
    return float((np.sum(comparisons > 0) + 0.5 * np.sum(comparisons == 0)) / comparisons.size)


def subject_mse(
    rows: Sequence[dict[str, str]], indices: np.ndarray, truth: np.ndarray, estimate: np.ndarray
) -> tuple[list[str], np.ndarray]:
    grouped: dict[str, list[float]] = {}
    for position, index in enumerate(indices):
        grouped.setdefault(rows[int(index)]["subject_id"], []).append(
            float((estimate[position] - truth[position]) ** 2)
        )
    subjects = sorted(grouped)
    return subjects, np.asarray([np.mean(grouped[subject]) for subject in subjects])


def bootstrap_skill_interval(
    rows: Sequence[dict[str, str]],
    indices: np.ndarray,
    truth: np.ndarray,
    estimate: np.ndarray,
    baseline: np.ndarray,
    seed: int,
) -> tuple[float, float]:
    subjects, model_mse = subject_mse(rows, indices, truth, estimate)
    baseline_subjects, baseline_mse = subject_mse(rows, indices, truth, baseline)
    if subjects != baseline_subjects:
        raise AssertionError("Bootstrap subject ordering differs.")
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(subjects), size=(BOOTSTRAP_REPLICATES, len(subjects)))
    denominator = baseline_mse[draws].mean(axis=1)
    skills = 1.0 - model_mse[draws].mean(axis=1) / denominator
    return tuple(float(value) for value in np.quantile(skills, [0.025, 0.975]))


def evaluate_feature_space(
    *,
    features: np.ndarray,
    targets: np.ndarray,
    target_names: Sequence[str],
    target_types: Sequence[str],
    rows: Sequence[dict[str, str]],
    train_indices: np.ndarray,
    validation_indices: np.ndarray,
    development_indices: np.ndarray,
    test_indices: np.ndarray,
    alphas: Sequence[float],
    layer: str,
    representation: str,
    seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, dict[str, np.ndarray]]]:
    train_weights = subject_weights(rows, train_indices)
    validation_weights = subject_weights(rows, validation_indices)
    development_weights = subject_weights(rows, development_indices)
    test_weights = subject_weights(rows, test_indices)
    train_kernel, validation_cross = prepare_kernel(
        features, train_indices, validation_indices, train_weights
    )

    validation_records: list[dict[str, Any]] = []
    validation_predictions: dict[float, np.ndarray] = {}
    for alpha in alphas:
        prediction, _ = fit_from_kernel(
            train_kernel, validation_cross, targets[train_indices], train_weights, alpha
        )
        validation_predictions[float(alpha)] = prediction
        for column, name in enumerate(target_names):
            valid = np.isfinite(targets[validation_indices, column])
            if not np.any(valid):
                continue
            score = weighted_mse(
                targets[validation_indices, column][valid], prediction[:, column][valid],
                validation_weights[valid],
            )
            validation_records.append(
                {
                    "layer": layer,
                    "representation": representation,
                    "target": name,
                    "alpha": float(alpha),
                    "validation_mse": score,
                }
            )

    best_alpha: dict[str, float] = {}
    for name in target_names:
        candidates = [row for row in validation_records if row["target"] == name]
        if not candidates:
            raise ConceptAuditError(f"No finite validation data for {name}.")
        best_alpha[name] = float(min(candidates, key=lambda row: (row["validation_mse"], row["alpha"]))["alpha"])

    development_kernel, test_cross = prepare_kernel(
        features, development_indices, test_indices, development_weights
    )
    results: list[dict[str, Any]] = []
    stored: dict[str, dict[str, np.ndarray]] = {}
    for column, (name, target_type) in enumerate(zip(target_names, target_types)):
        alpha = best_alpha[name]
        prediction, target_mean = fit_from_kernel(
            development_kernel,
            test_cross,
            targets[development_indices, column : column + 1],
            development_weights,
            alpha,
        )
        truth = targets[test_indices, column]
        estimate = prediction[:, 0]
        baseline = np.full_like(truth, target_mean[0], dtype=np.float64)
        valid = np.isfinite(truth) & np.isfinite(estimate)
        if not np.any(valid):
            raise ConceptAuditError(f"No finite test data for {name}.")
        truth = truth[valid]
        estimate = estimate[valid]
        baseline = baseline[valid]
        valid_indices = test_indices[valid]
        valid_weights = test_weights[valid]
        model_mse = weighted_mse(truth, estimate, valid_weights)
        baseline_mse = weighted_mse(truth, baseline, valid_weights)
        development_truth = targets[development_indices, column]
        development_sd = math.sqrt(
            weighted_mse(
                development_truth,
                np.full_like(development_truth, target_mean[0]),
                development_weights,
            )
        )
        skill = 1.0 - model_mse / baseline_mse
        ci_low, ci_high = bootstrap_skill_interval(
            rows, valid_indices, truth, estimate, baseline,
            seed + int(hashlib.sha256(f"{layer}:{representation}:{name}".encode()).hexdigest()[:8], 16),
        )
        row: dict[str, Any] = {
            "layer": layer,
            "representation": representation,
            "target": name,
            "target_type": target_type,
            "selected_alpha": alpha,
            "development_n_subjects": len({rows[int(i)]["subject_id"] for i in development_indices}),
            "test_n_subjects": len({rows[int(i)]["subject_id"] for i in valid_indices}),
            "test_n_visits": len(valid_indices),
            "rmse": math.sqrt(model_mse),
            "normalized_rmse_development_sd": math.sqrt(model_mse) / development_sd,
            "correlation": weighted_correlation(truth, estimate, valid_weights),
            "skill_vs_development_mean": skill,
            "skill_ci_2.5_percent": ci_low,
            "skill_ci_97.5_percent": ci_high,
        }
        if target_type == "binary":
            row["auc"] = auc(truth, estimate)
            row["brier_score"] = model_mse
        else:
            row["auc"] = math.nan
            row["brier_score"] = math.nan
        results.append(row)
        stored[name] = {
            "indices": valid_indices,
            "truth": truth,
            "estimate": estimate,
            "baseline": baseline,
        }
    return results, validation_records, stored


def evaluate_split_skill(
    *,
    features: np.ndarray,
    target: np.ndarray,
    rows: Sequence[dict[str, str]],
    train_indices: np.ndarray,
    validation_indices: np.ndarray,
    development_indices: np.ndarray,
    test_indices: np.ndarray,
    alphas: Sequence[float],
) -> dict[str, float]:
    """A compact no-bootstrap probe used only for repeated-split stability."""

    train_weights = subject_weights(rows, train_indices)
    validation_weights = subject_weights(rows, validation_indices)
    development_weights = subject_weights(rows, development_indices)
    test_weights = subject_weights(rows, test_indices)
    train_kernel, validation_cross = prepare_kernel(
        features, train_indices, validation_indices, train_weights
    )
    tuning = []
    for alpha in alphas:
        prediction, _ = fit_from_kernel(
            train_kernel, validation_cross, target[train_indices, None], train_weights, alpha
        )
        tuning.append(
            (
                weighted_mse(target[validation_indices], prediction[:, 0], validation_weights),
                float(alpha),
            )
        )
    _, selected_alpha = min(tuning)
    development_kernel, test_cross = prepare_kernel(
        features, development_indices, test_indices, development_weights
    )
    prediction, target_mean = fit_from_kernel(
        development_kernel,
        test_cross,
        target[development_indices, None],
        development_weights,
        selected_alpha,
    )
    truth = target[test_indices]
    estimate = prediction[:, 0]
    baseline = np.full_like(truth, target_mean[0])
    model_mse = weighted_mse(truth, estimate, test_weights)
    baseline_mse = weighted_mse(truth, baseline, test_weights)
    return {
        "selected_alpha": selected_alpha,
        "skill_vs_development_mean": 1.0 - model_mse / baseline_mse,
        "correlation": weighted_correlation(truth, estimate, test_weights),
    }


def summarize_stability(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault((row["layer"], row["target"], row["representation"]), []).append(row)
    output = []
    for (layer, target, representation), values in sorted(grouped.items()):
        skill = np.asarray([row["skill_vs_development_mean"] for row in values], dtype=float)
        correlation = np.asarray([row["correlation"] for row in values], dtype=float)
        output.append(
            {
                "layer": layer,
                "target": target,
                "representation": representation,
                "splits": len(values),
                "median_skill": float(np.median(skill)),
                "skill_q25": float(np.quantile(skill, 0.25)),
                "skill_q75": float(np.quantile(skill, 0.75)),
                "fraction_skill_positive": float(np.mean(skill > 0)),
                "median_correlation": float(np.nanmedian(correlation)),
            }
        )
    return output


def summarize_specificity_stability(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    lookup = {
        (row["layer"], row["target"], row["repeat"], row["representation"]): row
        for row in rows
    }
    output = []
    for layer in sorted({row["layer"] for row in rows}):
        for target in LOCAL_TARGETS:
            deltas = []
            repeats = sorted(
                {row["repeat"] for row in rows if row["layer"] == layer and row["target"] == target}
            )
            for repeat in repeats:
                roi = lookup[(layer, target, repeat, f"{target}_roi_tokens")]
                remote = lookup[(layer, target, repeat, f"{target}_remote_tokens")]
                deltas.append(roi["skill_vs_development_mean"] - remote["skill_vs_development_mean"])
            values = np.asarray(deltas, dtype=float)
            output.append(
                {
                    "layer": layer,
                    "target": target,
                    "splits": len(values),
                    "median_delta_skill_roi_minus_remote": float(np.median(values)),
                    "delta_q25": float(np.quantile(values, 0.25)),
                    "delta_q75": float(np.quantile(values, 0.75)),
                    "fraction_roi_better_than_remote": float(np.mean(values > 0)),
                }
            )
    return output


def log_jacobian_path(nifti_path: str, nonlinear_dir: Path) -> Path:
    name = Path(nifti_path).name
    suffix = "_Warped.nii.gz"
    if not name.endswith(suffix):
        raise ConceptAuditError(f"Cannot derive log-Jacobian name from {name}")
    return nonlinear_dir / f"{name[:-len(suffix)]}_LogJacobian_masked.nii.gz"


def valid_mask_path(log_jacobian: Path) -> Path:
    suffix = "_LogJacobian_masked.nii.gz"
    return log_jacobian.with_name(
        f"{log_jacobian.name[:-len(suffix)]}_analysisValidMask_templateSpace.nii.gz"
    )


def load_physical_targets(
    rows: Sequence[dict[str, str]], atlas_path: Path, nonlinear_dir: Path
) -> tuple[np.ndarray, Any, list[dict[str, Any]]]:
    runtime = extractor.import_runtime()
    nibabel = runtime["nibabel"]
    atlas_image = nibabel.load(str(atlas_path))
    atlas = np.asanyarray(atlas_image.dataobj).astype(np.int16)
    masks = {name: np.isin(atlas, labels) for name, labels in ROI_LABELS.items()}
    targets = np.empty((len(rows), len(ROI_LABELS)), dtype=np.float64)
    coverage = np.empty_like(targets)
    for position, row in enumerate(rows):
        log_path = log_jacobian_path(row["nifti_path"], nonlinear_dir)
        mask_path = valid_mask_path(log_path)
        if not log_path.is_file() or not mask_path.is_file():
            missing = log_path if not log_path.is_file() else mask_path
            raise ConceptAuditError(f"Missing deformation-map input: {missing}")
        image = nibabel.load(str(log_path))
        valid_image = nibabel.load(str(mask_path))
        if image.shape != atlas_image.shape or valid_image.shape != atlas_image.shape:
            raise ConceptAuditError(f"Atlas/deformation shape mismatch for {log_path}")
        if not np.allclose(image.affine, atlas_image.affine, atol=1e-5) or not np.allclose(
            valid_image.affine, atlas_image.affine, atol=1e-5
        ):
            raise ConceptAuditError(f"Atlas/deformation affine mismatch for {log_path}")
        values = np.asanyarray(image.dataobj).astype(np.float32)
        valid_map = np.asanyarray(valid_image.dataobj) > 0.5
        valid_map &= np.isfinite(values)
        for column, (name, roi_mask) in enumerate(masks.items()):
            valid_roi = roi_mask & valid_map
            coverage[position, column] = valid_roi.sum() / roi_mask.sum()
            if valid_roi.sum() < 100:
                raise ConceptAuditError(f"Too few valid voxels for {name} in {log_path}")
            targets[position, column] = float(values[valid_roi].mean(dtype=np.float64))
        if (position + 1) % 20 == 0 or position + 1 == len(rows):
            log(f"loaded deformation targets {position + 1}/{len(rows)}")
    qc = []
    for column, name in enumerate(ROI_LABELS):
        qc.append(
            {
                "target": name,
                "atlas_voxels": int(masks[name].sum()),
                "mean_valid_fraction": float(coverage[:, column].mean()),
                "minimum_valid_fraction": float(coverage[:, column].min()),
                "target_mean": float(targets[:, column].mean()),
                "target_sd": float(targets[:, column].std(ddof=1)),
                "target_minimum": float(targets[:, column].min()),
                "target_maximum": float(targets[:, column].max()),
            }
        )
    return targets, atlas_image, qc


def resize_atlas(atlas_image: Any) -> np.ndarray:
    runtime = extractor.import_runtime()
    torch = runtime["torch"]
    atlas = np.asanyarray(atlas_image.dataobj).astype(np.float32)
    tensor = torch.as_tensor(atlas)[None, None]
    resized = torch.nn.functional.interpolate(
        tensor, size=extractor.IMAGE_SIZE, mode="nearest"
    )[0, 0]
    return resized.numpy().astype(np.int16)


def patch_counts(mask: np.ndarray) -> np.ndarray:
    grid = extractor.PATCH_GRID
    patch = extractor.PATCH_SIZE
    return (
        mask.reshape(grid[0], patch[0], grid[1], patch[1], grid[2], patch[2])
        .transpose(0, 2, 4, 1, 3, 5)
        .reshape(extractor.EXPECTED_TOKENS, -1)
        .sum(axis=1)
    )


def select_tokens(atlas_96: np.ndarray) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], list[dict[str, Any]]]:
    brain_counts = patch_counts(atlas_96 > 0)
    brain_tokens = np.flatnonzero(brain_counts > 0)
    grid_coordinates = np.stack(
        np.unravel_index(np.arange(extractor.EXPECTED_TOKENS), extractor.PATCH_GRID), axis=1
    )
    roi_tokens: dict[str, np.ndarray] = {}
    remote_tokens: dict[str, np.ndarray] = {}
    definitions: list[dict[str, Any]] = []
    for name in LOCAL_TARGETS:
        counts = patch_counts(np.isin(atlas_96, ROI_LABELS[name]))
        candidates = np.flatnonzero(counts > 0)
        ordered = candidates[np.argsort(-counts[candidates], kind="stable")]
        selected = np.sort(ordered[:MAX_LOCAL_TOKENS])
        centroid = np.average(grid_coordinates[candidates], axis=0, weights=counts[candidates])
        remote_candidates = np.setdiff1d(brain_tokens, candidates, assume_unique=True)
        distances = np.sum((grid_coordinates[remote_candidates] - centroid) ** 2, axis=1)
        remote_order = np.lexsort((remote_candidates, -distances))
        remote = np.sort(remote_candidates[remote_order[: len(selected)]])
        if len(remote) != len(selected):
            raise ConceptAuditError(f"Could not match remote tokens for {name}")
        roi_tokens[name] = selected
        remote_tokens[name] = remote
        for group, tokens in (("roi", selected), ("remote", remote)):
            for token in tokens:
                coordinate = grid_coordinates[token]
                definitions.append(
                    {
                        "target": name,
                        "group": group,
                        "token_index_zero_based": int(token),
                        "grid_i": int(coordinate[0]),
                        "grid_j": int(coordinate[1]),
                        "grid_k": int(coordinate[2]),
                        "roi_voxels_in_token": int(counts[token]),
                        "brain_voxels_in_token": int(brain_counts[token]),
                    }
                )
    return roi_tokens, remote_tokens, definitions


def baseline_indices(rows: Sequence[dict[str, str]]) -> np.ndarray:
    first: dict[str, tuple[float, int, int]] = {}
    for index, row in enumerate(rows):
        key = (float(row["time_month"]), int(row["image_id"]), index)
        subject = row["subject_id"]
        if subject not in first or key < first[subject]:
            first[subject] = key
    return np.asarray([first[subject][2] for subject in sorted(first)], dtype=np.int64)


def clinical_matrix(rows: Sequence[dict[str, str]]) -> tuple[np.ndarray, list[str], list[str]]:
    values = np.empty((len(rows), len(CLINICAL_TARGETS)), dtype=np.float64)
    names = list(CLINICAL_TARGETS)
    types = [CLINICAL_TARGETS[name][0] for name in names]
    for column, name in enumerate(names):
        target_type, source = CLINICAL_TARGETS[name]
        values[:, column] = [
            parse_binary(row.get(source, "")) if target_type == "binary" else finite_float(row.get(source, ""))
            for row in rows
        ]
    return values, names, types


def subset_indices(indices: np.ndarray, allowed: set[int]) -> np.ndarray:
    return np.asarray([value for value in indices if int(value) in allowed], dtype=np.int64)


def local_specificity_rows(
    stored: dict[tuple[str, str, str], dict[str, np.ndarray]],
    rows: Sequence[dict[str, str]],
    seed: int,
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for layer in sorted({key[0] for key in stored}):
        for target in LOCAL_TARGETS:
            roi = stored[(layer, f"{target}_roi_tokens", target)]
            remote = stored[(layer, f"{target}_remote_tokens", target)]
            if not np.array_equal(roi["indices"], remote["indices"]):
                raise AssertionError("ROI and remote probes have different test visits.")
            indices = roi["indices"]
            subjects, roi_mse = subject_mse(rows, indices, roi["truth"], roi["estimate"])
            remote_subjects, remote_mse = subject_mse(
                rows, indices, remote["truth"], remote["estimate"]
            )
            baseline_subjects, baseline_mse = subject_mse(
                rows, indices, roi["truth"], roi["baseline"]
            )
            if subjects != remote_subjects or subjects != baseline_subjects:
                raise AssertionError("Specificity comparison subject mismatch.")
            delta = float((remote_mse.mean() - roi_mse.mean()) / baseline_mse.mean())
            rng = np.random.default_rng(
                seed + int(hashlib.sha256(f"specificity:{layer}:{target}".encode()).hexdigest()[:8], 16)
            )
            draws = rng.integers(0, len(subjects), size=(BOOTSTRAP_REPLICATES, len(subjects)))
            boot = (
                remote_mse[draws].mean(axis=1) - roi_mse[draws].mean(axis=1)
            ) / baseline_mse[draws].mean(axis=1)
            low, high = np.quantile(boot, [0.025, 0.975])
            output.append(
                {
                    "layer": layer,
                    "target": target,
                    "delta_skill_roi_minus_remote": delta,
                    "delta_ci_2.5_percent": float(low),
                    "delta_ci_97.5_percent": float(high),
                    "spatially_specific_positive": bool(delta > 0),
                }
            )
    return output


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    reconstruction.write_csv(path, rows)


def draw_token_map(
    path: Path, atlas_96: np.ndarray, roi_tokens: dict[str, np.ndarray], remote_tokens: dict[str, np.ndarray]
) -> None:
    from PIL import Image, ImageDraw, ImageFont

    colors = {
        "hippocampus": (214, 39, 40),
        "lateral_ventricles": (31, 119, 180),
        "amygdala": (255, 127, 14),
        "thalamus": (44, 160, 44),
    }
    font = ImageFont.load_default()
    tile = 230
    canvas = Image.new("RGB", (tile * 2, tile * 2 + 45), "white")
    draw = ImageDraw.Draw(canvas)
    draw.text((8, 7), "Atlas ROIs and selected 16^3 BrainIAC tokens (thick=ROI, thin=remote)", fill="black", font=font)
    for panel, name in enumerate(LOCAL_TARGETS):
        mask = np.isin(atlas_96, ROI_LABELS[name])
        coordinates = np.argwhere(mask)
        slice_index = int(np.rint(coordinates[:, 2].mean()))
        background = np.asarray(np.rot90(atlas_96[:, :, slice_index] > 0), dtype=np.uint8) * 225
        rgb = np.repeat(background[:, :, None], 3, axis=2)
        roi_slice = np.rot90(mask[:, :, slice_index])
        rgb[roi_slice] = colors[name]
        image = Image.fromarray(rgb).resize((192, 192), Image.Resampling.NEAREST)
        x0 = (panel % 2) * tile + 19
        y0 = (panel // 2) * tile + 39
        canvas.paste(image, (x0, y0))
        draw.text((x0, y0 - 13), name.replace("_", " "), fill=colors[name], font=font)
        # Axial slice: token grid i/j is drawn after np.rot90, with k selecting the slice.
        for group, tokens in (("roi", roi_tokens[name]), ("remote", remote_tokens[name])):
            for token in tokens:
                i, j, k = np.unravel_index(int(token), extractor.PATCH_GRID)
                if k != min(5, max(0, slice_index // 16)):
                    continue
                left = x0 + i * 32
                top = y0 + (5 - j) * 32
                width = 3 if group == "roi" else 1
                draw.rectangle((left, top, left + 31, top + 31), outline=colors[name], width=width)
    canvas.save(path, optimize=True)


def draw_performance_plot(path: Path, results: Sequence[dict[str, Any]], physical: bool) -> None:
    from PIL import Image, ImageDraw, ImageFont

    if physical:
        targets = list(ROI_LABELS)
        representations = ("mean_pool", "all_tokens")
        title = "Held-out physical concept skill (1 - MSE / development-mean MSE)"
    else:
        targets = list(CLINICAL_TARGETS)
        representations = ("mean_pool", "all_tokens", "hippocampus_roi_tokens")
        title = "Exploratory held-out clinical probe skill (30-subject pilot)"
    layers = ("block_3", "block_6", "final")
    colors = {"block_3": (31, 119, 180), "block_6": (255, 127, 14), "final": (44, 160, 44)}
    width, row_height = 960, 30
    height = 75 + len(targets) * len(representations) * row_height
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()
    draw.text((10, 8), title, fill="black", font=font)
    draw.text((10, 23), "Dots: layer; horizontal lines: subject-bootstrap 95% interval; zero = mean-only baseline", fill="black", font=font)
    x_left, x_right = 310, 940
    low_limit, high_limit = -2.0, 1.0
    def xpos(value: float) -> int:
        value = max(low_limit, min(high_limit, value))
        return int(x_left + (value - low_limit) / (high_limit - low_limit) * (x_right - x_left))
    for tick in (-2, -1, 0, 0.5, 1):
        x = xpos(tick)
        draw.line((x, 48, x, height - 10), fill=(190, 190, 190) if tick else (70, 70, 70), width=1)
        draw.text((x - 9, 37), str(tick), fill="black", font=font)
    lookup = {(row["layer"], row["representation"], row["target"]): row for row in results}
    y = 58
    for target in targets:
        for representation in representations:
            label = f"{target.replace('_', ' ')} | {representation.replace('_', ' ')}"
            draw.text((8, y + 7), label, fill="black", font=font)
            for offset, layer in enumerate(layers):
                row = lookup.get((layer, representation, target))
                if row is None:
                    continue
                yy = y + 8 + (offset - 1) * 5
                left = xpos(float(row["skill_ci_2.5_percent"]))
                right = xpos(float(row["skill_ci_97.5_percent"]))
                center = xpos(float(row["skill_vs_development_mean"]))
                draw.line((left, yy, right, yy), fill=colors[layer], width=2)
                draw.ellipse((center - 3, yy - 3, center + 3, yy + 3), fill=colors[layer])
            y += row_height
    canvas.save(path, optimize=True)


def draw_specificity_plot(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    from PIL import Image, ImageDraw, ImageFont

    layers = ("block_3", "block_6", "final")
    colors = {"block_3": (31, 119, 180), "block_6": (255, 127, 14), "final": (44, 160, 44)}
    width, height = 880, 90 + 58 * len(LOCAL_TARGETS)
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()
    draw.text((10, 8), "Local spatial specificity: ROI-token skill minus distant-token skill", fill="black", font=font)
    draw.text((10, 23), "Positive values support location-specific information; intervals resample held-out subjects", fill="black", font=font)
    x_left, x_right = 220, 850
    low_limit, high_limit = -2.0, 2.0
    def xpos(value: float) -> int:
        value = max(low_limit, min(high_limit, value))
        return int(x_left + (value - low_limit) / (high_limit - low_limit) * (x_right - x_left))
    for tick in (-2, -1, 0, 1, 2):
        x = xpos(tick)
        draw.line((x, 50, x, height - 12), fill=(70, 70, 70) if tick == 0 else (200, 200, 200))
        draw.text((x - 8, 38), str(tick), fill="black", font=font)
    lookup = {(row["layer"], row["target"]): row for row in rows}
    for target_position, target in enumerate(LOCAL_TARGETS):
        y = 68 + target_position * 58
        draw.text((8, y + 15), target.replace("_", " "), fill="black", font=font)
        for layer_position, layer in enumerate(layers):
            row = lookup[(layer, target)]
            yy = y + layer_position * 13
            left = xpos(float(row["delta_ci_2.5_percent"]))
            right = xpos(float(row["delta_ci_97.5_percent"]))
            center = xpos(float(row["delta_skill_roi_minus_remote"]))
            draw.line((left, yy, right, yy), fill=colors[layer], width=2)
            draw.ellipse((center - 4, yy - 4, center + 4, yy + 4), fill=colors[layer])
            draw.text((x_right - 65, yy - 5), layer.replace("block_", "B"), fill=colors[layer], font=font)
    canvas.save(path, optimize=True)


def draw_local_absolute_plot(path: Path, results: Sequence[dict[str, Any]]) -> None:
    from PIL import Image, ImageDraw, ImageFont

    layers = ("block_3", "block_6", "final")
    colors = {"block_3": (31, 119, 180), "block_6": (255, 127, 14), "final": (44, 160, 44)}
    width, row_height = 920, 32
    height = 75 + len(LOCAL_TARGETS) * 2 * row_height
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()
    draw.text((10, 8), "Absolute local-target skill: anatomical ROI tokens versus distant tokens", fill="black", font=font)
    draw.text((10, 23), "A local claim needs positive ROI skill as well as an ROI-over-remote advantage", fill="black", font=font)
    x_left, x_right = 280, 890
    low_limit, high_limit = -3.0, 1.0
    def xpos(value: float) -> int:
        value = max(low_limit, min(high_limit, value))
        return int(x_left + (value - low_limit) / (high_limit - low_limit) * (x_right - x_left))
    for tick in (-3, -2, -1, 0, 0.5, 1):
        x = xpos(tick)
        draw.line((x, 48, x, height - 10), fill=(70, 70, 70) if tick == 0 else (200, 200, 200))
        draw.text((x - 9, 37), str(tick), fill="black", font=font)
    lookup = {(row["layer"], row["representation"], row["target"]): row for row in results}
    y = 59
    for target in LOCAL_TARGETS:
        for group in ("roi", "remote"):
            representation = f"{target}_{group}_tokens"
            draw.text((8, y + 8), f"{target.replace('_', ' ')} | {group}", fill="black", font=font)
            for offset, layer in enumerate(layers):
                row = lookup[(layer, representation, target)]
                yy = y + 9 + (offset - 1) * 6
                left = xpos(float(row["skill_ci_2.5_percent"]))
                right = xpos(float(row["skill_ci_97.5_percent"]))
                center = xpos(float(row["skill_vs_development_mean"]))
                draw.line((left, yy, right, yy), fill=colors[layer], width=2)
                draw.ellipse((center - 3, yy - 3, center + 3, yy + 3), fill=colors[layer])
            y += row_height
    canvas.save(path, optimize=True)


def parse_layer(value: str) -> tuple[str, Path]:
    return reconstruction.parse_layer(value)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--layer", action="append", type=parse_layer)
    parser.add_argument("--atlas", type=Path, default=Path(DEFAULT_ATLAS))
    parser.add_argument("--nonlinear-dir", type=Path, default=Path(DEFAULT_NONLINEAR_DIR))
    parser.add_argument("--output-dir", type=Path, default=Path(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--test-fraction", type=float, default=0.20)
    parser.add_argument("--validation-fraction", type=float, default=0.20)
    parser.add_argument("--alpha-ratios", type=float, nargs="+", default=list(DEFAULT_ALPHA_RATIOS))
    parser.add_argument("--stability-repeats", type=int, default=DEFAULT_STABILITY_REPEATS)
    args = parser.parse_args()
    if not 0 < args.test_fraction < 0.5 or not 0 < args.validation_fraction < 0.5:
        parser.error("split fractions must lie strictly between 0 and 0.5")
    if args.stability_repeats < 0:
        parser.error("--stability-repeats cannot be negative")
    return args


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parents[2]
    layer_args = args.layer or [
        (name, (project_root / relative).resolve())
        for name, relative in reconstruction.DEFAULT_LAYER_PATHS
    ]
    layers = [(name, path, reconstruction.read_metadata(path)) for name, path in layer_args]
    reconstruction.validate_layer_contracts(layers)
    rows = layers[0][2]
    train_subjects, validation_subjects, test_subjects = reconstruction.make_subject_split(
        rows, args.seed, args.test_fraction, args.validation_fraction
    )
    development_subjects = train_subjects | validation_subjects
    train_indices = reconstruction.indices_for(rows, train_subjects)
    validation_indices = reconstruction.indices_for(rows, validation_subjects)
    development_indices = reconstruction.indices_for(rows, development_subjects)
    test_indices = reconstruction.indices_for(rows, test_subjects)
    log(
        f"subject split: {len(train_subjects)} tune-train, {len(validation_subjects)} "
        f"tune-validation, {len(test_subjects)} final-test"
    )

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    atlas_path = args.atlas.expanduser().resolve()
    nonlinear_dir = args.nonlinear_dir.expanduser().resolve()
    if not atlas_path.is_file() or not nonlinear_dir.is_dir():
        raise ConceptAuditError("Atlas or nonlinear deformation-map directory is missing.")

    log("loading atlas-defined regional deformation targets")
    physical_targets, atlas_image, target_qc = load_physical_targets(rows, atlas_path, nonlinear_dir)
    atlas_96 = resize_atlas(atlas_image)
    roi_tokens, remote_tokens, token_definitions = select_tokens(atlas_96)
    write_csv(output_dir / "physical_target_qc.csv", target_qc)
    write_csv(output_dir / "token_definitions.csv", token_definitions)
    draw_token_map(output_dir / "anatomical_targets_and_tokens.png", atlas_96, roi_tokens, remote_tokens)

    baseline = baseline_indices(rows)
    baseline_set = set(int(value) for value in baseline)
    clinical_targets, clinical_names, clinical_types = clinical_matrix(rows)
    clinical_train = subset_indices(train_indices, baseline_set)
    clinical_validation = subset_indices(validation_indices, baseline_set)
    clinical_development = subset_indices(development_indices, baseline_set)
    clinical_test = subset_indices(test_indices, baseline_set)
    for column, name in enumerate(clinical_names):
        for split_name, indices in (("train", clinical_train), ("validation", clinical_validation), ("test", clinical_test)):
            if not np.isfinite(clinical_targets[indices, column]).all():
                raise ConceptAuditError(f"Missing {name} in the {split_name} baseline subset.")

    physical_results: list[dict[str, Any]] = []
    clinical_results: list[dict[str, Any]] = []
    tuning_rows: list[dict[str, Any]] = []
    stability_rows: list[dict[str, Any]] = []
    stored_local: dict[tuple[str, str, str], dict[str, np.ndarray]] = {}
    for layer_name, _, layer_rows in layers:
        log(f"loading {layer_name} token matrices")
        tokens = reconstruction.load_features(layer_rows).reshape(
            len(rows), extractor.EXPECTED_TOKENS, extractor.EXPECTED_FEATURES
        )
        common_spaces = {
            "mean_pool": tokens.mean(axis=1, dtype=np.float32),
            "all_tokens": tokens.reshape(len(rows), -1),
        }
        for representation, features in common_spaces.items():
            log(f"probing {layer_name} / {representation} physical targets")
            results, tuning, _ = evaluate_feature_space(
                features=features,
                targets=physical_targets,
                target_names=list(ROI_LABELS),
                target_types=["continuous"] * len(ROI_LABELS),
                rows=rows,
                train_indices=train_indices,
                validation_indices=validation_indices,
                development_indices=development_indices,
                test_indices=test_indices,
                alphas=args.alpha_ratios,
                layer=layer_name,
                representation=representation,
                seed=args.seed,
            )
            physical_results.extend(results)
            tuning_rows.extend({"audit": "physical", **row} for row in tuning)

            log(f"probing {layer_name} / {representation} clinical targets")
            results, tuning, _ = evaluate_feature_space(
                features=features,
                targets=clinical_targets,
                target_names=clinical_names,
                target_types=clinical_types,
                rows=rows,
                train_indices=clinical_train,
                validation_indices=clinical_validation,
                development_indices=clinical_development,
                test_indices=clinical_test,
                alphas=args.alpha_ratios,
                layer=layer_name,
                representation=representation,
                seed=args.seed,
            )
            clinical_results.extend(results)
            tuning_rows.extend({"audit": "clinical", **row} for row in tuning)

        for target in LOCAL_TARGETS:
            column = list(ROI_LABELS).index(target)
            for group, selected in (("roi", roi_tokens[target]), ("remote", remote_tokens[target])):
                representation = f"{target}_{group}_tokens"
                features = tokens[:, selected, :].reshape(len(rows), -1)
                log(f"probing {layer_name} / {representation}")
                results, tuning, stored = evaluate_feature_space(
                    features=features,
                    targets=physical_targets[:, column : column + 1],
                    target_names=[target],
                    target_types=["continuous"],
                    rows=rows,
                    train_indices=train_indices,
                    validation_indices=validation_indices,
                    development_indices=development_indices,
                    test_indices=test_indices,
                    alphas=args.alpha_ratios,
                    layer=layer_name,
                    representation=representation,
                    seed=args.seed,
                )
                physical_results.extend(results)
                tuning_rows.extend({"audit": "physical", **row} for row in tuning)
                stored_local[(layer_name, representation, target)] = stored[target]

            # A hippocampal token probe is a useful local clinical comparator.
            if target == "hippocampus":
                representation = "hippocampus_roi_tokens"
                features = tokens[:, roi_tokens[target], :].reshape(len(rows), -1)
                results, tuning, _ = evaluate_feature_space(
                    features=features,
                    targets=clinical_targets,
                    target_names=clinical_names,
                    target_types=clinical_types,
                    rows=rows,
                    train_indices=clinical_train,
                    validation_indices=clinical_validation,
                    development_indices=clinical_development,
                    test_indices=clinical_test,
                    alphas=args.alpha_ratios,
                    layer=layer_name,
                    representation=representation,
                    seed=args.seed,
                )
                clinical_results.extend(results)
                tuning_rows.extend({"audit": "clinical", **row} for row in tuning)

        if args.stability_repeats:
            log(f"running {args.stability_repeats} repeated subject splits for {layer_name}")
            stability_spaces: dict[str, tuple[np.ndarray, Sequence[str]]] = {
                "mean_pool": (common_spaces["mean_pool"], tuple(ROI_LABELS)),
                "all_tokens": (common_spaces["all_tokens"], tuple(ROI_LABELS)),
            }
            for target in LOCAL_TARGETS:
                stability_spaces[f"{target}_roi_tokens"] = (
                    tokens[:, roi_tokens[target], :].reshape(len(rows), -1), (target,)
                )
                stability_spaces[f"{target}_remote_tokens"] = (
                    tokens[:, remote_tokens[target], :].reshape(len(rows), -1), (target,)
                )
            for repeat in range(args.stability_repeats):
                repeat_seed = args.seed + 1000 + repeat
                split_train, split_validation, split_test = reconstruction.make_subject_split(
                    rows, repeat_seed, args.test_fraction, args.validation_fraction
                )
                split_development = split_train | split_validation
                split_indices = (
                    reconstruction.indices_for(rows, split_train),
                    reconstruction.indices_for(rows, split_validation),
                    reconstruction.indices_for(rows, split_development),
                    reconstruction.indices_for(rows, split_test),
                )
                for representation, (features, names) in stability_spaces.items():
                    for target in names:
                        column = list(ROI_LABELS).index(target)
                        result = evaluate_split_skill(
                            features=features,
                            target=physical_targets[:, column],
                            rows=rows,
                            train_indices=split_indices[0],
                            validation_indices=split_indices[1],
                            development_indices=split_indices[2],
                            test_indices=split_indices[3],
                            alphas=args.alpha_ratios,
                        )
                        stability_rows.append(
                            {
                                "layer": layer_name,
                                "repeat": repeat + 1,
                                "split_seed": repeat_seed,
                                "target": target,
                                "representation": representation,
                                **result,
                            }
                        )
        del tokens

    specificity = local_specificity_rows(stored_local, rows, args.seed)
    write_csv(output_dir / "physical_probe_metrics.csv", physical_results)
    write_csv(output_dir / "clinical_probe_metrics.csv", clinical_results)
    write_csv(output_dir / "local_spatial_specificity.csv", specificity)
    write_csv(output_dir / "ridge_tuning.csv", tuning_rows)
    if stability_rows:
        write_csv(output_dir / "split_stability_metrics.csv", stability_rows)
        write_csv(output_dir / "split_stability_summary.csv", summarize_stability(stability_rows))
        write_csv(
            output_dir / "local_specificity_stability.csv",
            summarize_specificity_stability(stability_rows),
        )
    split_rows = [
        {
            "anonymous_subject": f"subject_{position + 1:02d}",
            "split": "tune_train" if subject in train_subjects else (
                "tune_validation" if subject in validation_subjects else "final_test"
            ),
            "visits": sum(row["subject_id"] == subject for row in rows),
        }
        for position, subject in enumerate(sorted(train_subjects | validation_subjects | test_subjects))
    ]
    write_csv(output_dir / "subject_split_anonymous.csv", split_rows)
    draw_performance_plot(output_dir / "physical_probe_performance.png", physical_results, True)
    draw_performance_plot(output_dir / "clinical_probe_performance.png", clinical_results, False)
    draw_specificity_plot(output_dir / "local_spatial_specificity.png", specificity)
    draw_local_absolute_plot(output_dir / "local_probe_absolute_performance.png", physical_results)

    manifest = {
        "script_version": SCRIPT_VERSION,
        "seed": args.seed,
        "matrixfpca_used": False,
        "cohort_subjects": len(train_subjects | validation_subjects | test_subjects),
        "cohort_visits": len(rows),
        "tune_train_subjects": len(train_subjects),
        "tune_validation_subjects": len(validation_subjects),
        "final_test_subjects": len(test_subjects),
        "atlas": str(atlas_path),
        "atlas_sha256": reconstruction.sha256_file(atlas_path),
        "nonlinear_dir": str(nonlinear_dir),
        "layers": [{"name": name, "metadata": str(path)} for name, path, _ in layers],
        "roi_labels": {name: list(labels) for name, labels in ROI_LABELS.items()},
        "local_token_limit": MAX_LOCAL_TOKENS,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "stability_repeats": args.stability_repeats,
        "alpha_ratios": list(args.alpha_ratios),
        "physical_target": "regional mean visit-to-template log-Jacobian over valid atlas voxels",
        "clinical_scope": "exploratory; one earliest scan per subject",
    }
    with (output_dir / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    log(f"complete: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ConceptAuditError, reconstruction.AuditError, extractor.ExtractionError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        raise SystemExit(2)
