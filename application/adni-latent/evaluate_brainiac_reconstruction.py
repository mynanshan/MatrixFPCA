#!/usr/bin/env python3
"""Audit image information retained by frozen BrainIAC features.

BrainIAC itself has no image decoder.  This script fits a deliberately simple
linear kernel-ridge *probe* from a complete 216 x 768 token matrix to the
normalized 96 x 96 x 96 model input.  Subjects, rather than visits, are split
between decoder development and test sets.  MatrixFPCA is never loaded.

The probe answers a limited question: is image-space information linearly
decodable from a frozen BrainIAC layer on held-out subjects?  It is not an
official BrainIAC decoder, does not recover native NIfTI resolution or native
intensity units, and does not establish that arbitrary MatFPCA-produced latent
matrices lie on the decoder's training distribution.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import sys
import time
from typing import Any, Iterable, Sequence

import numpy as np

import extract_brainiac_latents as extractor


SCRIPT_VERSION = "1.0.0"
DEFAULT_SEED = 20260819
DEFAULT_LAYER_PATHS = (
    ("block_3", "application/adni-latent/latents/brainiac_patch_layer03/metadata.csv"),
    ("block_6", "application/adni-latent/latents/brainiac_patch_layer06/metadata.csv"),
    ("final", "application/adni-latent/latents/brainiac_patch/metadata.csv"),
)
DEFAULT_OUTPUT_DIR = "application/adni-latent/reconstruction_probe"
DEFAULT_ALPHA_RATIOS = (
    1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0,
)


class AuditError(RuntimeError):
    """A user-actionable reconstruction-audit failure."""


def log(message: str) -> None:
    print(time.strftime("%Y-%m-%d %H:%M:%S"), "|", message, flush=True)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_path(value: str, metadata_path: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = metadata_path.parent / path
    return path.resolve()


def read_metadata(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise AuditError(f"Metadata does not exist: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise AuditError(f"Metadata has no header: {path}")
        required = {
            "subject_id", "image_id", "time_month", "nifti_path", "latent_path",
            "latent_rows", "latent_cols", "checkpoint_sha256", "encoder_layer",
        }
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise AuditError(f"{path} is missing columns: {', '.join(missing)}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise AuditError(f"Metadata contains no visits: {path}")
    for row in rows:
        row["subject_id"] = row["subject_id"].strip()
        row["image_id"] = row["image_id"].strip()
        row["nifti_path"] = str(resolve_path(row["nifti_path"], path))
        row["latent_path"] = str(resolve_path(row["latent_path"], path))
        if (int(row["latent_rows"]), int(row["latent_cols"])) != (
            extractor.EXPECTED_TOKENS,
            extractor.EXPECTED_FEATURES,
        ):
            raise AuditError(
                f"Unexpected latent shape for image {row['image_id']} in {path}."
            )
    rows.sort(key=lambda row: (row["subject_id"], float(row["time_month"]), int(row["image_id"])))
    return rows


def validate_layer_contracts(
    layers: Sequence[tuple[str, Path, list[dict[str, str]]]],
) -> None:
    reference = layers[0][2]
    reference_keys = [(row["subject_id"], row["image_id"]) for row in reference]
    reference_niftis = [row["nifti_path"] for row in reference]
    checkpoint_hashes: set[str] = set()
    for name, path, rows in layers:
        keys = [(row["subject_id"], row["image_id"]) for row in rows]
        if keys != reference_keys:
            raise AuditError(f"Visit order/cohort for {name} does not match {layers[0][0]}.")
        if [row["nifti_path"] for row in rows] != reference_niftis:
            raise AuditError(f"NIfTI paths for {name} do not match {layers[0][0]}.")
        layer_values = {row["encoder_layer"] for row in rows}
        if len(layer_values) != 1:
            raise AuditError(f"{path} mixes encoder layers: {sorted(layer_values)}")
        checkpoint_hashes.update(row["checkpoint_sha256"] for row in rows)
        for row in rows:
            if not Path(row["latent_path"]).is_file():
                raise AuditError(f"Missing latent: {row['latent_path']}")
    if len(checkpoint_hashes) != 1:
        raise AuditError("Layer metadata do not use one common BrainIAC checkpoint.")
    missing_niftis = [value for value in reference_niftis if not Path(value).is_file()]
    if missing_niftis:
        raise AuditError(f"Missing source NIfTI: {missing_niftis[0]}")


def make_subject_split(
    rows: Sequence[dict[str, str]], seed: int, test_fraction: float, validation_fraction: float
) -> tuple[set[str], set[str], set[str]]:
    subjects = sorted({row["subject_id"] for row in rows})
    if len(subjects) < 5:
        raise AuditError("At least five subjects are required for a nested subject split.")
    rng = random.Random(seed)
    n_test = max(1, round(test_fraction * len(subjects)))
    test = set(rng.sample(subjects, n_test))
    development = [subject for subject in subjects if subject not in test]
    n_validation = max(1, round(validation_fraction * len(development)))
    validation = set(rng.sample(development, n_validation))
    train = set(development) - validation
    if not train:
        raise AuditError("Nested split left no decoder-training subjects.")
    return train, validation, test


def indices_for(rows: Sequence[dict[str, str]], subjects: set[str]) -> np.ndarray:
    return np.asarray(
        [index for index, row in enumerate(rows) if row["subject_id"] in subjects],
        dtype=np.int64,
    )


def load_features(rows: Sequence[dict[str, str]]) -> np.ndarray:
    feature_count = extractor.EXPECTED_TOKENS * extractor.EXPECTED_FEATURES
    features = np.empty((len(rows), feature_count), dtype=np.float32)
    for index, row in enumerate(rows):
        latent = np.load(row["latent_path"], allow_pickle=False)
        if latent.shape != (extractor.EXPECTED_TOKENS, extractor.EXPECTED_FEATURES):
            raise AuditError(f"Bad latent shape: {row['latent_path']} {latent.shape}")
        if latent.dtype != np.float32 or not np.isfinite(latent).all():
            raise AuditError(f"Latent must be finite float32: {row['latent_path']}")
        features[index] = latent.reshape(-1)
    return features


def load_model_inputs(rows: Sequence[dict[str, str]]) -> np.ndarray:
    runtime = extractor.import_runtime()
    load_resize, normalize = extractor.build_preprocessors(runtime)
    voxel_count = math.prod(extractor.IMAGE_SIZE)
    images = np.empty((len(rows), voxel_count), dtype=np.float32)
    for index, row in enumerate(rows):
        tensor = extractor.preprocess_volume(
            Path(row["nifti_path"]), load_resize, normalize, runtime
        )
        images[index] = np.asarray(tensor, dtype=np.float32).reshape(-1)
        if (index + 1) % 20 == 0 or index + 1 == len(rows):
            log(f"preprocessed {index + 1}/{len(rows)} model inputs")
    return images


def centered_kernel(
    train: np.ndarray, other: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray, float]:
    mean = train.mean(axis=0, dtype=np.float64).astype(np.float32)
    train_centered = train - mean
    scale = float(np.einsum("ij,ij->", train_centered, train_centered, dtype=np.float64))
    scale /= train_centered.shape[0]
    if not math.isfinite(scale) or scale <= 0:
        raise AuditError("Feature kernel has zero or invalid scale.")
    train_kernel = (train_centered @ train_centered.T).astype(np.float64) / scale
    if other is None:
        return train_kernel, mean, scale
    cross_kernel = ((other - mean) @ train_centered.T).astype(np.float64) / scale
    return train_kernel, cross_kernel, scale


def ridge_weights(kernel: np.ndarray, cross_kernel: np.ndarray, alpha: float) -> np.ndarray:
    regularized = kernel.copy()
    regularized.flat[:: regularized.shape[0] + 1] += alpha
    try:
        return np.linalg.solve(regularized, cross_kernel.T).T
    except np.linalg.LinAlgError as error:
        raise AuditError(f"Kernel-ridge solve failed at alpha={alpha:g}: {error}") from error


def squared_error_from_grams(
    weights: np.ndarray,
    train_image_gram: np.ndarray,
    validation_train_image_gram: np.ndarray,
    validation_squared_norm: float,
) -> float:
    cross_term = float(np.einsum("ij,ij->", weights, validation_train_image_gram))
    fitted_norm = float(np.einsum("ij,ij->", weights @ train_image_gram, weights))
    return validation_squared_norm - 2.0 * cross_term + fitted_norm


def tune_alpha(
    features: np.ndarray,
    images: np.ndarray,
    train_indices: np.ndarray,
    validation_indices: np.ndarray,
    alpha_ratios: Sequence[float],
) -> tuple[float, list[dict[str, float]]]:
    x_train = features[train_indices]
    x_validation = features[validation_indices]
    kernel, cross, _ = centered_kernel(x_train, x_validation)

    y_train = images[train_indices]
    y_validation = images[validation_indices]
    image_mean = y_train.mean(axis=0, dtype=np.float64).astype(np.float32)
    y_train_centered = y_train - image_mean
    y_validation_centered = y_validation - image_mean
    train_gram = (y_train_centered @ y_train_centered.T).astype(np.float64)
    validation_train_gram = (
        y_validation_centered @ y_train_centered.T
    ).astype(np.float64)
    validation_norm = float(
        np.einsum(
            "ij,ij->", y_validation_centered, y_validation_centered, dtype=np.float64
        )
    )
    denominator = y_validation_centered.size

    tuning: list[dict[str, float]] = []
    for alpha in alpha_ratios:
        weights = ridge_weights(kernel, cross, alpha)
        sse = max(
            0.0,
            squared_error_from_grams(
                weights, train_gram, validation_train_gram, validation_norm
            ),
        )
        tuning.append({"alpha": float(alpha), "validation_rmse": math.sqrt(sse / denominator)})
    best = min(tuning, key=lambda row: (row["validation_rmse"], row["alpha"]))
    return best["alpha"], tuning


def fit_predict(
    features: np.ndarray,
    images: np.ndarray,
    development_indices: np.ndarray,
    test_indices: np.ndarray,
    alpha: float,
) -> tuple[np.ndarray, np.ndarray]:
    x_development = features[development_indices]
    x_test = features[test_indices]
    kernel, cross, _ = centered_kernel(x_development, x_test)
    weights = ridge_weights(kernel, cross, alpha)
    y_development = images[development_indices]
    image_mean = y_development.mean(axis=0, dtype=np.float64).astype(np.float32)
    predictions = image_mean + weights.astype(np.float32) @ (y_development - image_mean)
    return predictions.astype(np.float32), image_mean


def safe_correlation(truth: np.ndarray, estimate: np.ndarray) -> float:
    truth_centered = truth - truth.mean(dtype=np.float64)
    estimate_centered = estimate - estimate.mean(dtype=np.float64)
    denominator = math.sqrt(
        float(np.dot(truth_centered, truth_centered))
        * float(np.dot(estimate_centered, estimate_centered))
    )
    return float(np.dot(truth_centered, estimate_centered) / denominator) if denominator else math.nan


def scan_metrics(truth: np.ndarray, estimate: np.ndarray) -> dict[str, float]:
    error = estimate - truth
    truth_norm = math.sqrt(float(np.dot(truth, truth)))
    foreground = truth != 0
    foreground_error = error[foreground]
    foreground_truth = truth[foreground]
    foreground_norm = math.sqrt(float(np.dot(foreground_truth, foreground_truth)))
    return {
        "rmse": math.sqrt(float(np.mean(error * error, dtype=np.float64))),
        "relative_rmse": math.sqrt(float(np.dot(error, error))) / truth_norm,
        "mae": float(np.mean(np.abs(error), dtype=np.float64)),
        "correlation": safe_correlation(truth, estimate),
        "foreground_rmse": math.sqrt(
            float(np.mean(foreground_error * foreground_error, dtype=np.float64))
        ),
        "foreground_relative_rmse": math.sqrt(float(np.dot(foreground_error, foreground_error)))
        / foreground_norm,
    }


def summarize(values: Sequence[dict[str, Any]], keys: Iterable[str]) -> dict[str, float]:
    result: dict[str, float] = {}
    for key in keys:
        column = np.asarray([row[key] for row in values], dtype=np.float64)
        result[f"mean_{key}"] = float(np.nanmean(column))
        result[f"median_{key}"] = float(np.nanmedian(column))
    return result


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise AuditError(f"Refusing to write an empty CSV: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def grayscale(array: np.ndarray, low: float, high: float) -> Any:
    from PIL import Image

    scaled = np.clip((array - low) / (high - low), 0.0, 1.0)
    return Image.fromarray(np.asarray(np.rint(255 * scaled), dtype=np.uint8), mode="L").convert("RGB")


def error_rgb(array: np.ndarray, limit: float) -> Any:
    from PIL import Image

    scaled = np.clip(array / limit, -1.0, 1.0)
    red = np.where(scaled >= 0, 255, 255 * (1 + scaled))
    blue = np.where(scaled <= 0, 255, 255 * (1 - scaled))
    green = 255 * (1 - np.abs(scaled))
    rgb = np.stack((red, green, blue), axis=-1)
    return Image.fromarray(np.asarray(np.rint(rgb), dtype=np.uint8), mode="RGB")


def orthogonal_slices(volume: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    center = tuple(size // 2 for size in volume.shape)
    return (
        np.rot90(volume[center[0], :, :]),
        np.rot90(volume[:, center[1], :]),
        np.rot90(volume[:, :, center[2]]),
    )


def draw_examples(
    path: Path,
    truth: np.ndarray,
    prediction: np.ndarray,
    mean_image: np.ndarray,
    selected: Sequence[int],
    labels: Sequence[str],
    layer_name: str,
) -> None:
    from PIL import Image, ImageDraw, ImageFont

    font = ImageFont.load_default()
    scale = 2
    tile = extractor.IMAGE_SIZE[0] * scale
    left = 86
    top = 44
    gap = 6
    column_labels = [
        f"truth {axis}" for axis in ("axis 0", "axis 1", "axis 2")
    ] + [f"mean {axis}" for axis in ("axis 0", "axis 1", "axis 2")] + [
        f"probe {axis}" for axis in ("axis 0", "axis 1", "axis 2")
    ] + [
        f"error {axis}" for axis in ("axis 0", "axis 1", "axis 2")
    ]
    width = left + len(column_labels) * (tile + gap) + gap
    height = top + len(selected) * (tile + 28) + gap
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    draw.text((8, 6), f"{layer_name}: held-out linear-probe reconstruction", fill="black", font=font)
    draw.text((8, 20), "gray: normalized model input/probe; error: blue = low, red = high", fill="black", font=font)
    for column, label in enumerate(column_labels):
        draw.text((left + column * (tile + gap), top - 14), label, fill="black", font=font)

    for row_position, selected_position in enumerate(selected):
        y = top + row_position * (tile + 28)
        draw.text((8, y + 4), labels[row_position], fill="black", font=font)
        truth_volume = truth[selected_position].reshape(extractor.IMAGE_SIZE)
        prediction_volume = prediction[selected_position].reshape(extractor.IMAGE_SIZE)
        mean_volume = mean_image.reshape(extractor.IMAGE_SIZE)
        panels = [grayscale(value, -2.0, 3.0) for value in orthogonal_slices(truth_volume)]
        panels += [grayscale(value, -2.0, 3.0) for value in orthogonal_slices(mean_volume)]
        panels += [grayscale(value, -2.0, 3.0) for value in orthogonal_slices(prediction_volume)]
        panels += [error_rgb(value, 2.0) for value in orthogonal_slices(prediction_volume - truth_volume)]
        for column, panel in enumerate(panels):
            panel = panel.resize((tile, tile), resample=Image.Resampling.NEAREST)
            canvas.paste(panel, (left + column * (tile + gap), y))
        metrics = scan_metrics(truth[selected_position], prediction[selected_position])
        draw.text(
            (8, y + 20),
            f"rRMSE {metrics['relative_rmse']:.3f}\nr {metrics['correlation']:.3f}",
            fill="black",
            font=font,
        )
    canvas.save(path, optimize=True)


def parse_layer(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("--layer must have NAME=METADATA.csv form")
    name, raw_path = value.split("=", 1)
    name = name.strip()
    if not name or not raw_path.strip():
        raise argparse.ArgumentTypeError("--layer must have NAME=METADATA.csv form")
    return name, Path(raw_path).expanduser().resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--layer", action="append", type=parse_layer,
        help="Layer and extraction metadata as NAME=PATH; repeat for several layers.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--test-fraction", type=float, default=0.20)
    parser.add_argument("--validation-fraction", type=float, default=0.20)
    parser.add_argument(
        "--alpha-ratios", type=float, nargs="+", default=list(DEFAULT_ALPHA_RATIOS)
    )
    parser.add_argument("--examples", type=int, default=3)
    args = parser.parse_args()
    if not 0 < args.test_fraction < 0.5 or not 0 < args.validation_fraction < 0.5:
        parser.error("split fractions must lie strictly between 0 and 0.5")
    if args.examples < 1 or any(value <= 0 or not math.isfinite(value) for value in args.alpha_ratios):
        parser.error("examples and alpha ratios must be positive")
    return args


def main() -> int:
    args = parse_args()
    project_root = Path(__file__).resolve().parents[2]
    layer_args = args.layer or [
        (name, (project_root / relative).resolve()) for name, relative in DEFAULT_LAYER_PATHS
    ]
    if len({name for name, _ in layer_args}) != len(layer_args):
        raise AuditError("Layer names must be unique.")
    layers = [(name, path, read_metadata(path)) for name, path in layer_args]
    validate_layer_contracts(layers)
    rows = layers[0][2]

    train_subjects, validation_subjects, test_subjects = make_subject_split(
        rows, args.seed, args.test_fraction, args.validation_fraction
    )
    development_subjects = train_subjects | validation_subjects
    train_indices = indices_for(rows, train_subjects)
    validation_indices = indices_for(rows, validation_subjects)
    development_indices = indices_for(rows, development_subjects)
    test_indices = indices_for(rows, test_subjects)
    log(
        f"subject split: {len(train_subjects)} tune-train, {len(validation_subjects)} "
        f"tune-validation, {len(test_subjects)} final-test"
    )

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    log("loading and reproducing normalized 96^3 BrainIAC inputs")
    images = load_model_inputs(rows)

    metric_keys = (
        "rmse", "relative_rmse", "mae", "correlation", "foreground_rmse",
        "foreground_relative_rmse", "relative_rmse_reduction_vs_mean_percent",
    )
    all_metrics: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    tuning_rows: list[dict[str, Any]] = []
    split_rows = [
        {
            "anonymous_subject": f"subject_{position + 1:02d}",
            "split": "tune_train" if subject in train_subjects else (
                "tune_validation" if subject in validation_subjects else "final_test"
            ),
        }
        for position, subject in enumerate(sorted(train_subjects | validation_subjects | test_subjects))
    ]
    write_csv(output_dir / "subject_split_anonymous.csv", split_rows)

    final_predictions: dict[str, np.ndarray] = {}
    test_truth = images[test_indices]
    for layer_name, metadata_path, layer_rows in layers:
        log(f"loading {layer_name} features")
        features = load_features(layer_rows)
        log(f"tuning subject-held-out kernel ridge for {layer_name}")
        best_alpha, tuning = tune_alpha(
            features, images, train_indices, validation_indices, args.alpha_ratios
        )
        for result in tuning:
            tuning_rows.append({"layer": layer_name, **result})
        log(f"{layer_name}: selected alpha={best_alpha:g}")
        prediction, mean_image = fit_predict(
            features, images, development_indices, test_indices, best_alpha
        )
        final_predictions[layer_name] = prediction

        probe_rows: list[dict[str, Any]] = []
        baseline_rows: list[dict[str, Any]] = []
        for position, visit_index in enumerate(test_indices):
            anonymous_visit = f"test_visit_{position + 1:02d}"
            probe_metric = scan_metrics(test_truth[position], prediction[position])
            baseline_metric = scan_metrics(test_truth[position], mean_image)
            probe_metric["relative_rmse_reduction_vs_mean_percent"] = 100.0 * (
                baseline_metric["relative_rmse"] - probe_metric["relative_rmse"]
            ) / baseline_metric["relative_rmse"]
            baseline_metric["relative_rmse_reduction_vs_mean_percent"] = 0.0
            probe_row = {
                "layer": layer_name, "method": "linear_probe",
                "anonymous_visit": anonymous_visit, **probe_metric,
            }
            baseline_row = {
                "layer": layer_name, "method": "training_mean_baseline",
                "anonymous_visit": anonymous_visit, **baseline_metric,
            }
            probe_rows.append(probe_row)
            baseline_rows.append(baseline_row)
            all_metrics.extend((probe_row, baseline_row))
        summaries.append(
            {
                "layer": layer_name, "method": "linear_probe", "selected_alpha": best_alpha,
                "test_subjects": len(test_subjects), "test_visits": len(test_indices),
                **summarize(probe_rows, metric_keys),
            }
        )
        summaries.append(
            {
                "layer": layer_name, "method": "training_mean_baseline", "selected_alpha": "",
                "test_subjects": len(test_subjects), "test_visits": len(test_indices),
                **summarize(baseline_rows, metric_keys),
            }
        )

        errors = np.asarray([row["relative_rmse"] for row in probe_rows])
        order = np.argsort(errors)
        quantile_positions = np.linspace(0, len(order) - 1, min(args.examples, len(order)))
        selected = [int(order[round(value)]) for value in quantile_positions]
        labels = [f"example {position + 1}" for position in range(len(selected))]
        draw_examples(
            output_dir / f"reconstruction_examples_{layer_name}.png",
            test_truth, prediction, mean_image, selected, labels, layer_name,
        )
        del features, prediction

    write_csv(output_dir / "reconstruction_metrics.csv", all_metrics)
    write_csv(output_dir / "reconstruction_summary.csv", summaries)
    write_csv(output_dir / "ridge_tuning.csv", tuning_rows)
    configuration = {
        "script_version": SCRIPT_VERSION,
        "script_sha256": sha256_file(Path(__file__)),
        "seed": args.seed,
        "test_fraction": args.test_fraction,
        "validation_fraction": args.validation_fraction,
        "alpha_ratios": args.alpha_ratios,
        "layers": [
            {
                "name": name,
                "metadata": str(path),
                "encoder_layer": rows_for_layer[0]["encoder_layer"],
                "checkpoint_sha256": rows_for_layer[0]["checkpoint_sha256"],
            }
            for name, path, rows_for_layer in layers
        ],
        "subjects": {
            "tune_train": len(train_subjects),
            "tune_validation": len(validation_subjects),
            "final_test": len(test_subjects),
        },
        "visits": {
            "tune_train": len(train_indices),
            "tune_validation": len(validation_indices),
            "development_refit": len(development_indices),
            "final_test": len(test_indices),
        },
        "target": "MONAI-resized and nonzero-z-normalized 96x96x96 model input",
        "interpretation": (
            "Subject-held-out linear decodability audit; not an official BrainIAC decoder, "
            "native-space reconstruction, or validation of off-manifold MatFPCA latents."
        ),
    }
    with (output_dir / "run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(configuration, handle, indent=2)
        handle.write("\n")
    log(f"wrote reconstruction audit to {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AuditError, extractor.ExtractionError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
