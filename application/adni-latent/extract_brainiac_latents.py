#!/usr/bin/env python3
"""Extract longitudinal ADNI BrainIAC patch-token matrices.

This script intentionally separates frozen representation extraction from the
R MatFPCA analysis.  It follows the deterministic validation transform in the
official BrainIAC repository:

    LoadImaged -> EnsureChannelFirstd -> fixed common foreground crop
    -> Resized(96^3, trilinear)
    -> NormalizeIntensityd(nonzero=True, channel_wise=True) -> ToTensord

BrainIAC uses ``monai.networks.nets.ViT`` with ``classification=False``.  There
is consequently no CLS token: a 96^3 input and 16^3 patches produce all 216
patch tokens, each with 768 features.  The upstream feature-extraction wrapper
calls token zero a CLS token, but it is the first spatial patch token.

The command is safe to resume.  An existing latent is skipped only when the
run configuration matches and the file is a finite float32 array with the
expected shape.  Writes are atomic.  No encoder is trained or fine-tuned.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import datetime as dt
import hashlib
import json
import logging
import math
import os
from pathlib import Path
import random
import re
import subprocess
import sys
import tempfile
from typing import Any, Iterable, Mapping, Sequence


SCHEMA_VERSION = 1
SCRIPT_VERSION = "1.1.0"
MODEL_NAME = "BrainIAC"
MODEL_ARCHITECTURE = "MONAI ViT-B/16^3 (3D, classification=False)"
# These affine-registered ADNI inputs share one grid.  Use one common crop:
# independently cropping each scan would move anatomy between token indices.
SOURCE_IMAGE_SIZE = (193, 229, 193)
FOREGROUND_ROI_START = (16, 16, 0)
FOREGROUND_ROI_END = (177, 217, 177)  # Python-exclusive upper bounds.
CROPPED_IMAGE_SIZE = tuple(
    end - start for start, end in zip(FOREGROUND_ROI_START, FOREGROUND_ROI_END)
)
IMAGE_SIZE = (96, 96, 96)
PATCH_SIZE = (16, 16, 16)
PATCH_GRID = (6, 6, 6)
EXPECTED_TOKENS = 216
EXPECTED_FEATURES = 768
NUM_LAYERS = 12
MLP_DIM = 3072
NUM_HEADS = 12
REPRESENTATION_TYPE = "3d_patch_token_by_feature"
PREPROCESSING_DESCRIPTION = (
    "MONAI LoadImaged; EnsureChannelFirstd; fixed common SpatialCropd "
    "roi_start=(16,16,0) roi_end=(177,217,177) on the affine-registered "
    "193x229x193 grid; Resized spatial_size=(96,96,96) mode=trilinear; "
    "NormalizeIntensityd nonzero=True channel_wise=True; ToTensord"
)

REQUIRED_MANIFEST_COLUMNS = {
    "PTID",
    "image_id",
    "image_time_month",
    "mri_VISCODE2",
    "image_retained_for_analysis",
}

DERIVED_METADATA_COLUMNS = [
    "subject_id",
    "visit_id",
    "time_month",
    "viscode",
    "time_year",
    "nifti_path",
    "nifti_size_bytes",
    "nifti_mtime_ns",
    "latent_path",
    "latent_sha256",
    "extraction_status",
    "encoder_name",
    "encoder_architecture",
    "checkpoint_sha256",
    "representation_type",
    "encoder_layer",
    "patch_grid",
    "latent_rows",
    "latent_cols",
    "latent_dtype",
    "preprocessing",
]

TRUE_VALUES = {"1", "true", "t", "yes", "y"}
FALSE_VALUES = {"", "0", "false", "f", "no", "n", "na", "nan", "null"}
NIFTI_PATTERN = re.compile(r"_I([0-9]+)(?:_|(?=[.]nii(?:[.]gz)?$))", re.IGNORECASE)
SAFE_PATH_COMPONENT = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")

LOGGER = logging.getLogger("brainiac_latents")


class ExtractionError(RuntimeError):
    """A user-actionable extraction or validation failure."""


@dataclasses.dataclass(frozen=True)
class Visit:
    """One retained longitudinal scan and its original manifest row."""

    original: dict[str, str]
    source_row: int
    subject_id: str
    image_id: int
    visit_code: str
    image_time_month: float
    model_time_month: float
    time_year: float
    nifti_path: Path | None = None
    nifti_size_bytes: int | None = None
    nifti_mtime_ns: int | None = None


@dataclasses.dataclass(frozen=True)
class Cohort:
    fieldnames: tuple[str, ...]
    visits: tuple[Visit, ...]
    manifest_rows: int
    retained_within_window: int
    subjects_before_min_visits: int
    eligible_subjects: int
    selected_subjects: tuple[str, ...]


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def sha256_json(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def analysis_time_tick(model_time_month: float) -> int:
    """Match the day-resolution convention used by adni_latent_matfpca.R."""

    return round(model_time_month * 365.25 / 12.0)


def parse_retained(value: str | None, row_number: int) -> bool:
    normalized = "" if value is None else value.strip().lower()
    if normalized in TRUE_VALUES:
        return True
    if normalized in FALSE_VALUES:
        return False
    raise ExtractionError(
        "Manifest row "
        f"{row_number}: image_retained_for_analysis has unrecognized value "
        f"{value!r}; expected TRUE/FALSE."
    )


def parse_finite_float(value: str | None, field: str, row_number: int) -> float:
    try:
        parsed = float(value) if value is not None else math.nan
    except (TypeError, ValueError) as error:
        raise ExtractionError(
            f"Manifest row {row_number}: {field}={value!r} is not numeric."
        ) from error
    if not math.isfinite(parsed):
        raise ExtractionError(
            f"Manifest row {row_number}: {field}={value!r} is not finite."
        )
    return parsed


def parse_image_id(value: str | None, row_number: int) -> int:
    normalized = "" if value is None else value.strip()
    if not re.fullmatch(r"[0-9]+", normalized):
        raise ExtractionError(
            f"Manifest row {row_number}: image_id={value!r} is not an integer."
        )
    return int(normalized)


def validate_subject_id(subject_id: str, row_number: int) -> None:
    if not SAFE_PATH_COMPONENT.fullmatch(subject_id):
        raise ExtractionError(
            f"Manifest row {row_number}: PTID={subject_id!r} is not a safe "
            "directory name. Expected only letters, digits, '.', '_' or '-'."
        )


def read_cohort(
    manifest_path: Path,
    *,
    max_months: float,
    min_visits: int,
    max_subjects: int | None,
    seed: int,
) -> Cohort:
    """Read and deterministically select a cohort using only stdlib CSV."""

    if not manifest_path.is_file():
        raise ExtractionError(f"Manifest does not exist: {manifest_path}")
    if not math.isfinite(max_months):
        raise ExtractionError("--max-months must be finite.")
    if min_visits < 1:
        raise ExtractionError("--min-visits must be at least one.")
    if max_subjects is not None and max_subjects < 1:
        raise ExtractionError("--max-subjects must be positive or 'all'.")

    visits: list[Visit] = []
    manifest_rows = 0
    with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ExtractionError(f"Manifest has no header: {manifest_path}")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ExtractionError("Manifest contains duplicate column names.")
        missing = sorted(REQUIRED_MANIFEST_COLUMNS - set(reader.fieldnames))
        if missing:
            raise ExtractionError(
                "Manifest is missing required columns: " + ", ".join(missing)
            )
        collisions = sorted(set(reader.fieldnames) & set(DERIVED_METADATA_COLUMNS))
        if collisions:
            raise ExtractionError(
                "Manifest already contains reserved output columns: "
                + ", ".join(collisions)
            )

        for row_number, raw_row in enumerate(reader, start=2):
            manifest_rows += 1
            if None in raw_row:
                raise ExtractionError(
                    f"Manifest row {row_number} has more fields than the header."
                )
            row = {key: "" if value is None else value for key, value in raw_row.items()}
            if not parse_retained(row.get("image_retained_for_analysis"), row_number):
                continue

            image_time_month = parse_finite_float(
                row.get("image_time_month"), "image_time_month", row_number
            )
            if image_time_month > max_months:
                continue

            subject_id = row.get("PTID", "").strip()
            if not subject_id:
                raise ExtractionError(f"Manifest row {row_number}: PTID is empty.")
            validate_subject_id(subject_id, row_number)
            image_id = parse_image_id(row.get("image_id"), row_number)
            visit_code = row.get("mri_VISCODE2", "").strip()
            model_time_month = 0.0 if visit_code.lower() == "bl" else image_time_month
            visits.append(
                Visit(
                    original=row,
                    source_row=row_number,
                    subject_id=subject_id,
                    image_id=image_id,
                    visit_code=visit_code,
                    image_time_month=image_time_month,
                    model_time_month=model_time_month,
                    time_year=model_time_month / 12.0,
                )
            )

        fieldnames = tuple(reader.fieldnames)

    if not visits:
        raise ExtractionError(
            "No retained visits remain at or before --max-months."
        )

    image_rows: dict[int, int] = {}
    duplicate_images: list[tuple[int, int, int]] = []
    for visit in visits:
        previous = image_rows.get(visit.image_id)
        if previous is not None:
            duplicate_images.append((visit.image_id, previous, visit.source_row))
        else:
            image_rows[visit.image_id] = visit.source_row
    if duplicate_images:
        examples = ", ".join(
            f"I{image_id} (rows {first}/{second})"
            for image_id, first, second in duplicate_images[:10]
        )
        raise ExtractionError(
            f"Retained manifest image_id values are not unique: {examples}"
        )

    by_subject: dict[str, list[Visit]] = {}
    for visit in visits:
        by_subject.setdefault(visit.subject_id, []).append(visit)
    subjects_before_min_visits = len(by_subject)
    eligible = {
        subject: subject_visits
        for subject, subject_visits in by_subject.items()
        if len(subject_visits) >= min_visits
    }
    if not eligible:
        raise ExtractionError(
            f"No subjects have at least {min_visits} retained visits."
        )

    eligible_subject_ids = sorted(eligible)
    if max_subjects is not None and max_subjects < len(eligible_subject_ids):
        rng = random.Random(seed)
        selected_subjects = sorted(rng.sample(eligible_subject_ids, max_subjects))
    else:
        selected_subjects = eligible_subject_ids

    selected_set = set(selected_subjects)
    selected_visits = sorted(
        (visit for visit in visits if visit.subject_id in selected_set),
        key=lambda visit: (visit.subject_id, visit.time_year, visit.image_id),
    )

    subject_times: dict[tuple[str, int], Visit] = {}
    collisions: list[tuple[Visit, Visit]] = []
    for visit in selected_visits:
        key = (visit.subject_id, analysis_time_tick(visit.model_time_month))
        previous = subject_times.get(key)
        if previous is not None:
            collisions.append((previous, visit))
        else:
            subject_times[key] = visit
    if collisions:
        examples = ", ".join(
            f"{first.subject_id} on analysis day "
            f"{analysis_time_tick(first.model_time_month)} "
            f"(I{first.image_id}/I{second.image_id})"
            for first, second in collisions[:10]
        )
        raise ExtractionError(
            "Selected visits contain duplicate subject/day locations after "
            f"forcing baseline to zero and rounding to analysis days: {examples}"
        )

    return Cohort(
        fieldnames=fieldnames,
        visits=tuple(selected_visits),
        manifest_rows=manifest_rows,
        retained_within_window=len(visits),
        subjects_before_min_visits=subjects_before_min_visits,
        eligible_subjects=len(eligible),
        selected_subjects=tuple(selected_subjects),
    )


def is_nifti(path: Path) -> bool:
    lowered = path.name.lower()
    return lowered.endswith(".nii") or lowered.endswith(".nii.gz")


def index_nifti_files(nifti_dir: Path, *, recursive: bool) -> dict[int, Path]:
    if not nifti_dir.is_dir():
        raise ExtractionError(f"NIfTI directory does not exist: {nifti_dir}")
    paths = nifti_dir.rglob("*") if recursive else nifti_dir.iterdir()
    matches: dict[int, Path] = {}
    duplicates: dict[int, list[Path]] = {}
    for path in paths:
        if not path.is_file() or not is_nifti(path):
            continue
        match = NIFTI_PATTERN.search(path.name)
        if match is None:
            continue
        image_id = int(match.group(1))
        resolved = path.resolve()
        if image_id in matches:
            duplicates.setdefault(image_id, [matches[image_id]]).append(resolved)
        else:
            matches[image_id] = resolved
    if duplicates:
        examples = "; ".join(
            f"I{image_id}: " + ", ".join(str(path) for path in paths[:3])
            for image_id, paths in sorted(duplicates.items())[:10]
        )
        raise ExtractionError(
            "NIfTI filenames do not provide unique _I<image_id>_ matches: "
            + examples
        )
    if not matches:
        raise ExtractionError(
            f"No .nii/.nii.gz filenames with _I<image_id>_ were found in {nifti_dir}."
        )
    return matches


def attach_nifti_paths(cohort: Cohort, nifti_index: Mapping[int, Path]) -> Cohort:
    missing = [visit for visit in cohort.visits if visit.image_id not in nifti_index]
    if missing:
        examples = ", ".join(
            f"I{visit.image_id} ({visit.subject_id})" for visit in missing[:20]
        )
        raise ExtractionError(
            f"{len(missing)} selected visits have no unique NIfTI match: {examples}"
        )
    attached: list[Visit] = []
    for visit in cohort.visits:
        path = nifti_index[visit.image_id]
        stat = path.stat()
        if stat.st_size <= 0:
            raise ExtractionError(f"NIfTI file is empty: {path}")
        attached.append(
            dataclasses.replace(
                visit,
                nifti_path=path,
                nifti_size_bytes=stat.st_size,
                nifti_mtime_ns=stat.st_mtime_ns,
            )
        )
    return dataclasses.replace(cohort, visits=tuple(attached))


def parse_max_subjects(value: str) -> int | None:
    if value.strip().lower() in {"all", "none", "inf", "infinity"}:
        return None
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "must be a positive integer or 'all'"
        ) from error
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be positive or 'all'")
    return parsed


def parse_layer(value: str) -> str | int:
    normalized = value.strip().lower()
    if normalized == "final":
        return "final"
    try:
        layer = int(normalized)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "must be 'final' or an integer from 1 through 12"
        ) from error
    if not 1 <= layer <= NUM_LAYERS:
        raise argparse.ArgumentTypeError("must lie between 1 and 12")
    return layer


def parse_seed(value: str) -> int:
    try:
        seed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be an integer") from error
    if not 0 <= seed <= 2**32 - 1:
        raise argparse.ArgumentTypeError("must lie between 0 and 2^32 - 1")
    return seed


def layer_label(layer: str | int) -> str:
    if layer == "final":
        return "final_normalized"
    return f"block_{layer:02d}_pre_final_norm"


def latent_relative_path(visit: Visit) -> Path:
    return Path("latents") / visit.subject_id / f"{visit.image_id}.npy"


def cohort_fingerprint(cohort: Cohort) -> str:
    rows = [
        {
            "subject_id": visit.subject_id,
            "image_id": visit.image_id,
            "model_time_month": visit.model_time_month,
            "nifti_path": str(visit.nifti_path),
            "nifti_size_bytes": visit.nifti_size_bytes,
            "nifti_mtime_ns": visit.nifti_mtime_ns,
        }
        for visit in cohort.visits
    ]
    return sha256_json(rows)


def git_provenance(repo_root: Path) -> dict[str, Any]:
    def git(*arguments: str) -> str | None:
        try:
            result = subprocess.run(
                ["git", *arguments],
                cwd=repo_root,
                check=False,
                capture_output=True,
                text=True,
                timeout=10,
            )
        except (OSError, subprocess.SubprocessError):
            return None
        if result.returncode != 0:
            return None
        return result.stdout.strip()

    commit = git("rev-parse", "HEAD")
    status = git("status", "--porcelain", "--untracked-files=all")
    return {
        "repository_root": str(repo_root),
        "commit": commit,
        "dirty": None if status is None else bool(status),
        "status_porcelain": None if status is None else status.splitlines(),
    }


def atomic_write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def atomic_write_metadata(
    path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def quarantine_existing_metadata(metadata_path: Path) -> Path | None:
    """Hide stale metadata before any latent can be replaced in a new run."""

    if not metadata_path.exists():
        return None
    if not metadata_path.is_file():
        raise ExtractionError(
            f"Expected metadata path to be a regular file: {metadata_path}"
        )
    previous_path = metadata_path.with_name("metadata.previous.csv")
    try:
        os.replace(metadata_path, previous_path)
    except OSError as error:
        raise ExtractionError(
            f"Could not quarantine existing metadata before extraction: {error}"
        ) from error
    return previous_path


def atomic_save_npy(path: Path, array: Any, np_module: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w+b",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            np_module.save(handle, array, allow_pickle=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def validate_cached_latent(
    path: Path, np_module: Any, expected_shape: tuple[int, int]
) -> tuple[bool, str]:
    if not path.is_file():
        return False, "missing"
    try:
        array = np_module.load(path, mmap_mode="r", allow_pickle=False)
    except Exception as error:  # NumPy raises several format/IO exception types.
        return False, f"unreadable: {error}"
    if tuple(array.shape) != expected_shape:
        return False, f"shape {tuple(array.shape)} != {expected_shape}"
    if array.dtype != np_module.dtype("float32"):
        return False, f"dtype {array.dtype} != float32"
    if not bool(np_module.isfinite(array).all()):
        return False, "contains non-finite values"
    return True, "valid"


def import_runtime() -> dict[str, Any]:
    try:
        import numpy as np
        import nibabel
        import monai
        import torch
        from monai.networks.nets import ViT
        from monai.transforms import (
            Compose,
            EnsureChannelFirstd,
            LoadImaged,
            NormalizeIntensityd,
            Resized,
            SpatialCropd,
            ToTensord,
        )
    except ImportError as error:
        raise ExtractionError(
            "Missing BrainIAC extraction dependency. Install "
            "application/adni-latent/requirements-brainiac.txt in a dedicated "
            f"environment. Original import error: {error}"
        ) from error

    return {
        "np": np,
        "nibabel": nibabel,
        "monai": monai,
        "torch": torch,
        "ViT": ViT,
        "Compose": Compose,
        "LoadImaged": LoadImaged,
        "EnsureChannelFirstd": EnsureChannelFirstd,
        "Resized": Resized,
        "SpatialCropd": SpatialCropd,
        "NormalizeIntensityd": NormalizeIntensityd,
        "ToTensord": ToTensord,
    }


def runtime_versions(runtime: Mapping[str, Any]) -> dict[str, str]:
    return {
        "python": sys.version.split()[0],
        "numpy": str(runtime["np"].__version__),
        "torch": str(runtime["torch"].__version__),
        "monai": str(runtime["monai"].__version__),
        "nibabel": str(runtime["nibabel"].__version__),
    }


def configure_determinism(runtime: Mapping[str, Any], seed: int) -> None:
    np = runtime["np"]
    torch = runtime["torch"]
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    if hasattr(torch.backends, "cudnn"):
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
    torch.use_deterministic_algorithms(True)


def resolve_device(requested: str, torch: Any) -> Any:
    normalized = requested.strip().lower()
    if normalized == "auto":
        normalized = "cuda" if torch.cuda.is_available() else "cpu"
    if normalized == "cpu":
        return torch.device("cpu")
    if normalized == "cuda":
        normalized = "cuda:0"
    if not re.fullmatch(r"cuda:[0-9]+", normalized):
        raise ExtractionError(
            "--device must be 'auto', 'cpu', 'cuda', or 'cuda:<index>'."
        )
    if not torch.cuda.is_available():
        raise ExtractionError(
            f"CUDA device {normalized!r} was requested, but torch reports no CUDA."
        )
    index = int(normalized.split(":", 1)[1])
    if index >= torch.cuda.device_count():
        raise ExtractionError(
            f"CUDA device index {index} is unavailable; found "
            f"{torch.cuda.device_count()} device(s)."
        )
    return torch.device(normalized)


def create_model(runtime: Mapping[str, Any]) -> Any:
    model = runtime["ViT"](
        in_channels=1,
        img_size=IMAGE_SIZE,
        patch_size=PATCH_SIZE,
        hidden_size=EXPECTED_FEATURES,
        mlp_dim=MLP_DIM,
        num_layers=NUM_LAYERS,
        num_heads=NUM_HEADS,
        classification=False,
        dropout_rate=0.0,
        spatial_dims=3,
        save_attn=True,
    )
    if hasattr(model, "cls_token"):
        raise ExtractionError(
            "Constructed MONAI ViT unexpectedly has a CLS token; refusing to "
            "silently change the 216-patch representation."
        )
    return model


def load_checkpoint_mapping(
    checkpoint_path: Path, runtime: Mapping[str, Any]
) -> Mapping[str, Any]:
    torch = runtime["torch"]
    if checkpoint_path.suffix.lower() == ".safetensors":
        try:
            from safetensors.torch import load_file
        except ImportError as error:
            raise ExtractionError(
                "A .safetensors checkpoint requires the safetensors package."
            ) from error
        loaded: Any = load_file(str(checkpoint_path), device="cpu")
    else:
        LOGGER.warning(
            "Loading a pickle-based PyTorch checkpoint. Only use a checkpoint "
            "obtained from a trusted BrainIAC source."
        )
        try:
            loaded = torch.load(
                str(checkpoint_path), map_location="cpu", weights_only=False
            )
        except TypeError:
            loaded = torch.load(str(checkpoint_path), map_location="cpu")

    if not isinstance(loaded, Mapping):
        raise ExtractionError(
            f"Checkpoint root is {type(loaded).__name__}, not a state dictionary."
        )

    candidates: list[Mapping[str, Any]] = [loaded]
    for key in ("state_dict", "model_state_dict", "model", "backbone"):
        nested = loaded.get(key)
        if isinstance(nested, Mapping):
            candidates.insert(0, nested)
    for candidate in candidates:
        tensor_items = {
            str(key): value
            for key, value in candidate.items()
            if torch.is_tensor(value)
        }
        if tensor_items:
            return tensor_items
    raise ExtractionError("Checkpoint does not contain a tensor state dictionary.")


def load_model_weights(
    model: Any, raw_state: Mapping[str, Any]
) -> dict[str, Any]:
    expected = model.state_dict()
    expected_keys = set(expected)
    prefixes = {
        "",
        "backbone.",
        "model.",
        "module.",
        "model.backbone.",
        "module.backbone.",
        "backbone.backbone.",
        "module.model.backbone.",
    }
    anchor_keys = list(expected_keys)[: min(12, len(expected_keys))]
    for raw_key in raw_state:
        for expected_key in anchor_keys:
            if raw_key.endswith(expected_key):
                prefixes.add(raw_key[: -len(expected_key)])

    attempts: list[tuple[int, int, str, dict[str, Any], list[str]]] = []
    for prefix in prefixes:
        stripped = {
            key[len(prefix) :]: value
            for key, value in raw_state.items()
            if key.startswith(prefix)
        }
        matched = {key: value for key, value in stripped.items() if key in expected_keys}
        shape_errors = [
            key
            for key, value in matched.items()
            if tuple(value.shape) != tuple(expected[key].shape)
        ]
        attempts.append((len(matched), -len(shape_errors), prefix, matched, shape_errors))

    _, _, prefix, matched, shape_errors = max(
        attempts, key=lambda item: (item[0], item[1], -len(item[2]))
    )
    missing = sorted(expected_keys - set(matched))
    if missing or shape_errors:
        details = []
        if missing:
            details.append(
                f"missing {len(missing)} keys (examples: {', '.join(missing[:8])})"
            )
        if shape_errors:
            details.append(
                "shape mismatch for " + ", ".join(sorted(shape_errors)[:8])
            )
        raise ExtractionError(
            "Checkpoint is not an exact BrainIAC MONAI ViT backbone after best "
            f"prefix {prefix!r}: " + "; ".join(details)
        )

    model.load_state_dict(matched, strict=True)
    used_raw_keys = {
        prefix + key for key in matched if prefix + key in raw_state
    }
    ignored = sorted(set(raw_state) - used_raw_keys)
    return {
        "selected_prefix": prefix,
        "raw_tensor_keys": len(raw_state),
        "loaded_tensor_keys": len(matched),
        "ignored_tensor_keys": len(ignored),
        "ignored_key_examples": ignored[:20],
    }


def build_preprocessors(runtime: Mapping[str, Any]) -> tuple[Any, Any]:
    load_resize = runtime["Compose"](
        [
            runtime["LoadImaged"](keys=["image"]),
            runtime["EnsureChannelFirstd"](keys=["image"]),
            runtime["SpatialCropd"](
                keys=["image"],
                roi_start=FOREGROUND_ROI_START,
                roi_end=FOREGROUND_ROI_END,
            ),
            runtime["Resized"](
                keys=["image"], spatial_size=IMAGE_SIZE, mode="trilinear"
            ),
        ]
    )
    normalize = runtime["Compose"](
        [
            runtime["NormalizeIntensityd"](
                keys=["image"], nonzero=True, channel_wise=True
            ),
            runtime["ToTensord"](keys=["image"]),
        ]
    )
    return load_resize, normalize


def preprocess_volume(
    path: Path, load_resize: Any, normalize: Any, runtime: Mapping[str, Any]
) -> Any:
    torch = runtime["torch"]
    try:
        source_shape = tuple(runtime["nibabel"].load(str(path)).shape[:3])
    except Exception as error:
        raise ExtractionError(f"Failed to inspect NIfTI geometry {path}: {error}") from error
    if source_shape != SOURCE_IMAGE_SIZE:
        raise ExtractionError(
            f"Source {path} has spatial shape {source_shape}; the fixed common "
            f"crop requires the registered grid {SOURCE_IMAGE_SIZE}."
        )
    try:
        resized = load_resize({"image": str(path)})
    except Exception as error:
        raise ExtractionError(f"Failed to load/resize NIfTI {path}: {error}") from error

    resized_tensor = torch.as_tensor(resized["image"])
    if tuple(resized_tensor.shape) != (1, *IMAGE_SIZE):
        raise ExtractionError(
            f"Preprocessed {path} has shape {tuple(resized_tensor.shape)}; "
            f"expected {(1, *IMAGE_SIZE)}."
        )
    if not bool(torch.isfinite(resized_tensor).all()):
        raise ExtractionError(f"Resized volume contains non-finite values: {path}")
    nonzero = resized_tensor != 0
    nonzero_count = int(nonzero.sum().item())
    if nonzero_count < 2:
        raise ExtractionError(
            f"Resized volume has fewer than two nonzero voxels: {path}"
        )
    nonzero_values = resized_tensor[nonzero].float()
    if float(nonzero_values.std(unbiased=False).item()) <= 1e-8:
        raise ExtractionError(
            f"Nonzero voxels have zero or near-zero variance: {path}"
        )

    try:
        output = normalize(resized)["image"]
    except Exception as error:
        raise ExtractionError(f"Failed to normalize NIfTI {path}: {error}") from error
    tensor = torch.as_tensor(output, dtype=torch.float32)
    if tuple(tensor.shape) != (1, *IMAGE_SIZE):
        raise ExtractionError(
            f"Normalized {path} has shape {tuple(tensor.shape)}; expected "
            f"{(1, *IMAGE_SIZE)}."
        )
    if not bool(torch.isfinite(tensor).all()):
        raise ExtractionError(f"Normalized volume contains non-finite values: {path}")
    return tensor


def extract_tokens(
    model: Any,
    batch: Any,
    layer: str | int,
    runtime: Mapping[str, Any],
) -> Any:
    torch = runtime["torch"]
    with torch.inference_mode():
        output = model(batch)
    if not isinstance(output, (tuple, list)) or len(output) != 2:
        raise ExtractionError(
            "MONAI ViT did not return (final_tokens, hidden_states); check the "
            "pinned MONAI version."
        )
    final_tokens, hidden_states = output
    if layer == "final":
        tokens = final_tokens
    else:
        if len(hidden_states) != NUM_LAYERS:
            raise ExtractionError(
                f"MONAI ViT returned {len(hidden_states)} hidden layers, expected 12."
            )
        tokens = hidden_states[layer - 1]
    expected_shape = (batch.shape[0], EXPECTED_TOKENS, EXPECTED_FEATURES)
    if tuple(tokens.shape) != expected_shape:
        raise ExtractionError(
            f"Encoder returned token shape {tuple(tokens.shape)}, expected "
            f"{expected_shape}. This script never drops or invents a CLS token."
        )
    if not bool(torch.isfinite(tokens).all()):
        raise ExtractionError("Encoder returned non-finite patch tokens.")
    return tokens.detach().to(device="cpu", dtype=torch.float32).numpy()


def chunked(values: Sequence[Any], size: int) -> Iterable[Sequence[Any]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def metadata_row(
    visit: Visit,
    *,
    latent_sha256: str,
    status: str,
    checkpoint_sha256: str,
    layer: str | int,
) -> dict[str, Any]:
    if visit.nifti_path is None:
        raise AssertionError("NIfTI path was not attached.")
    row: dict[str, Any] = dict(visit.original)
    relative_latent = latent_relative_path(visit)
    row.update(
        {
            "subject_id": visit.subject_id,
            "visit_id": str(visit.image_id),
            "time_month": format(visit.model_time_month, ".17g"),
            "viscode": visit.visit_code,
            "time_year": format(visit.time_year, ".17g"),
            "nifti_path": str(visit.nifti_path),
            "nifti_size_bytes": visit.nifti_size_bytes,
            "nifti_mtime_ns": visit.nifti_mtime_ns,
            "latent_path": str(relative_latent),
            "latent_sha256": latent_sha256,
            "extraction_status": status,
            "encoder_name": MODEL_NAME,
            "encoder_architecture": MODEL_ARCHITECTURE,
            "checkpoint_sha256": checkpoint_sha256,
            "representation_type": REPRESENTATION_TYPE,
            "encoder_layer": layer_label(layer),
            "patch_grid": "x".join(str(value) for value in PATCH_GRID),
            "latent_rows": EXPECTED_TOKENS,
            "latent_cols": EXPECTED_FEATURES,
            "latent_dtype": "float32",
            "preprocessing": PREPROCESSING_DESCRIPTION,
        }
    )
    return row


def read_existing_manifest(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        raise ExtractionError(
            f"Existing extraction manifest is unreadable: {path}: {error}"
        ) from error
    if not isinstance(value, dict):
        raise ExtractionError(f"Existing extraction manifest is not an object: {path}")
    return value


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract deterministic 216x768 BrainIAC patch-token matrices from "
            "a longitudinal ADNI manifest."
        )
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--nifti-dir", type=Path, required=True)
    parser.add_argument(
        "--checkpoint",
        type=Path,
        help="Trusted BrainIAC .ckpt/.pt/.pth or backbone.safetensors file.",
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--max-months", type=float, default=48.0)
    parser.add_argument("--min-visits", type=int, default=3)
    parser.add_argument(
        "--max-subjects",
        type=parse_max_subjects,
        default=None,
        metavar="N|all",
        help="Optional reproducible subject cap (default: all eligible subjects).",
    )
    parser.add_argument("--seed", type=parse_seed, default=20260819)
    parser.add_argument(
        "--layer",
        type=parse_layer,
        default="final",
        metavar="final|1..12",
        help=(
            "Use normalized final tokens (default), or the output of a numbered "
            "transformer block before the global final LayerNorm."
        ),
    )
    parser.add_argument(
        "--device",
        default="auto",
        help="auto, cpu, cuda, or cuda:<index> (default: auto).",
    )
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search the NIfTI directory recursively.",
    )
    parser.add_argument(
        "--force-recompute",
        action="store_true",
        help=(
            "Allow a configuration mismatch and atomically replace selected "
            "latents. Without this flag, mismatched output directories fail."
        ),
    )
    parser.add_argument(
        "--plan-only",
        action="store_true",
        help="Validate cohort selection and NIfTI matching without loading a model.",
    )
    parser.add_argument("--progress-every", type=int, default=10)
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING"), default="INFO"
    )
    return parser


def run(args: argparse.Namespace) -> int:
    manifest_path = args.manifest.expanduser().resolve()
    nifti_dir = args.nifti_dir.expanduser().resolve()
    cohort = read_cohort(
        manifest_path,
        max_months=args.max_months,
        min_visits=args.min_visits,
        max_subjects=args.max_subjects,
        seed=args.seed,
    )
    nifti_index = index_nifti_files(nifti_dir, recursive=args.recursive)
    cohort = attach_nifti_paths(cohort, nifti_index)

    plan = {
        "manifest_rows": cohort.manifest_rows,
        "retained_visits_within_window": cohort.retained_within_window,
        "subjects_before_min_visits": cohort.subjects_before_min_visits,
        "eligible_subjects": cohort.eligible_subjects,
        "selected_subjects": len(cohort.selected_subjects),
        "selected_visits": len(cohort.visits),
        "first_time_year": min(visit.time_year for visit in cohort.visits),
        "last_time_year": max(visit.time_year for visit in cohort.visits),
        "nifti_files_indexed": len(nifti_index),
    }
    if args.plan_only:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return 0

    if args.checkpoint is None:
        raise ExtractionError("--checkpoint is required unless --plan-only is used.")
    if args.output_dir is None:
        raise ExtractionError("--output-dir is required unless --plan-only is used.")
    if args.batch_size < 1:
        raise ExtractionError("--batch-size must be at least one.")
    if args.progress_every < 1:
        raise ExtractionError("--progress-every must be at least one.")

    checkpoint_path = args.checkpoint.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not checkpoint_path.is_file():
        raise ExtractionError(f"Checkpoint does not exist: {checkpoint_path}")

    runtime = import_runtime()
    configure_determinism(runtime, args.seed)
    torch = runtime["torch"]
    np = runtime["np"]
    device = resolve_device(args.device, torch)
    versions = runtime_versions(runtime)
    if versions["monai"] != "1.3.2":
        raise ExtractionError(
            "Exact BrainIAC compatibility requires MONAI 1.3.2; found "
            f"{versions['monai']}."
        )

    checkpoint_sha256 = sha256_file(checkpoint_path)
    script_path = Path(__file__).resolve()
    script_sha256 = sha256_file(script_path)
    repo_root = script_path.parents[2]
    provenance = git_provenance(repo_root)
    configuration = {
        "schema_version": SCHEMA_VERSION,
        "script_version": SCRIPT_VERSION,
        "script_sha256": script_sha256,
        "manifest_sha256": sha256_file(manifest_path),
        "cohort_sha256": cohort_fingerprint(cohort),
        "checkpoint_sha256": checkpoint_sha256,
        "model": {
            "name": MODEL_NAME,
            "architecture": MODEL_ARCHITECTURE,
            "source_image_size": list(SOURCE_IMAGE_SIZE),
            "foreground_roi_start": list(FOREGROUND_ROI_START),
            "foreground_roi_end_exclusive": list(FOREGROUND_ROI_END),
            "cropped_image_size": list(CROPPED_IMAGE_SIZE),
            "image_size": list(IMAGE_SIZE),
            "patch_size": list(PATCH_SIZE),
            "patch_grid": list(PATCH_GRID),
            "num_patches": EXPECTED_TOKENS,
            "hidden_size": EXPECTED_FEATURES,
            "num_layers": NUM_LAYERS,
            "num_heads": NUM_HEADS,
            "mlp_dim": MLP_DIM,
            "classification": False,
            "save_attn": True,
        },
        "representation": {
            "type": REPRESENTATION_TYPE,
            "layer": layer_label(args.layer),
            "shape": [EXPECTED_TOKENS, EXPECTED_FEATURES],
            "dtype": "float32",
            "cls_token_removed": False,
            "cls_token_present": False,
        },
        "preprocessing": PREPROCESSING_DESCRIPTION,
        "selection": {
            "max_months": args.max_months,
            "min_visits": args.min_visits,
            "max_subjects": args.max_subjects,
            "seed": args.seed,
            "baseline_visit_code": "bl",
            "baseline_forced_to_month": 0.0,
        },
        "execution": {
            "requested_device": args.device,
            "resolved_device": str(device),
            "cuda_device_name": (
                torch.cuda.get_device_name(device)
                if device.type == "cuda"
                else None
            ),
            "batch_size": args.batch_size,
            "inference_mode": True,
            "model_eval": True,
            "frozen_parameters": True,
            "deterministic_algorithms": True,
            "cublas_workspace_config": os.environ.get("CUBLAS_WORKSPACE_CONFIG"),
        },
        "runtime_versions": versions,
        "git_commit": provenance["commit"],
    }
    configuration_sha256 = sha256_json(configuration)

    extraction_manifest_path = output_dir / "extraction_manifest.json"
    metadata_path = output_dir / "metadata.csv"
    existing_manifest = read_existing_manifest(extraction_manifest_path)
    configuration_matches = (
        existing_manifest is not None
        and existing_manifest.get("configuration_sha256") == configuration_sha256
    )
    if existing_manifest is not None and not configuration_matches and not args.force_recompute:
        raise ExtractionError(
            "Output directory contains an extraction_manifest.json for a "
            "different configuration. Use a distinct --output-dir or explicitly "
            "pass --force-recompute."
        )
    cache_allowed = configuration_matches and not args.force_recompute

    output_dir.mkdir(parents=True, exist_ok=True)
    started_at = (
        existing_manifest.get("started_at")
        if configuration_matches and existing_manifest is not None
        else utc_now()
    )
    extraction_manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "status": "in_progress",
        "started_at": started_at,
        "updated_at": utc_now(),
        "configuration_sha256": configuration_sha256,
        "configuration": configuration,
        "paths": {
            "manifest": str(manifest_path),
            "nifti_dir": str(nifti_dir),
            "checkpoint": str(checkpoint_path),
            "output_dir": str(output_dir),
            "metadata": str(metadata_path),
        },
        "cohort": {**plan, "selected_subject_ids": list(cohort.selected_subjects)},
        "git": provenance,
        "device": {
            "requested": args.device,
            "resolved": str(device),
            "cuda_available": bool(torch.cuda.is_available()),
            "cuda_device_name": (
                torch.cuda.get_device_name(device)
                if device.type == "cuda"
                else None
            ),
        },
    }
    atomic_write_json(extraction_manifest_path, extraction_manifest)

    try:
        LOGGER.info("Loading %s checkpoint %s", MODEL_NAME, checkpoint_path)
        model = create_model(runtime)
        raw_state = load_checkpoint_mapping(checkpoint_path, runtime)
        checkpoint_load = load_model_weights(model, raw_state)
        model.requires_grad_(False)
        model.eval()
        model.to(device)
        if model.training or any(parameter.requires_grad for parameter in model.parameters()):
            raise ExtractionError("Failed to freeze the BrainIAC encoder.")
        extraction_manifest["checkpoint_load"] = checkpoint_load
        extraction_manifest["updated_at"] = utc_now()
        atomic_write_json(extraction_manifest_path, extraction_manifest)

        load_resize, normalize = build_preprocessors(runtime)
        expected_shape = (EXPECTED_TOKENS, EXPECTED_FEATURES)
        cached: dict[int, tuple[str, str]] = {}
        pending: list[Visit] = []
        for visit in cohort.visits:
            latent_path = output_dir / latent_relative_path(visit)
            if cache_allowed:
                valid, reason = validate_cached_latent(latent_path, np, expected_shape)
                if valid:
                    cached[visit.image_id] = (sha256_file(latent_path), "cached")
                    continue
                if latent_path.exists():
                    LOGGER.warning("Recomputing invalid cache %s: %s", latent_path, reason)
            pending.append(visit)

        previous_metadata_path = (
            quarantine_existing_metadata(metadata_path) if pending else None
        )
        extraction_manifest["paths"]["previous_metadata"] = (
            str(previous_metadata_path) if previous_metadata_path is not None else None
        )
        extraction_manifest["updated_at"] = utc_now()
        atomic_write_json(extraction_manifest_path, extraction_manifest)

        completed: dict[int, tuple[str, str]] = dict(cached)
        processed_count = 0
        for batch_visits in chunked(pending, args.batch_size):
            tensors = []
            for visit in batch_visits:
                if visit.nifti_path is None:
                    raise AssertionError("NIfTI path was not attached.")
                tensors.append(
                    preprocess_volume(
                        visit.nifti_path, load_resize, normalize, runtime
                    )
                )
            batch = torch.stack(tensors, dim=0).to(device=device, dtype=torch.float32)
            token_batch = extract_tokens(model, batch, args.layer, runtime)
            if tuple(token_batch.shape) != (
                len(batch_visits), EXPECTED_TOKENS, EXPECTED_FEATURES
            ):
                raise ExtractionError(
                    f"Unexpected NumPy token batch shape: {token_batch.shape}"
                )
            for visit, tokens in zip(batch_visits, token_batch, strict=True):
                tokens = np.ascontiguousarray(tokens, dtype=np.float32)
                if tuple(tokens.shape) != expected_shape or not bool(np.isfinite(tokens).all()):
                    raise ExtractionError(
                        f"Invalid latent matrix for I{visit.image_id}: "
                        f"shape={tokens.shape}, dtype={tokens.dtype}."
                    )
                latent_path = output_dir / latent_relative_path(visit)
                atomic_save_npy(latent_path, tokens, np)
                valid, reason = validate_cached_latent(latent_path, np, expected_shape)
                if not valid:
                    raise ExtractionError(
                        f"Atomic output validation failed for {latent_path}: {reason}"
                    )
                completed[visit.image_id] = (sha256_file(latent_path), "extracted")
                processed_count += 1
                if processed_count % args.progress_every == 0:
                    LOGGER.info(
                        "Extracted %d / %d pending visits (%d cached)",
                        processed_count,
                        len(pending),
                        len(cached),
                    )
            del batch, token_batch, tensors

        if len(completed) != len(cohort.visits):
            raise ExtractionError(
                f"Internal completion mismatch: {len(completed)} outputs for "
                f"{len(cohort.visits)} visits."
            )

        metadata_rows = [
            metadata_row(
                visit,
                latent_sha256=completed[visit.image_id][0],
                status=completed[visit.image_id][1],
                checkpoint_sha256=checkpoint_sha256,
                layer=args.layer,
            )
            for visit in cohort.visits
        ]
        atomic_write_metadata(
            metadata_path,
            [*cohort.fieldnames, *DERIVED_METADATA_COLUMNS],
            metadata_rows,
        )

        extraction_manifest.update(
            {
                "status": "complete",
                "updated_at": utc_now(),
                "completed_at": utc_now(),
                "outputs": {
                    "visits": len(cohort.visits),
                    "extracted": len(pending),
                    "cached": len(cached),
                    "metadata_sha256": sha256_file(metadata_path),
                },
            }
        )
        atomic_write_json(extraction_manifest_path, extraction_manifest)
        LOGGER.info(
            "Complete: %d visits (%d extracted, %d cached) in %s",
            len(cohort.visits),
            len(pending),
            len(cached),
            output_dir,
        )
        return 0
    except Exception as error:
        extraction_manifest.update(
            {
                "status": "failed",
                "updated_at": utc_now(),
                "failure": {
                    "type": type(error).__name__,
                    "message": str(error),
                },
            }
        )
        atomic_write_json(extraction_manifest_path, extraction_manifest)
        raise


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    try:
        return run(args)
    except ExtractionError as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
