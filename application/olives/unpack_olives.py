#!/usr/bin/env python3
"""Extract a deterministic, QC-only image sample from local OLIVES parquet.

The Hugging Face parquet files embed each image in ``Image.bytes`` and retain
only a basename in ``Image.path``.  They do *not* contain a visit/week field.
Consequently this utility deliberately does not infer OCT volumes, B-scan
order, visits, or longitudinal time from row order or filenames.

The extractor reads only selected parquet row groups, decodes every selected
payload with Pillow, writes the original encoded bytes, and records provenance
and available scalar labels in a CSV manifest.  Results are built in a staging
directory and published with one directory rename; an existing output path is
never overwritten.

Example::

    application/olives/.venv/bin/python \
        application/olives/unpack_olives.py --samples 16 --contact-sheet
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import io
import json
import math
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Any, Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = PROJECT_ROOT / "data" / "OLIVES" / "disease_classification"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "output" / "olives_qc"
MANIFEST_NAME = "manifest.csv"
NOTICE_NAME = "QC_ONLY.txt"
CONTACT_SHEET_NAME = "contact_sheet.png"

QC_NOTICE = """OLIVES PARQUET QC SAMPLE -- NOT LONGITUDINAL DATA

These files were extracted only for visual and metadata quality control.
This tool does not infer visits, weeks, OCT volumes, B-scan position, treatment
time, or longitudinal order. In the local Hugging Face parquet export,
Image.path is generally a basename rather than the original full path, no
visit/week column is present, Scan (n/49) is missing for most rows, and row
order is not a safe longitudinal key.

Do not treat sample_index, dataset_row, shard_row, filename order, or contact-
sheet position as time. A scientifically usable longitudinal manifest must be
joined from the original OLIVES full File_Path/Week metadata.

sample_index is one-based. dataset_row, shard_index, and shard_row are
zero-based provenance indices in the sorted list of parquet shard filenames.
"""

CORE_MANIFEST_FIELDS = [
    "sample_index",
    "dataset_row",
    "shard_index",
    "shard",
    "shard_row",
    "original_path",
    "image_sha256",
    "encoded_bytes",
    "output_path",
    "source_format",
    "width",
    "height",
    "mode",
    "longitudinal_grouping_inferred",
]

FORMAT_EXTENSIONS = {
    "BMP": ".bmp",
    "GIF": ".gif",
    "JPEG": ".jpg",
    "PNG": ".png",
    "TIFF": ".tif",
    "WEBP": ".webp",
}


@dataclass(frozen=True)
class ShardInfo:
    index: int
    path: Path
    first_global_row: int
    row_count: int
    metadata_fields: tuple[str, ...]


@dataclass(frozen=True)
class Selection:
    sample_index: int
    dataset_row: int
    shard_index: int
    shard_row: int


@dataclass(frozen=True)
class ExtractedSample:
    selection: Selection
    output_path: Path
    original_path: str
    image_sha256: str
    encoded_bytes: int
    source_format: str
    width: int
    height: int
    mode: str
    metadata: dict[str, Any]


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("value must be at least 1")
    return parsed


def resolve_from_project(path: Path) -> Path:
    """Resolve relative CLI paths from the repository root."""

    if not path.is_absolute():
        path = PROJECT_ROOT / path
    return path.resolve()


def display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def require_dependencies() -> tuple[Any, Any, Any, Any, Any, Any]:
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
        from PIL import Image, ImageDraw, ImageFont, ImageOps
    except ImportError as error:
        raise RuntimeError(
            "unpack_olives.py requires pyarrow and Pillow. Install "
            "application/olives/requirements-qc.txt in a dedicated environment."
        ) from error
    return pa, pq, Image, ImageDraw, ImageFont, ImageOps


def is_scalar_manifest_type(data_type: Any, pa: Any) -> bool:
    checks = (
        pa.types.is_boolean,
        pa.types.is_integer,
        pa.types.is_floating,
        pa.types.is_decimal,
        pa.types.is_string,
        pa.types.is_large_string,
        pa.types.is_date,
        pa.types.is_time,
        pa.types.is_timestamp,
        pa.types.is_duration,
    )
    return any(check(data_type) for check in checks)


def validate_image_field(schema: Any, shard: Path, pa: Any) -> None:
    index = schema.get_field_index("Image")
    if index < 0:
        raise ValueError(f"Missing required Image column: {shard}")
    image_type = schema.field(index).type
    if not pa.types.is_struct(image_type):
        raise ValueError(f"Expected Image to be a struct column in {shard}")
    child_names = {field.name for field in image_type}
    if "bytes" not in child_names:
        raise ValueError(f"Image struct has no embedded bytes field: {shard}")
    bytes_type = image_type.field("bytes").type
    if not (pa.types.is_binary(bytes_type) or pa.types.is_large_binary(bytes_type)):
        raise ValueError(f"Image.bytes is not a binary field: {shard}")


def discover_shards(
    input_dir: Path, pattern: str, pq: Any, pa: Any
) -> list[ShardInfo]:
    paths = sorted(
        (path.resolve() for path in input_dir.glob(pattern) if path.is_file()),
        key=lambda path: path.name,
    )
    if not paths:
        raise FileNotFoundError(
            f"No parquet files matching {pattern!r} in input directory: {input_dir}"
        )

    shards: list[ShardInfo] = []
    first_global_row = 0
    for index, path in enumerate(paths):
        parquet = pq.ParquetFile(path)
        schema = parquet.schema_arrow
        validate_image_field(schema, path, pa)
        row_count = parquet.metadata.num_rows
        if row_count < 1:
            raise ValueError(f"Parquet shard contains no rows: {path}")
        metadata_fields = tuple(
            field.name
            for field in schema
            if field.name != "Image" and is_scalar_manifest_type(field.type, pa)
        )
        shards.append(
            ShardInfo(
                index=index,
                path=path,
                first_global_row=first_global_row,
                row_count=row_count,
                metadata_fields=metadata_fields,
            )
        )
        first_global_row += row_count
    return shards


def ordered_metadata_fields(shards: Iterable[ShardInfo]) -> list[str]:
    fields: list[str] = []
    seen: set[str] = set()
    for shard in shards:
        for field in shard.metadata_fields:
            if field not in seen:
                seen.add(field)
                fields.append(field)
    collisions = set(fields).intersection(CORE_MANIFEST_FIELDS)
    if collisions:
        joined = ", ".join(sorted(collisions))
        raise ValueError(f"Parquet metadata conflicts with manifest fields: {joined}")
    return fields


def evenly_spaced_indices(total_rows: int, requested: int) -> list[int]:
    """Return deterministic midpoint-quantile row indices."""

    count = min(total_rows, requested)
    return [((2 * index + 1) * total_rows) // (2 * count) for index in range(count)]


def choose_samples(shards: list[ShardInfo], requested: int) -> list[Selection]:
    total_rows = sum(shard.row_count for shard in shards)
    starts = [shard.first_global_row for shard in shards]
    selections: list[Selection] = []
    for sample_index, dataset_row in enumerate(
        evenly_spaced_indices(total_rows, requested), start=1
    ):
        shard_index = bisect.bisect_right(starts, dataset_row) - 1
        shard = shards[shard_index]
        selections.append(
            Selection(
                sample_index=sample_index,
                dataset_row=dataset_row,
                shard_index=shard_index,
                shard_row=dataset_row - shard.first_global_row,
            )
        )
    return selections


def row_group_offsets(parquet: Any) -> list[int]:
    offsets: list[int] = []
    current = 0
    for index in range(parquet.metadata.num_row_groups):
        offsets.append(current)
        current += parquet.metadata.row_group(index).num_rows
    return offsets


def group_selections_by_row_group(
    selections: list[Selection], parquet: Any
) -> dict[int, list[tuple[int, Selection]]]:
    offsets = row_group_offsets(parquet)
    grouped: dict[int, list[tuple[int, Selection]]] = {}
    for selection in selections:
        group_index = bisect.bisect_right(offsets, selection.shard_row) - 1
        row_in_group = selection.shard_row - offsets[group_index]
        grouped.setdefault(group_index, []).append((row_in_group, selection))
    for values in grouped.values():
        values.sort(key=lambda value: value[0])
    return grouped


def atomic_write_bytes(destination: Path, content: bytes) -> None:
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            temporary.write(content)
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def canonical_extension(source_format: str) -> str:
    return FORMAT_EXTENSIONS.get(source_format.upper(), ".img")


def normalize_original_path(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def extract_one(
    image_value: Any,
    metadata: dict[str, Any],
    selection: Selection,
    image_dir: Path,
    Image: Any,
) -> ExtractedSample:
    if not isinstance(image_value, dict):
        raise ValueError(
            f"Image row is not a bytes/path struct at dataset row "
            f"{selection.dataset_row}"
        )
    payload = image_value.get("bytes")
    if isinstance(payload, memoryview):
        payload = payload.tobytes()
    elif isinstance(payload, bytearray):
        payload = bytes(payload)
    if not isinstance(payload, bytes) or not payload:
        raise ValueError(
            f"Image.bytes is empty at dataset row {selection.dataset_row}"
        )

    try:
        with Image.open(io.BytesIO(payload)) as image:
            source_format = (image.format or "UNKNOWN").upper()
            width, height = image.size
            mode = image.mode
            image.load()
    except Exception as error:
        raise ValueError(
            f"Pillow could not decode Image.bytes at dataset row "
            f"{selection.dataset_row}: {error}"
        ) from error

    digest = hashlib.sha256(payload).hexdigest()
    filename = (
        f"sample_{selection.sample_index:04d}_"
        f"shard_{selection.shard_index:04d}_"
        f"row_{selection.shard_row:07d}{canonical_extension(source_format)}"
    )
    output_path = image_dir / filename
    atomic_write_bytes(output_path, payload)
    return ExtractedSample(
        selection=selection,
        output_path=output_path,
        original_path=normalize_original_path(image_value.get("path")),
        image_sha256=digest,
        encoded_bytes=len(payload),
        source_format=source_format,
        width=width,
        height=height,
        mode=mode,
        metadata=metadata,
    )


def extract_selected_rows(
    shards: list[ShardInfo],
    selections: list[Selection],
    staging_dir: Path,
    batch_size: int,
    pq: Any,
    Image: Any,
) -> list[ExtractedSample]:
    image_dir = staging_dir / "images"
    image_dir.mkdir()
    by_shard: dict[int, list[Selection]] = {}
    for selection in selections:
        by_shard.setdefault(selection.shard_index, []).append(selection)

    extracted: dict[int, ExtractedSample] = {}
    for shard_index in sorted(by_shard):
        shard = shards[shard_index]
        parquet = pq.ParquetFile(shard.path)
        grouped = group_selections_by_row_group(by_shard[shard_index], parquet)
        columns = ["Image", *shard.metadata_fields]

        for group_index in sorted(grouped):
            targets = grouped[group_index]
            target_index = 0
            batch_start = 0
            for batch in parquet.iter_batches(
                batch_size=batch_size,
                row_groups=[group_index],
                columns=columns,
                use_threads=False,
            ):
                batch_end = batch_start + batch.num_rows
                while (
                    target_index < len(targets)
                    and targets[target_index][0] < batch_end
                ):
                    row_in_group, selection = targets[target_index]
                    row_in_batch = row_in_group - batch_start
                    image_column = batch.schema.get_field_index("Image")
                    image_value = batch.column(image_column)[row_in_batch].as_py()
                    metadata: dict[str, Any] = {}
                    for field in shard.metadata_fields:
                        column_index = batch.schema.get_field_index(field)
                        metadata[field] = batch.column(column_index)[
                            row_in_batch
                        ].as_py()
                    extracted[selection.sample_index] = extract_one(
                        image_value=image_value,
                        metadata=metadata,
                        selection=selection,
                        image_dir=image_dir,
                        Image=Image,
                    )
                    target_index += 1
                batch_start = batch_end
                if target_index == len(targets):
                    break
            if target_index != len(targets):
                raise RuntimeError(
                    f"Failed to read all selected rows from row group {group_index} "
                    f"of {shard.path}"
                )

    if len(extracted) != len(selections):
        raise RuntimeError(
            f"Expected {len(selections)} extracted rows, got {len(extracted)}"
        )
    return [extracted[index] for index in sorted(extracted)]


def csv_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, (date, datetime)):
        return value.isoformat()
    if isinstance(value, bytes):
        return value.hex()
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def manifest_rows(
    samples: list[ExtractedSample],
    shards: list[ShardInfo],
    staging_dir: Path,
    metadata_fields: list[str],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for sample in samples:
        selection = sample.selection
        shard = shards[selection.shard_index]
        row: dict[str, Any] = {
            "sample_index": selection.sample_index,
            "dataset_row": selection.dataset_row,
            "shard_index": selection.shard_index,
            "shard": display_path(shard.path),
            "shard_row": selection.shard_row,
            "original_path": sample.original_path,
            "image_sha256": sample.image_sha256,
            "encoded_bytes": sample.encoded_bytes,
            "output_path": sample.output_path.relative_to(staging_dir).as_posix(),
            "source_format": sample.source_format,
            "width": sample.width,
            "height": sample.height,
            "mode": sample.mode,
            "longitudinal_grouping_inferred": "false",
        }
        row.update(
            {
                field: csv_value(sample.metadata.get(field))
                for field in metadata_fields
            }
        )
        rows.append(row)
    return rows


def atomic_write_csv(
    destination: Path, rows: list[dict[str, Any]], fieldnames: list[str]
) -> None:
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            writer = csv.DictWriter(
                temporary, fieldnames=fieldnames, lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def shortened_label(draw: Any, label: str, font: Any, max_width: int) -> str:
    if draw.textbbox((0, 0), label, font=font)[2] <= max_width:
        return label
    suffix = "..."
    while label and draw.textbbox((0, 0), label + suffix, font=font)[2] > max_width:
        label = label[:-1]
    return label + suffix


def write_contact_sheet(
    samples: list[ExtractedSample],
    destination: Path,
    columns: int,
    tile_size: int,
    Image: Any,
    ImageDraw: Any,
    ImageFont: Any,
    ImageOps: Any,
) -> None:
    if not samples:
        raise ValueError("Cannot create a contact sheet without samples")
    rows = math.ceil(len(samples) / columns)
    padding = 5
    label_height = 20
    cell_width = tile_size + 2 * padding
    cell_height = tile_size + label_height + 2 * padding
    sheet = Image.new("L", (columns * cell_width, rows * cell_height), color=0)
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.load_default()

    for index, sample in enumerate(samples):
        row, column = divmod(index, columns)
        left = column * cell_width + padding
        top = row * cell_height + padding
        with Image.open(sample.output_path) as image:
            image.load()
            tile = ImageOps.contain(image.convert("L"), (tile_size, tile_size))
        tile_left = left + (tile_size - tile.width) // 2
        tile_top = top + (tile_size - tile.height) // 2
        sheet.paste(tile, (tile_left, tile_top))
        label = (
            f"{sample.selection.sample_index}: "
            f"{sample.original_path or '<no path>'}"
        )
        label = shortened_label(draw, label, font, tile_size)
        draw.text((left, top + tile_size + 2), label, fill=255, font=font)

    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
        sheet.save(temporary_path, format="PNG", optimize=False, compress_level=9)
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def validate_paths(input_dir: Path, output_dir: Path) -> None:
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input directory does not exist: {input_dir}")
    if output_dir.exists() or output_dir.is_symlink():
        raise FileExistsError(
            f"Output path already exists; choose a new directory so prior QC "
            f"results cannot be overwritten: {output_dir}"
        )
    if output_dir == input_dir or input_dir in output_dir.parents:
        raise ValueError(
            "Output directory must not be the parquet input directory or one of "
            "its descendants"
        )
    if output_dir in input_dir.parents:
        raise ValueError(
            "Output directory must not contain the parquet input directory"
        )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract deterministic OLIVES parquet images for QC only. No visit, "
            "volume, scan-order, or longitudinal grouping is inferred. Relative "
            "paths are resolved from the repository root."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=(
            "directory containing parquet shards (default: "
            f"{DEFAULT_INPUT_DIR.relative_to(PROJECT_ROOT)})"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "new output directory; it must not already exist (default: "
            f"{DEFAULT_OUTPUT_DIR.relative_to(PROJECT_ROOT)})"
        ),
    )
    parser.add_argument(
        "--pattern",
        default="*.parquet",
        help="input filename glob within --input-dir (default: *.parquet)",
    )
    parser.add_argument(
        "--samples",
        "--max-images",
        dest="samples",
        type=positive_int,
        default=16,
        help=(
            "number of midpoint-quantile rows sampled across sorted shards "
            "(default: 16)"
        ),
    )
    parser.add_argument(
        "--batch-size",
        type=positive_int,
        default=128,
        help=(
            "maximum Arrow record-batch size while reading selected row groups "
            "(default: 128)"
        ),
    )
    parser.add_argument(
        "--contact-sheet",
        action="store_true",
        help=f"also write a grayscale {CONTACT_SHEET_NAME}",
    )
    parser.add_argument(
        "--contact-sheet-columns",
        type=positive_int,
        default=4,
        help="number of contact-sheet columns (default: 4)",
    )
    parser.add_argument(
        "--contact-sheet-tile-size",
        type=positive_int,
        default=256,
        help="square contact-sheet tile size in pixels (default: 256)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    input_dir = resolve_from_project(args.input_dir)
    output_dir = resolve_from_project(args.output_dir)
    validate_paths(input_dir, output_dir)

    pa, pq, Image, ImageDraw, ImageFont, ImageOps = require_dependencies()
    shards = discover_shards(input_dir, args.pattern, pq, pa)
    metadata_fields = ordered_metadata_fields(shards)
    selections = choose_samples(shards, args.samples)
    selected_shard_count = len({selection.shard_index for selection in selections})

    print(
        "WARNING: QC only; no longitudinal/visit grouping is inferred.",
        file=sys.stderr,
    )
    if len(selections) < args.samples:
        print(
            f"Requested {args.samples} samples but the dataset has only "
            f"{len(selections)} rows; extracting every row.",
            file=sys.stderr,
        )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.staging.", dir=output_dir.parent)
    )
    published = False
    try:
        samples = extract_selected_rows(
            shards=shards,
            selections=selections,
            staging_dir=staging_dir,
            batch_size=args.batch_size,
            pq=pq,
            Image=Image,
        )
        rows = manifest_rows(samples, shards, staging_dir, metadata_fields)
        atomic_write_csv(
            staging_dir / MANIFEST_NAME,
            rows,
            [*CORE_MANIFEST_FIELDS, *metadata_fields],
        )
        atomic_write_bytes(staging_dir / NOTICE_NAME, QC_NOTICE.encode("utf-8"))
        if args.contact_sheet:
            write_contact_sheet(
                samples=samples,
                destination=staging_dir / CONTACT_SHEET_NAME,
                columns=args.contact_sheet_columns,
                tile_size=args.contact_sheet_tile_size,
                Image=Image,
                ImageDraw=ImageDraw,
                ImageFont=ImageFont,
                ImageOps=ImageOps,
            )

        if output_dir.exists() or output_dir.is_symlink():
            raise FileExistsError(
                f"Output path appeared during extraction; refusing to overwrite it: "
                f"{output_dir}"
            )
        staging_dir.rename(output_dir)
        published = True
    finally:
        if not published and staging_dir.exists():
            shutil.rmtree(staging_dir)

    total_rows = sum(shard.row_count for shard in shards)
    print(
        f"Extracted {len(selections)} deterministic QC sample(s) from "
        f"{selected_shard_count}/{len(shards)} shard(s) and {total_rows} total rows."
    )
    print(f"Output: {output_dir}")
    print(f"Manifest: {output_dir / MANIFEST_NAME}")
    if args.contact_sheet:
        print(f"Contact sheet: {output_dir / CONTACT_SHEET_NAME}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        FileExistsError,
        FileNotFoundError,
        NotADirectoryError,
        PermissionError,
        RuntimeError,
        ValueError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2) from error
