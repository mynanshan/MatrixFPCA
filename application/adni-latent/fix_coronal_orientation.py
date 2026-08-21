#!/usr/bin/env python3
"""Create upright copies of the ADNI coronal PNG slices.

The source exporter placed increasing RAS ``z`` coordinates from the top to
the bottom of each PNG.  Correcting the display orientation therefore needs
one operation only: a top-to-bottom (vertical) flip.  This utility always
writes to a separate directory and records checksums for reproducibility.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import sys
import tempfile
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont, ImageOps


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = PROJECT_ROOT / "data" / "ADNI" / "affine_coronal_slice110"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "data" / "ADNI" / "affine_coronal_slice110_upright"
MANIFEST_NAME = "orientation_manifest.csv"
OPERATION = "vertical_flip_top_bottom"
GRAYSCALE_MODES = {"1", "L", "I", "I;16"}


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


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def display_path(path: Path) -> str:
    """Prefer portable repository-relative paths in the manifest."""

    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(path)


def png_sources(input_dir: Path, max_images: int | None) -> list[Path]:
    sources = sorted(
        (
            path
            for path in input_dir.iterdir()
            if path.is_file() and path.suffix.lower() == ".png"
        ),
        key=lambda path: path.name,
    )
    if max_images is not None:
        sources = sources[:max_images]
    return sources


def save_png_atomically(image: Image.Image, destination: Path) -> None:
    """Write a deterministic PNG, replacing the destination only on success."""

    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
        image.save(temporary_path, format="PNG", optimize=False, compress_level=9)
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def flip_one(source: Path, destination: Path, overwrite: bool) -> dict[str, str]:
    source_resolved = source.resolve(strict=True)
    destination_resolved = destination.resolve(strict=False)
    if source_resolved == destination_resolved:
        raise ValueError(f"Refusing to overwrite source image: {source}")

    source_hash = sha256_file(source)
    destination_existed = destination.exists()
    if destination_existed and not overwrite:
        return {
            "source_path": display_path(source_resolved),
            "output_path": display_path(destination_resolved),
            "source_sha256": source_hash,
            "output_sha256": sha256_file(destination),
            "operation": OPERATION,
            "status": "skipped_existing",
        }

    with Image.open(source) as image:
        image.load()
        if image.format != "PNG":
            raise ValueError(f"Expected a PNG image, got {image.format!r}: {source}")
        if image.mode not in GRAYSCALE_MODES:
            raise ValueError(
                f"Expected a grayscale PNG, got mode {image.mode!r}: {source}"
            )
        upright = image.transpose(Image.Transpose.FLIP_TOP_BOTTOM)
        if upright.mode != image.mode:
            raise RuntimeError(
                f"Vertical flip changed image mode from {image.mode!r} "
                f"to {upright.mode!r}: {source}"
            )
        save_png_atomically(upright, destination)

    return {
        "source_path": display_path(source_resolved),
        "output_path": display_path(destination_resolved),
        "source_sha256": source_hash,
        "output_sha256": sha256_file(destination),
        "operation": OPERATION,
        "status": "overwritten" if destination_existed else "written",
    }


def write_manifest(rows: Iterable[dict[str, str]], manifest_path: Path) -> None:
    fieldnames = [
        "source_path",
        "output_path",
        "source_sha256",
        "output_sha256",
        "operation",
        "status",
    ]
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            dir=manifest_path.parent,
            prefix=f".{manifest_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            writer = csv.DictWriter(temporary, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_path, manifest_path)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def shortened_label(draw: ImageDraw.ImageDraw, label: str, max_width: int) -> str:
    font = ImageFont.load_default()
    if draw.textbbox((0, 0), label, font=font)[2] <= max_width:
        return label
    suffix = "..."
    while label and draw.textbbox((0, 0), label + suffix, font=font)[2] > max_width:
        label = label[:-1]
    return label + suffix


def write_contact_sheet(
    image_paths: list[Path],
    destination: Path,
    columns: int,
    max_images: int,
    overwrite: bool,
) -> str:
    destination_existed = destination.exists()
    if destination_existed and not overwrite:
        return "skipped_existing"

    selected = image_paths[:max_images]
    if not selected:
        raise ValueError("Cannot create a contact sheet without output images")

    label_height = 18
    tile_padding = 4
    with Image.open(selected[0]) as first:
        first.load()
        tile_width, tile_height = first.size

    rows = math.ceil(len(selected) / columns)
    sheet = Image.new(
        "L",
        (
            columns * (tile_width + 2 * tile_padding),
            rows * (tile_height + label_height + 2 * tile_padding),
        ),
        color=0,
    )
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.load_default()

    for index, path in enumerate(selected):
        row, column = divmod(index, columns)
        left = column * (tile_width + 2 * tile_padding) + tile_padding
        top = row * (tile_height + label_height + 2 * tile_padding) + tile_padding
        with Image.open(path) as image:
            image.load()
            if image.mode not in GRAYSCALE_MODES:
                raise ValueError(f"Expected a grayscale contact-sheet input: {path}")
            tile = ImageOps.contain(image, (tile_width, tile_height))
            if tile.mode != "L":
                tile = tile.convert("L")
            sheet.paste(tile, (left, top))
        label = shortened_label(draw, path.stem, tile_width)
        draw.text((left, top + tile_height + 2), label, fill=255, font=font)

    save_png_atomically(sheet, destination)
    return "overwritten" if destination_existed else "written"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Vertically flip the ADNI coronal PNG slices into a separate output "
            "directory and record SHA-256 provenance. Relative paths are resolved "
            "from the repository root."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"source PNG directory (default: {DEFAULT_INPUT_DIR.relative_to(PROJECT_ROOT)})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"separate output directory (default: {DEFAULT_OUTPUT_DIR.relative_to(PROJECT_ROOT)})",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="replace existing output images and contact sheet; source files are never changed",
    )
    parser.add_argument(
        "--max-images",
        type=positive_int,
        default=None,
        help="process only the first N source filenames (useful for smoke tests)",
    )
    parser.add_argument(
        "--contact-sheet",
        action="store_true",
        help="write contact_sheet.png in the output directory",
    )
    parser.add_argument(
        "--contact-sheet-max-images",
        type=positive_int,
        default=16,
        help="maximum number of output images in the contact sheet (default: 16)",
    )
    parser.add_argument(
        "--contact-sheet-columns",
        type=positive_int,
        default=4,
        help="number of contact-sheet columns (default: 4)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    input_dir = resolve_from_project(args.input_dir)
    output_dir = resolve_from_project(args.output_dir)

    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    if input_dir == output_dir:
        raise ValueError(
            "Input and output directories resolve to the same path; refusing to "
            "risk overwriting source images"
        )

    sources = png_sources(input_dir, args.max_images)
    if not sources:
        raise FileNotFoundError(f"No PNG images found in input directory: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / MANIFEST_NAME
    rows: list[dict[str, str]] = []
    output_images: list[Path] = []
    for source in sources:
        destination = output_dir / source.name
        rows.append(flip_one(source, destination, args.overwrite))
        output_images.append(destination)

    write_manifest(rows, manifest_path)

    counts: dict[str, int] = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    summary = ", ".join(f"{status}={count}" for status, count in sorted(counts.items()))
    print(f"Processed {len(rows)} image(s): {summary}")
    print(f"Manifest: {manifest_path}")

    if args.contact_sheet:
        contact_sheet_path = output_dir / "contact_sheet.png"
        status = write_contact_sheet(
            output_images,
            contact_sheet_path,
            columns=args.contact_sheet_columns,
            max_images=args.contact_sheet_max_images,
            overwrite=args.overwrite,
        )
        print(f"Contact sheet ({status}): {contact_sheet_path}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, NotADirectoryError, PermissionError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2) from error
