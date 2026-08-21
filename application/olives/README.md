# OLIVES latent-matrix experiment

This directory is intentionally separate from `application/adni-latent`. It
contains OLIVES retinal OCT inspection and the preparation path for a future
frozen RETFound latent-matrix experiment.

## Current status

The local Hugging Face parquet export is sufficient for image and label quality
control, but not for defensible longitudinal modeling:

- `disease_classification` contains 78,822 rows across 32 parquet shards; the
  same image corpus is repeated in the biomarker-training export with a
  different label schema.
- Images are embedded in the parquet `Image.bytes` field. `Image.path` retains
  only a basename, not the original directory hierarchy.
- There are 96 eyes from 87 patients. Eye, rather than patient, is the natural
  longitudinal subject; patient can be retained as a clustering identifier.
- No week/date, visit, treatment time, or volume identifier survives in the
  local schema. `Scan (n/49)` is missing for 77.81% of disease rows, and row or
  lexicographic filename order cannot safely reconstruct visits.

The official OLIVES processing code derives `Week` from the full source path
and constructs volumes by `Eye_ID + Week`. That information was stripped from
the downloaded parquet, so this application deliberately does not invent a
time index. See the [official OLIVES repository](https://github.com/olivesgatech/OLIVES_Dataset).

## Visual QC / unpacking

`unpack_olives.py` streams selected parquet row groups, writes the exact
original encoded image bytes, validates them with Pillow, records hashes and
available scalar labels in `manifest.csv`, and can produce a grayscale contact
sheet. It labels every output as QC-only and infers no visit, volume, B-scan
position, or longitudinal order.

Create an isolated environment:

```bash
python -m venv application/olives/.venv
application/olives/.venv/bin/pip install \
  -r application/olives/requirements-qc.txt
```

Extract a fresh deterministic sample:

```bash
application/olives/.venv/bin/python \
  application/olives/unpack_olives.py \
  --samples=16 \
  --contact-sheet \
  --output-dir=application/olives/output/olives_qc_fresh
```

A validated 16-image sample, provenance manifest, QC warning, and contact sheet
already exist in `application/olives/output/olives_qc`. Generated output is
ignored by Git. The source parquet under `data/OLIVES` is never modified.

Do not globally deduplicate images by hash: distinct documented visits can
contain identical scans. Duplicate rows can be collapsed only after restoring
the original full path/week and showing that they refer to the same source
record and visit.

## What is needed for longitudinal OLIVES + RETFound

Please provide either:

1. the original Zenodo OLIVES directory plus `Clinical_Data_Images.xlsx`; or
2. a manifest with full file path, `Eye_ID`, week/date, visit/volume ID, scan
   position, and the desired clinical/treatment variables.

With that manifest, the next implementation can recover `(eye, week, scan)`,
validate 49-scan volumes, select or aggregate consistent B-scans, and encode
each visit as a fixed patch-token matrix.

The [official RETFound code](https://github.com/rmaphoh/RETFound) is public,
but the [official OCT checkpoint](https://huggingface.co/YukunZhou/RETFound_mae_natureOCT)
requires a Hugging Face account and acceptance of its access terms. Supply any
token through the environment and never commit it. RETFound is a ViT-Large
model, so a CUDA GPU is strongly recommended for full extraction.
