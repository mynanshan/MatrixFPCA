#!/usr/bin/env python3
"""Lightweight tests that do not require Torch, MONAI, or a checkpoint."""

from __future__ import annotations

import csv
import importlib.util
import numpy as np
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("extract_brainiac_latents.py")
SPEC = importlib.util.spec_from_file_location("extract_brainiac_latents", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
extractor = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = extractor
SPEC.loader.exec_module(extractor)


FIELDNAMES = [
    "PTID",
    "image_id",
    "image_time_month",
    "mri_VISCODE2",
    "image_retained_for_analysis",
    "baseline_age",
]


class CohortTests(unittest.TestCase):
    def test_command_line_default_uses_all_eligible_subjects(self) -> None:
        parser = extractor.build_argument_parser()
        options = parser.parse_args(["--manifest", "manifest.csv", "--nifti-dir", "nifti"])
        self.assertIsNone(options.max_subjects)

    def write_manifest(self, directory: Path, rows: list[dict[str, str]]) -> Path:
        path = directory / "manifest.csv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
            writer.writeheader()
            writer.writerows(rows)
        return path

    def make_niftis(self, directory: Path, image_ids: list[int]) -> Path:
        root = directory / "nifti"
        root.mkdir()
        for image_id in image_ids:
            (root / f"002_S_0001_2000-01-01_I{image_id}_Warped.nii.gz").write_bytes(
                b"not-a-real-nifti"
            )
        return root

    def test_selection_forces_baseline_and_filters_subjects(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = [
                {
                    "PTID": "A_S_0001",
                    "image_id": "1",
                    "image_time_month": "-0.25",
                    "mri_VISCODE2": "bl",
                    "image_retained_for_analysis": "TRUE",
                    "baseline_age": "70",
                },
                {
                    "PTID": "A_S_0001",
                    "image_id": "2",
                    "image_time_month": "6",
                    "mri_VISCODE2": "m06",
                    "image_retained_for_analysis": "TRUE",
                    "baseline_age": "70",
                },
                {
                    "PTID": "A_S_0001",
                    "image_id": "3",
                    "image_time_month": "60",
                    "mri_VISCODE2": "m60",
                    "image_retained_for_analysis": "TRUE",
                    "baseline_age": "70",
                },
                {
                    "PTID": "B_S_0002",
                    "image_id": "4",
                    "image_time_month": "0",
                    "mri_VISCODE2": "sc",
                    "image_retained_for_analysis": "TRUE",
                    "baseline_age": "71",
                },
                {
                    "PTID": "B_S_0002",
                    "image_id": "5",
                    "image_time_month": "12",
                    "mri_VISCODE2": "m12",
                    "image_retained_for_analysis": "FALSE",
                    "baseline_age": "71",
                },
            ]
            manifest = self.write_manifest(root, rows)
            cohort = extractor.read_cohort(
                manifest,
                max_months=48,
                min_visits=2,
                max_subjects=None,
                seed=1,
            )
            self.assertEqual(cohort.selected_subjects, ("A_S_0001",))
            self.assertEqual([visit.image_id for visit in cohort.visits], [1, 2])
            self.assertEqual(cohort.visits[0].image_time_month, -0.25)
            self.assertEqual(cohort.visits[0].model_time_month, 0.0)
            self.assertEqual(cohort.visits[0].time_year, 0.0)

    def test_nifti_matching_and_duplicate_detection(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            nifti_root = self.make_niftis(root, [12, 13])
            index = extractor.index_nifti_files(nifti_root, recursive=False)
            self.assertEqual(set(index), {12, 13})
            (nifti_root / "duplicate_I12_copy.nii").write_bytes(b"duplicate")
            with self.assertRaisesRegex(extractor.ExtractionError, "not provide unique"):
                extractor.index_nifti_files(nifti_root, recursive=False)

    def test_subject_cap_is_reproducible(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows: list[dict[str, str]] = []
            image_id = 0
            for subject in ("S1", "S2", "S3", "S4"):
                for month in (0, 12):
                    image_id += 1
                    rows.append(
                        {
                            "PTID": subject,
                            "image_id": str(image_id),
                            "image_time_month": str(month),
                            "mri_VISCODE2": "bl" if month == 0 else "m12",
                            "image_retained_for_analysis": "TRUE",
                            "baseline_age": "70",
                        }
                    )
            manifest = self.write_manifest(root, rows)
            first = extractor.read_cohort(
                manifest,
                max_months=48,
                min_visits=2,
                max_subjects=2,
                seed=99,
            )
            second = extractor.read_cohort(
                manifest,
                max_months=48,
                min_visits=2,
                max_subjects=2,
                seed=99,
            )
            self.assertEqual(first.selected_subjects, second.selected_subjects)
            self.assertEqual(len(first.selected_subjects), 2)

    def test_metadata_contract_for_r_pipeline(self) -> None:
        visit = extractor.Visit(
            original={field: "" for field in FIELDNAMES},
            source_row=2,
            subject_id="A_S_0001",
            image_id=42,
            visit_code="bl",
            image_time_month=-0.25,
            model_time_month=0.0,
            time_year=0.0,
            nifti_path=Path("/tmp/example_I42_Warped.nii.gz"),
            nifti_size_bytes=100,
            nifti_mtime_ns=200,
        )
        row = extractor.metadata_row(
            visit,
            latent_sha256="a" * 64,
            status="extracted",
            checkpoint_sha256="b" * 64,
            layer="final",
        )
        required = {
            "subject_id",
            "visit_id",
            "image_id",
            "time_month",
            "viscode",
            "latent_path",
            "nifti_path",
            "latent_rows",
            "latent_cols",
        }
        self.assertTrue(required.issubset(row))
        self.assertEqual(row["time_month"], "0")
        self.assertEqual(row["viscode"], "bl")
        self.assertEqual(row["latent_rows"], 216)
        self.assertEqual(row["latent_cols"], 768)

    def test_cache_requires_finite_float32_expected_shape(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "latent.npy"
            np.save(path, np.zeros((216, 768), dtype=np.float32))
            self.assertEqual(
                extractor.validate_cached_latent(path, np, (216, 768)),
                (True, "valid"),
            )

            np.save(path, np.zeros((216, 768), dtype=np.float64))
            valid, reason = extractor.validate_cached_latent(path, np, (216, 768))
            self.assertFalse(valid)
            self.assertIn("dtype", reason)

            invalid = np.zeros((216, 768), dtype=np.float32)
            invalid[0, 0] = np.nan
            np.save(path, invalid)
            valid, reason = extractor.validate_cached_latent(path, np, (216, 768))
            self.assertFalse(valid)
            self.assertIn("non-finite", reason)

    def test_near_duplicate_visits_use_r_day_resolution(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = []
            for image_id, month in enumerate((0.0, 1.0, 1.001), start=1):
                rows.append(
                    {
                        "PTID": "A_S_0001",
                        "image_id": str(image_id),
                        "image_time_month": str(month),
                        "mri_VISCODE2": "bl" if image_id == 1 else f"v{image_id}",
                        "image_retained_for_analysis": "TRUE",
                        "baseline_age": "70",
                    }
                )
            manifest = self.write_manifest(root, rows)
            with self.assertRaisesRegex(extractor.ExtractionError, "subject/day"):
                extractor.read_cohort(
                    manifest,
                    max_months=48,
                    min_visits=3,
                    max_subjects=None,
                    seed=1,
                )

    def test_existing_metadata_is_quarantined_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            metadata = Path(temporary) / "metadata.csv"
            metadata.write_text("old metadata\n", encoding="utf-8")
            previous = extractor.quarantine_existing_metadata(metadata)
            self.assertFalse(metadata.exists())
            self.assertEqual(previous, metadata.with_name("metadata.previous.csv"))
            self.assertEqual(previous.read_text(encoding="utf-8"), "old metadata\n")


if __name__ == "__main__":
    unittest.main()
