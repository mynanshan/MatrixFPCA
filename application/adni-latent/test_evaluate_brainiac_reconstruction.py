#!/usr/bin/env python3
"""Dependency-light tests for the BrainIAC reconstruction audit."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np


SCRIPT = Path(__file__).with_name("evaluate_brainiac_reconstruction.py")
SPEC = importlib.util.spec_from_file_location("evaluate_brainiac_reconstruction", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
audit = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)


class ReconstructionAuditTests(unittest.TestCase):
    def test_subject_split_is_disjoint_and_reproducible(self) -> None:
        rows = [
            {"subject_id": subject, "image_id": str(visit)}
            for subject in (f"S{value}" for value in range(10))
            for visit in range(3)
        ]
        first = audit.make_subject_split(rows, 42, 0.2, 0.25)
        second = audit.make_subject_split(rows, 42, 0.2, 0.25)
        self.assertEqual(first, second)
        train, validation, test = first
        self.assertFalse(train & validation)
        self.assertFalse(train & test)
        self.assertFalse(validation & test)
        self.assertEqual(train | validation | test, {f"S{value}" for value in range(10)})

    def test_kernel_ridge_recovers_linear_targets(self) -> None:
        rng = np.random.default_rng(7)
        features = rng.normal(size=(12, 5)).astype(np.float32)
        coefficient = rng.normal(size=(5, 8)).astype(np.float32)
        images = features @ coefficient + 2.0
        development = np.arange(10)
        test = np.arange(10, 12)
        prediction, _ = audit.fit_predict(features, images, development, test, 1e-8)
        np.testing.assert_allclose(prediction, images[test], atol=2e-4, rtol=2e-4)

    def test_relative_rmse_is_zero_for_exact_prediction(self) -> None:
        truth = np.asarray([0.0, 1.0, 2.0, -1.0], dtype=np.float32)
        metrics = audit.scan_metrics(truth, truth.copy())
        self.assertEqual(metrics["rmse"], 0.0)
        self.assertEqual(metrics["relative_rmse"], 0.0)
        self.assertAlmostEqual(metrics["correlation"], 1.0)


if __name__ == "__main__":
    unittest.main()
