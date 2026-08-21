#!/usr/bin/env python3
"""Small dependency-free tests for the BrainIAC global/local concept audit."""

from __future__ import annotations

import unittest

import numpy as np

import evaluate_brainiac_concepts as audit


class ConceptAuditTests(unittest.TestCase):
    def test_subject_weights_give_each_subject_equal_mass(self) -> None:
        rows = [{"subject_id": value} for value in ("a", "a", "a", "b")]
        weights = audit.subject_weights(rows, np.arange(4))
        self.assertAlmostEqual(float(weights[:3].sum()), float(weights[3]))

    def test_auc_handles_ties(self) -> None:
        truth = np.asarray([0, 0, 1, 1], dtype=float)
        score = np.asarray([0.0, 0.5, 0.5, 1.0])
        self.assertAlmostEqual(audit.auc(truth, score), 0.875)

    def test_patch_counts_follow_c_order_tokenization(self) -> None:
        mask = np.zeros((96, 96, 96), dtype=bool)
        mask[:16, :16, 16:32] = True
        counts = audit.patch_counts(mask)
        self.assertEqual(int(counts[1]), 16**3)
        self.assertEqual(int(np.count_nonzero(counts)), 1)

    def test_weighted_kernel_ridge_recovers_linear_signal(self) -> None:
        rng = np.random.default_rng(4)
        features = rng.normal(size=(18, 5)).astype(np.float32)
        targets = (features @ np.arange(1, 6, dtype=np.float32))[:, None]
        fit = np.arange(12)
        predict = np.arange(12, 18)
        rows = [{"subject_id": str(index)} for index in range(18)]
        weights = audit.subject_weights(rows, fit)
        kernel, cross = audit.prepare_kernel(features, fit, predict, weights)
        estimate, _ = audit.fit_from_kernel(kernel, cross, targets[fit], weights, 1e-8)
        self.assertLess(float(np.max(np.abs(estimate - targets[predict]))), 1e-4)


if __name__ == "__main__":
    unittest.main()
