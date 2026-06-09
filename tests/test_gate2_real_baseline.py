import unittest

import numpy as np

from pathscope_st.parity import build_parity_table
from pathscope_st.same_split_baseline import patch_grid_features, point_prediction_row


class Gate2RealBaselineTests(unittest.TestCase):
    def test_patch_grid_features_include_rgb_and_grid(self):
        images = np.zeros((3, 8, 8, 3), dtype=np.uint8)
        features = patch_grid_features(images, grid_size=4)
        self.assertEqual(features.shape, (3, 22))

    def test_point_prediction_rows_mark_uncalibrated_same_split(self):
        truth = np.array([[0.0, 1.0], [1.0, 2.0], [2.0, 3.0]], dtype=np.float32)
        pred = truth.copy()
        row = point_prediction_row("patch_feature_ridge", truth, pred, "unit", "unit")
        self.assertEqual(row["provenance"], "RAN")
        self.assertTrue(row["same_split"])
        self.assertFalse(row["calibrated_interval"])
        self.assertIsNone(row["interval_coverage"])
        self.assertAlmostEqual(row["heldout_gene_pearson_mean"], 1.0)

    def test_gate3_complete_parity_leaves_only_license_review_missing(self):
        table = build_parity_table(
            {
                "heldout_gene_pearson_mean": 0.269,
                "heldout_spot_pearson_mean": 0.874,
                "interval_coverage": 0.947,
                "nominal_interval_coverage": 0.95,
            },
            [
                {
                    "method_id": "patch_feature_ridge",
                    "provenance": "RAN",
                    "heldout_gene_pearson_mean": 0.18,
                    "heldout_spot_pearson_mean": 0.86,
                    "same_split": True,
                    "calibrated_interval": False,
                }
            ],
            gate3_complete=True,
        )
        self.assertEqual(table["claim_status"], "locked")
        self.assertEqual(table["missing_claim_evidence"], ["license_review"])
        self.assertEqual(table["rows"][1]["provenance"], "RAN")


if __name__ == "__main__":
    unittest.main()
