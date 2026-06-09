import unittest

from pathscope_st.gate3 import build_gate3_table


class Gate3AnalysisTests(unittest.TestCase):
    def test_calibration_ablation_and_failure_modes_leave_only_license_review(self):
        table = build_gate3_table(
            diagnostics={
                "calibrated_interval_coverage": 0.947,
                "uncalibrated_point_interval_coverage": 0.0,
                "heldout_gene_pearson_mean": 0.269,
                "hard_gene_bottom_quartile_pearson_mean": 0.202,
                "hard_gene_pearson_floor": 0.17,
                "low_information_patch_pearson_floor": 0.657,
            }
        )
        self.assertEqual(table["claim_status"], "locked")
        self.assertEqual(table["missing_claim_evidence"], ["license_review"])
        self.assertGreater(table["ablation"]["coverage_error_increase"], 0.9)
        self.assertEqual(table["failure_mode"]["hard_gene_pearson_floor"], 0.17)


if __name__ == "__main__":
    unittest.main()
