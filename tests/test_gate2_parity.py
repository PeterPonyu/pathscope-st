import unittest

from pathscope_st.parity import build_parity_table


class Gate2ParityTests(unittest.TestCase):
    def test_gate2_table_keeps_claims_locked_and_records_calibration(self):
        table = build_parity_table(
            {
                "heldout_gene_pearson_mean": 0.269,
                "heldout_spot_pearson_mean": 0.874,
                "interval_coverage": 0.947,
                "nominal_interval_coverage": 0.95,
            },
            [
                {
                    "method_id": "reported_reference",
                    "provenance": "REFERENCE_REPORTED",
                    "heldout_gene_pearson_mean": 0.48,
                    "note": "not run on this held-out split",
                }
            ],
        )
        self.assertEqual(table["claim_status"], "locked")
        self.assertEqual(table["rows"][0]["provenance"], "RAN")
        self.assertAlmostEqual(table["differentiator"]["absolute_coverage_error"], 0.003)
        self.assertIn("ablation", table["missing_claim_evidence"])


if __name__ == "__main__":
    unittest.main()
