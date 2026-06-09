import unittest
from pathscope_st import ExternalRunSummary, run_synthetic_smoke
from pathscope_st.metrics import interval_coverage

class SmokeTests(unittest.TestCase):
    def test_synthetic_prediction_is_calibrated_contract_only(self):
        report = run_synthetic_smoke()
        self.assertEqual(report.claim_status.value, "locked")
        self.assertGreater(report.metrics["heldout_pearson_gene_a"], 0.99)
        self.assertGreaterEqual(report.metrics["interval_coverage_gene_a"], 0.9)
    def test_intervals_required_for_coverage(self):
        self.assertIsNone(interval_coverage([1.0], None, None))
    def test_external_adapter_is_file_boundary(self):
        summary = ExternalRunSummary("r", ("python", "runner.py"), ("patch_id", "gene", "prediction"), {"sha": "x"})
        summary.validate_file_boundary()

if __name__ == "__main__":
    unittest.main()
