import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from pathscope_st.contracts import ClaimGateEvidence
from pathscope_st.real_smoke import RealSmokeConfig, RealSmokeResult
from pathscope_st.result_wiring import emit_gate2_results


class ResultsContractWiringTests(unittest.TestCase):
    def test_vendored_contract_matches_recorded_sha256(self):
        package_dir = Path(__file__).resolve().parents[1] / "src" / "pathscope_st"
        contract_path = package_dir / "results_contract.py"
        expected = (package_dir / "results_contract.sha256").read_text(encoding="utf-8").strip().split()[0]
        observed = hashlib.sha256(contract_path.read_bytes()).hexdigest()
        self.assertEqual(observed, expected)

    def test_gate2_emits_contract_files(self):
        report = RealSmokeResult(
            metrics={
                "paired_spot_count": 10.0,
                "heldout_spot_count": 2.0,
                "gene_count": 4.0,
                "heldout_gene_pearson_mean": 0.269,
                "heldout_spot_pearson_mean": 0.874,
                "interval_coverage": 0.947,
                "nominal_interval_coverage": 0.95,
            },
            evidence=ClaimGateEvidence(public_data_smoke=True),
        )
        parity = {
            "rows": [
                {
                    "method_id": "local_calibrated_smoke",
                    "heldout_gene_pearson_mean": 0.269,
                    "heldout_spot_pearson_mean": 0.874,
                    "interval_coverage": 0.947,
                },
                {
                    "method_id": "simple_ridge",
                    "heldout_gene_pearson_mean": 0.180,
                    "heldout_spot_pearson_mean": 0.866,
                    "interval_coverage": None,
                },
            ],
            "differentiator": {
                "observed_interval_coverage": 0.947,
                "nominal_interval_coverage": 0.95,
                "absolute_coverage_error": 0.003,
            },
        }
        config = RealSmokeConfig(
            expression_path=Path("data/processed/example/expression.h5ad"),
            patches_path=Path("data/processed/example/patches.h5"),
            max_spots=30,
            max_genes=4,
            max_gene_scan=8,
        )
        with tempfile.TemporaryDirectory() as tmp:
            paths = emit_gate2_results(report, parity, config, results_dir=Path(tmp))
            metrics = json.loads(paths["metrics"].read_text(encoding="utf-8"))
            self.assertEqual(metrics["project"], "pathscope-st")
            self.assertEqual(metrics["metrics"]["gate2.observed_interval_coverage"], 0.947)
            self.assertEqual(metrics["metrics"]["gate2.simple_ridge.interval_coverage"], None)
            self.assertTrue(paths["run_metadata"].is_file())


if __name__ == "__main__":
    unittest.main()
