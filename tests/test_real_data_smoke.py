import tempfile
import unittest
from pathlib import Path

from pathscope_st.contracts import ClaimGateEvidence
from pathscope_st.data_paths import find_repo_root, processed_data_path
from pathscope_st.real_smoke import RealSmokeConfig, RealSmokeResult, _patch_features


class RealDataSmokeUnitTests(unittest.TestCase):
    def test_result_keeps_claim_locked(self):
        result = RealSmokeResult({"heldout_gene_pearson_mean": 0.1}, ClaimGateEvidence(public_data_smoke=True))
        self.assertEqual(result.claim_status.value, "locked")
        self.assertEqual(
            result.to_jsonable()["missing_claim_evidence"],
            ["baseline_comparison", "ablation", "failure_modes", "license_review"],
        )

    def test_config_rejects_tiny_spot_count(self):
        with self.assertRaises(ValueError):
            RealSmokeConfig(expression_path=Path("expression.h5ad"), patches_path=Path("patches.h5"), max_spots=3).validate()

    def test_repo_root_path_resolution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "data" / "processed").mkdir(parents=True)
            anchor = root / "src" / "package" / "module.py"
            anchor.parent.mkdir(parents=True)
            anchor.touch()
            self.assertEqual(find_repo_root(anchor), root)
            self.assertEqual(
                processed_data_path("fixture_card", anchor=anchor),
                root / "data" / "processed" / "fixture_card",
            )

    def test_patch_feature_shape(self):
        import numpy as np

        images = np.zeros((2, 4, 4, 3), dtype=np.uint8)
        self.assertEqual(_patch_features(images).shape, (2, 8))


if __name__ == "__main__":
    unittest.main()
