import unittest
from pathscope_st import ClaimGateEvidence, ClaimStatus, PairedPatchExpressionRecord, evaluate_claim_gate

class ContractTests(unittest.TestCase):
    def test_record_schema_fails_loudly(self):
        rec = PairedPatchExpressionRecord("s", "p", "uri", 0, 0, ("g1",), (1.0, 2.0), "test", {})
        with self.assertRaises(ValueError):
            rec.validate()
    def test_claim_lock(self):
        self.assertEqual(evaluate_claim_gate(ClaimGateEvidence(public_data_smoke=True)), ClaimStatus.LOCKED)

if __name__ == "__main__":
    unittest.main()
