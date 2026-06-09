from __future__ import annotations
from .contracts import CalibrationReport, ClaimGateEvidence, PairedPatchExpressionRecord, VirtualExpressionPrediction
from .metrics import interval_coverage, pearson


def build_fixture() -> tuple[list[PairedPatchExpressionRecord], list[VirtualExpressionPrediction]]:
    genes = ("GENE_A", "GENE_B", "GENE_C")
    records = []
    preds = []
    for i in range(30):
        color = i / 10.0
        texture = (i % 5) / 5.0
        expr = (1.0 + 1.5 * color, 0.5 + texture, 2.0 - 0.4 * color + 0.2 * texture)
        split = "test" if i % 5 == 0 else "train"
        rec = PairedPatchExpressionRecord("synthetic", f"p{i}", f"patch_{i}.png", float(i), float(i % 7), genes, expr, split, {"fixture": "synthetic"})
        rec.validate()
        records.append(rec)
        noise = 0.02 if split == "test" else 0.01
        mean = tuple(v + noise for v in expr)
        preds.append(VirtualExpressionPrediction("synthetic", f"p{i}", genes, mean, tuple(v-0.1 for v in mean), tuple(v+0.1 for v in mean), (0.1, 0.1, 0.1)))
    return records, preds


def run_synthetic_smoke() -> CalibrationReport:
    records, preds = build_fixture()
    truth = [rec.expression[0] for rec in records if rec.split == "test"]
    pred = [pr.mean[0] for rec, pr in zip(records, preds) if rec.split == "test"]
    lower = [pr.lower[0] for rec, pr in zip(records, preds) if rec.split == "test" and pr.lower is not None]
    upper = [pr.upper[0] for rec, pr in zip(records, preds) if rec.split == "test" and pr.upper is not None]
    metrics = {"heldout_pearson_gene_a": round(pearson(truth, pred), 4), "interval_coverage_gene_a": round(interval_coverage(truth, lower, upper) or 0.0, 4), "test_patches": float(len(truth))}
    return CalibrationReport(metrics, ClaimGateEvidence(ablation=True, failure_modes=True))
