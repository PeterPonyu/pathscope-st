from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .contracts import ClaimGateEvidence, ClaimStatus, evaluate_claim_gate

PROVENANCE_RAN = "RAN"
PROVENANCE_REFERENCE_REPORTED = "REFERENCE_REPORTED"


@dataclass(frozen=True)
class ReferenceMetricRow:
    method_id: str
    provenance: str
    heldout_gene_pearson_mean: float | None = None
    heldout_spot_pearson_mean: float | None = None
    interval_coverage: float | None = None
    nominal_interval_coverage: float | None = None
    reported_gene_pearson_range: str | None = None
    same_split: bool | None = None
    calibrated_interval: bool | None = None
    source: str | None = None
    note: str | None = None

    def validate(self) -> None:
        if not self.method_id:
            raise ValueError("method_id is required")
        if self.provenance not in {PROVENANCE_RAN, PROVENANCE_REFERENCE_REPORTED}:
            raise ValueError("provenance must be RAN or REFERENCE_REPORTED")

    def to_jsonable(self) -> dict[str, object]:
        self.validate()
        return {
            "method_id": self.method_id,
            "provenance": self.provenance,
            "heldout_gene_pearson_mean": self.heldout_gene_pearson_mean,
            "heldout_spot_pearson_mean": self.heldout_spot_pearson_mean,
            "interval_coverage": self.interval_coverage,
            "nominal_interval_coverage": self.nominal_interval_coverage,
            "reported_gene_pearson_range": self.reported_gene_pearson_range,
            "same_split": self.same_split,
            "calibrated_interval": self.calibrated_interval,
            "source": self.source,
            "note": self.note,
        }


def _as_reference_rows(rows: list[dict[str, Any]]) -> list[ReferenceMetricRow]:
    parsed: list[ReferenceMetricRow] = []
    for row in rows:
        parsed.append(
            ReferenceMetricRow(
                method_id=str(row["method_id"]),
                provenance=str(row["provenance"]),
                heldout_gene_pearson_mean=row.get("heldout_gene_pearson_mean"),
                heldout_spot_pearson_mean=row.get("heldout_spot_pearson_mean"),
                interval_coverage=row.get("interval_coverage"),
                nominal_interval_coverage=row.get("nominal_interval_coverage"),
                reported_gene_pearson_range=row.get("reported_gene_pearson_range"),
                same_split=row.get("same_split"),
                calibrated_interval=row.get("calibrated_interval"),
                source=row.get("source"),
                note=row.get("note"),
            )
        )
    return parsed


def build_parity_table(
    real_metrics: dict[str, float],
    reference_rows: list[dict[str, Any]] | None = None,
    *,
    gate3_complete: bool = False,
) -> dict[str, object]:
    evidence = ClaimGateEvidence(
        public_data_smoke=True,
        baseline_comparison=True,
        ablation=gate3_complete,
        failure_modes=gate3_complete,
    )
    observed_coverage = float(real_metrics.get("interval_coverage", 0.0))
    nominal_coverage = float(real_metrics.get("nominal_interval_coverage", 0.95))
    local_row = ReferenceMetricRow(
        method_id="local_calibrated_smoke",
        provenance=PROVENANCE_RAN,
        heldout_gene_pearson_mean=float(real_metrics["heldout_gene_pearson_mean"]),
        heldout_spot_pearson_mean=float(real_metrics["heldout_spot_pearson_mean"]),
        interval_coverage=observed_coverage,
        nominal_interval_coverage=nominal_coverage,
        same_split=True,
        calibrated_interval=True,
        source="smoke-real",
        note="same held-out spots as the gate-1 real-data smoke",
    )
    parsed_rows = [local_row, *_as_reference_rows(reference_rows or [])]
    rows = [row.to_jsonable() for row in parsed_rows]
    return {
        "rows": rows,
        "differentiator": {
            "observed_interval_coverage": round(observed_coverage, 6),
            "nominal_interval_coverage": round(nominal_coverage, 6),
            "absolute_coverage_error": round(abs(observed_coverage - nominal_coverage), 6),
            "reference_interval_coverage_available": any(
                row.interval_coverage is not None for row in parsed_rows[1:]
            ),
        },
        "claim_status": evaluate_claim_gate(evidence).value,
        "missing_claim_evidence": list(evidence.missing()),
    }


def assert_locked_after_gate2(table: dict[str, object]) -> ClaimStatus:
    status = ClaimStatus(str(table["claim_status"]))
    if status != ClaimStatus.LOCKED:
        raise AssertionError("gate-2 parity must not unlock claims")
    return status
