from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .contracts import ClaimGateEvidence
from .real_smoke import RealSmokeConfig, RealSmokeResult
from .results_contract import dataset_card_id, write_results

PROJECT_ID = "pathscope-st"


def _numeric_items(payload: Mapping[str, Any], prefix: str = "") -> dict[str, float | None]:
    metrics: dict[str, float | None] = {}
    for key, value in payload.items():
        name = f"{prefix}.{key}" if prefix else str(key)
        if value is None:
            metrics[name] = None
        elif isinstance(value, bool):
            continue
        elif isinstance(value, (int, float)):
            metrics[name] = float(value)
        elif isinstance(value, Mapping):
            metrics.update(_numeric_items(value, name))
    return metrics


def _config_paths(config: RealSmokeConfig) -> list[str]:
    return [str(config.expression_path), str(config.patches_path)]


def _metadata(config: RealSmokeConfig, report: RealSmokeResult | None, *, notes: str) -> dict[str, Any]:
    metrics = report.metrics if report is not None else {}
    return {
        "dataset_paths": _config_paths(config),
        "n_obs": metrics.get("paired_spot_count"),
        "n_vars": metrics.get("gene_count"),
        "seed": config.seed,
        "deterministic": True,
        "num_threads": 1,
        "reproducibility_level": "seeded",
        "normalization": {"applied": True, "method": "log1p"},
        "interpretability": {
            "claim_scope": "calibrated_uncertainty_not_best_accuracy",
            "calibration_target": "nominal_interval_coverage_0.95",
        },
        "notes": notes,
        "provenance": {
            "max_spots": config.max_spots,
            "max_genes": config.max_genes,
            "max_gene_scan": config.max_gene_scan,
            "test_stride": config.test_stride,
            "ridge_alpha": config.ridge_alpha,
        },
    }


def emit_real_smoke_results(
    report: RealSmokeResult,
    config: RealSmokeConfig,
    *,
    results_dir: Path | None = None,
    outputs: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    return write_results(
        PROJECT_ID,
        dataset_card_id(_config_paths(config)),
        report.metrics,
        outputs=outputs,
        run_metadata=_metadata(config, report, notes="real-data calibrated smoke emitted through the vendored results contract"),
        results_dir=results_dir,
    )


def emit_gate2_results(
    report: RealSmokeResult,
    parity_table: Mapping[str, Any],
    config: RealSmokeConfig,
    *,
    results_dir: Path | None = None,
    outputs: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    metrics = {f"real_smoke.{key}": value for key, value in report.metrics.items()}
    metrics.update(_numeric_items(parity_table.get("differentiator", {}), "gate2"))
    for row in parity_table.get("rows", []):
        if not isinstance(row, Mapping):
            continue
        method = str(row.get("method_id", "method")).replace(".", "_")
        for key in ("heldout_gene_pearson_mean", "heldout_spot_pearson_mean", "interval_coverage"):
            value = row.get(key)
            if value is None or isinstance(value, (int, float)):
                metrics[f"gate2.{method}.{key}"] = value
    return write_results(
        PROJECT_ID,
        dataset_card_id(_config_paths(config)),
        metrics,
        outputs=outputs,
        run_metadata=_metadata(config, report, notes="gate-2 same-split parity emitted through the vendored results contract"),
        results_dir=results_dir,
    )


def emit_gate3_results(
    payload: Mapping[str, Any],
    config: RealSmokeConfig,
    *,
    results_dir: Path | None = None,
    outputs: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    real_payload = payload.get("real_smoke", {})
    real_metrics = real_payload.get("metrics", {}) if isinstance(real_payload, Mapping) else {}
    diagnostics = payload.get("diagnostics", {})
    if not real_metrics and isinstance(diagnostics, Mapping):
        real_metrics = diagnostics
    metrics = {f"real_smoke.{key}": value for key, value in real_metrics.items() if not isinstance(value, bool)}
    if isinstance(diagnostics, Mapping):
        metrics.update(_numeric_items(diagnostics, "gate3.diagnostics"))
    gate3_payload = payload.get("gate3", {})
    if isinstance(gate3_payload, Mapping):
        metrics.update(_numeric_items(gate3_payload.get("ablation", {}), "gate3.ablation"))
        metrics.update(_numeric_items(gate3_payload.get("failure_mode", {}), "gate3.failure_mode"))
    report = RealSmokeResult(dict(real_metrics), evidence=ClaimGateEvidence(public_data_smoke=True)) if isinstance(real_metrics, Mapping) else None
    return write_results(
        PROJECT_ID,
        dataset_card_id(_config_paths(config)),
        metrics,
        outputs=outputs,
        run_metadata=_metadata(config, report, notes="gate-3 calibration ablation/failure emitted through the vendored results contract"),
        results_dir=results_dir,
    )
