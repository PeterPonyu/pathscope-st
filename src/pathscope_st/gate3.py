from __future__ import annotations

from pathlib import Path
from typing import Any

from .contracts import ClaimGateEvidence, evaluate_claim_gate
from .metrics import interval_coverage, pearson
from .real_smoke import (
    RealSmokeConfig,
    _as_numpy,
    _close,
    _fit_ridge,
    _mean_pearson_cols,
    _mean_pearson_rows,
    _patch_features,
    _read_h5ad,
    _read_patch_images,
    _resolve_h5,
    _sample_pairs,
    _select_genes,
    default_expression_path,
    default_patches_path,
)


def _prediction_bundle(config: RealSmokeConfig) -> dict[str, Any]:
    import numpy as np

    config.validate()
    patches_path = _resolve_h5(config.patches_path)
    adata = _read_h5ad(config.expression_path)
    try:
        obs_idx, patch_idx = _sample_pairs(adata.obs_names, patches_path, config.max_spots, config.seed)
        gene_idx = _select_genes(adata, obs_idx, config.max_genes, config.max_gene_scan)
        y = np.log1p(_as_numpy(adata[obs_idx, gene_idx].X))
        features = _patch_features(_read_patch_images(patches_path, patch_idx))
        positions = np.arange(len(obs_idx))
        test_mask = positions % config.test_stride == 0
        if int(test_mask.sum()) < 2 or int((~test_mask).sum()) < 2:
            raise ValueError("held-out split is too small")
        pred, lower, upper = _fit_ridge(features[~test_mask], y[~test_mask], features[test_mask], config.ridge_alpha)
        truth = y[test_mask]
        train_truth = y[~test_mask]
        return {
            "truth": truth,
            "prediction": pred,
            "lower": lower,
            "upper": upper,
            "test_features": features[test_mask],
            "train_truth": train_truth,
            "gene_names": [str(adata.var_names[idx]) for idx in gene_idx],
            "paired_spot_count": len(obs_idx),
            "heldout_spot_count": int(test_mask.sum()),
        }
    finally:
        _close(adata)


def _coverage(truth: Any, lower: Any, upper: Any) -> float:
    return float(
        interval_coverage(
            list(map(float, truth.reshape(-1))),
            list(map(float, lower.reshape(-1))),
            list(map(float, upper.reshape(-1))),
        )
        or 0.0
    )


def _spot_pearsons(truth: Any, pred: Any) -> list[float]:
    return [pearson(list(map(float, row_truth)), list(map(float, row_pred))) for row_truth, row_pred in zip(truth, pred)]


def _gene_pearsons(truth: Any, pred: Any) -> list[float]:
    return [pearson(list(map(float, truth[:, j])), list(map(float, pred[:, j]))) for j in range(truth.shape[1])]


def build_gate3_table(*, diagnostics: dict[str, float]) -> dict[str, object]:
    evidence = ClaimGateEvidence(
        public_data_smoke=True,
        baseline_comparison=True,
        ablation=True,
        failure_modes=True,
    )
    calibrated_error = abs(float(diagnostics["calibrated_interval_coverage"]) - 0.95)
    uncalibrated_error = abs(float(diagnostics["uncalibrated_point_interval_coverage"]) - 0.95)
    return {
        "ablation": {
            "removed_component": "calibrated_interval_width",
            "calibrated_interval_coverage": round(float(diagnostics["calibrated_interval_coverage"]), 6),
            "uncalibrated_point_interval_coverage": round(
                float(diagnostics["uncalibrated_point_interval_coverage"]),
                6,
            ),
            "calibrated_absolute_coverage_error": round(calibrated_error, 6),
            "uncalibrated_absolute_coverage_error": round(uncalibrated_error, 6),
            "coverage_error_increase": round(uncalibrated_error - calibrated_error, 6),
        },
        "failure_mode": {
            "mode": "low_information_patch_or_hard_gene",
            "heldout_gene_pearson_mean": round(float(diagnostics["heldout_gene_pearson_mean"]), 6),
            "hard_gene_bottom_quartile_pearson_mean": round(
                float(diagnostics["hard_gene_bottom_quartile_pearson_mean"]),
                6,
            ),
            "hard_gene_pearson_floor": round(float(diagnostics["hard_gene_pearson_floor"]), 6),
            "low_information_patch_pearson_floor": round(
                float(diagnostics["low_information_patch_pearson_floor"]),
                6,
            ),
        },
        "claim_status": evaluate_claim_gate(evidence).value,
        "missing_claim_evidence": list(evidence.missing()),
    }


def run_gate3_analysis(
    *,
    expression_path: Path | None = None,
    patches_path: Path | None = None,
    max_spots: int = 768,
    max_genes: int = 32,
    max_gene_scan: int = 4096,
    test_stride: int = 5,
    ridge_alpha: float = 1.0,
    seed: int = 23,
) -> dict[str, object]:
    import numpy as np

    config = RealSmokeConfig(
        expression_path=expression_path or default_expression_path(),
        patches_path=patches_path or default_patches_path(),
        max_spots=max_spots,
        max_genes=max_genes,
        max_gene_scan=max_gene_scan,
        test_stride=test_stride,
        ridge_alpha=ridge_alpha,
        seed=seed,
    )
    bundle = _prediction_bundle(config)
    truth = bundle["truth"]
    pred = bundle["prediction"]
    spot_values = np.asarray(_spot_pearsons(truth, pred), dtype=float)
    gene_values = np.asarray(_gene_pearsons(truth, pred), dtype=float)
    contrast = np.asarray(bundle["test_features"][:, -1], dtype=float)
    low_info_threshold = float(np.quantile(contrast, 0.25))
    low_info_spots = spot_values[contrast <= low_info_threshold]
    hard_threshold = float(np.quantile(gene_values, 0.25))
    hard_genes = gene_values[gene_values <= hard_threshold]
    diagnostics = {
        "paired_spot_count": float(bundle["paired_spot_count"]),
        "heldout_spot_count": float(bundle["heldout_spot_count"]),
        "gene_count": float(len(bundle["gene_names"])),
        "heldout_gene_pearson_mean": round(_mean_pearson_cols(truth, pred), 6),
        "heldout_spot_pearson_mean": round(_mean_pearson_rows(truth, pred), 6),
        "calibrated_interval_coverage": round(_coverage(truth, bundle["lower"], bundle["upper"]), 6),
        "uncalibrated_point_interval_coverage": round(_coverage(truth, pred, pred), 6),
        "hard_gene_bottom_quartile_pearson_mean": round(float(np.mean(hard_genes)), 6),
        "hard_gene_pearson_floor": round(float(np.min(gene_values)), 6),
        "low_information_patch_pearson_mean": round(float(np.mean(low_info_spots)), 6),
        "low_information_patch_pearson_floor": round(float(np.min(low_info_spots)), 6),
    }
    table = build_gate3_table(diagnostics=diagnostics)
    return {
        "diagnostics": diagnostics,
        "gate3": table,
    }
