from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .metrics import pearson
from .real_smoke import (
    RealSmokeConfig,
    _as_numpy,
    _close,
    _mean_pearson_cols,
    _mean_pearson_rows,
    _read_h5ad,
    _read_patch_images,
    _resolve_h5,
    _sample_pairs,
    _select_genes,
    default_expression_path,
    default_patches_path,
    run_real_data_smoke,
)


def track_root(anchor: str | Path | None = None) -> Path:
    here = Path(anchor or __file__).resolve()
    if here.is_file():
        here = here.parent
    for candidate in (here, *here.parents):
        if (candidate / "CLAIM_LEDGER.md").is_file() and (candidate / "experiments").is_dir():
            return candidate
    raise FileNotFoundError("could not find pathscope track root")


def default_reference_context_path() -> Path:
    return track_root() / "experiments" / "gate2_baseline_comparison" / "reference_rows.json"


def default_real_baseline_dir() -> Path:
    return track_root() / "experiments" / "gate2_real_baseline"


def patch_grid_features(images: Any, grid_size: int = 6) -> np.ndarray:
    """Return deterministic, non-pretrained patch descriptors for point predictors."""
    if grid_size < 1:
        raise ValueError("grid_size must be positive")
    arr = np.asarray(images, dtype=np.float32) / 255.0
    if arr.ndim != 4 or arr.shape[-1] != 3:
        raise ValueError("images must have shape n x height x width x 3")
    rgb_mean = arr.mean(axis=(1, 2))
    rgb_std = arr.std(axis=(1, 2))
    gray = arr.mean(axis=3)
    h, w = gray.shape[1:]
    y_edges = np.linspace(0, h, grid_size + 1, dtype=int)
    x_edges = np.linspace(0, w, grid_size + 1, dtype=int)
    cells: list[np.ndarray] = []
    for y0, y1 in zip(y_edges[:-1], y_edges[1:]):
        for x0, x1 in zip(x_edges[:-1], x_edges[1:]):
            block = gray[:, y0:y1, x0:x1]
            cells.append(block.mean(axis=(1, 2)) if block.size else np.zeros(gray.shape[0], dtype=np.float32))
    grid = np.stack(cells, axis=1).astype(np.float32)
    return np.concatenate([rgb_mean, rgb_std, grid], axis=1).astype(np.float32)


def _safe_gene_names(adata: Any, gene_idx: list[int]) -> list[str]:
    return [str(adata.var_names[int(i)]) for i in gene_idx]


def _safe_spot_names(adata: Any, obs_idx: list[int], test_mask: np.ndarray) -> list[str]:
    return [str(adata.obs_names[int(obs_idx[int(i)])]) for i in np.where(test_mask)[0]]


def _fit_predictors(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray, ridge_alpha: float, knn_neighbors: int) -> dict[str, np.ndarray]:
    k = min(max(1, int(knn_neighbors)), int(train_x.shape[0]))
    predictors = {
        "patch_feature_ridge": make_pipeline(StandardScaler(), Ridge(alpha=float(ridge_alpha))),
        "patch_feature_knn": make_pipeline(StandardScaler(), KNeighborsRegressor(n_neighbors=k, weights="distance")),
    }
    predictions: dict[str, np.ndarray] = {}
    for method_id, model in predictors.items():
        model.fit(train_x, train_y)
        predictions[method_id] = np.asarray(model.predict(test_x), dtype=np.float32)
    return predictions


def point_prediction_row(method_id: str, truth: np.ndarray, pred: np.ndarray, source: str, note: str) -> dict[str, object]:
    return {
        "method_id": method_id,
        "provenance": "RAN",
        "heldout_gene_pearson_mean": round(_mean_pearson_cols(truth, pred), 6),
        "heldout_spot_pearson_mean": round(_mean_pearson_rows(truth, pred), 6),
        "interval_coverage": None,
        "nominal_interval_coverage": None,
        "reported_gene_pearson_range": None,
        "same_split": True,
        "calibrated_interval": False,
        "source": source,
        "note": note,
    }


def _per_gene_pearsons(truth: np.ndarray, pred: np.ndarray) -> list[float]:
    return [round(pearson(list(map(float, truth[:, j])), list(map(float, pred[:, j]))), 6) for j in range(truth.shape[1])]


def run_same_split_baseline(
    *,
    config: RealSmokeConfig,
    reference_context_path: Path | None = None,
    out_dir: Path | None = None,
    grid_size: int = 6,
    baseline_ridge_alpha: float = 10.0,
    knn_neighbors: int = 8,
    local_smoke_metrics: dict[str, float] | None = None,
) -> dict[str, object]:
    config.validate()
    patches_path = _resolve_h5(config.patches_path)
    adata = _read_h5ad(config.expression_path)
    try:
        obs_idx, patch_idx = _sample_pairs(adata.obs_names, patches_path, config.max_spots, config.seed)
        gene_idx = _select_genes(adata, obs_idx, config.max_genes, config.max_gene_scan)
        y = np.log1p(_as_numpy(adata[obs_idx, gene_idx].X))
        x = patch_grid_features(_read_patch_images(patches_path, patch_idx), grid_size=grid_size)
        positions = np.arange(len(obs_idx))
        test_mask = positions % config.test_stride == 0
        if int(test_mask.sum()) < 2 or int((~test_mask).sum()) < 2:
            raise ValueError("held-out split is too small")
        truth = y[test_mask]
        predictions = _fit_predictors(x[~test_mask], y[~test_mask], x[test_mask], baseline_ridge_alpha, knn_neighbors)
        source = "same_split_baseline"
        rows = [
            point_prediction_row(
                method_id=method_id,
                truth=truth,
                pred=pred,
                source=source,
                note="same train/held-out spots as the calibrated local smoke; point predictions only, no interval calibration emitted",
            )
            for method_id, pred in predictions.items()
        ]
        context_rows: list[dict[str, object]] = []
        if reference_context_path is not None and reference_context_path.exists():
            context_rows = json.loads(reference_context_path.read_text(encoding="utf-8"))
        local_metrics = local_smoke_metrics if local_smoke_metrics is not None else run_real_data_smoke(config).metrics
        diagnostics = {
            method_id: {
                "per_gene_pearson": _per_gene_pearsons(truth, pred),
                "mean_absolute_error": round(float(np.mean(np.abs(truth - pred))), 6),
            }
            for method_id, pred in predictions.items()
        }
        payload: dict[str, object] = {
            "provenance": "RAN",
            "same_split": True,
            "config": {
                "expression_path": str(config.expression_path),
                "patches_path": str(config.patches_path),
                "max_spots": config.max_spots,
                "max_genes": config.max_genes,
                "max_gene_scan": config.max_gene_scan,
                "test_stride": config.test_stride,
                "seed": config.seed,
                "grid_size": grid_size,
                "baseline_ridge_alpha": baseline_ridge_alpha,
                "knn_neighbors": min(max(1, int(knn_neighbors)), int((~test_mask).sum())),
            },
            "split": {
                "paired_spot_count": int(len(obs_idx)),
                "train_spot_count": int((~test_mask).sum()),
                "heldout_spot_count": int(test_mask.sum()),
                "gene_count": int(len(gene_idx)),
                "gene_names": _safe_gene_names(adata, gene_idx),
                "heldout_spot_names": _safe_spot_names(adata, obs_idx, test_mask),
            },
            "local_calibrated_smoke_metrics": local_metrics,
            "rows": rows,
            "reference_rows": [*rows, *context_rows],
            "diagnostics": diagnostics,
        }
        if out_dir is not None:
            out_dir.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(
                out_dir / "point_predictions.npz",
                truth=truth.astype(np.float32),
                **{method_id: pred.astype(np.float32) for method_id, pred in predictions.items()},
                obs_indices=np.asarray(obs_idx, dtype=np.int64),
                patch_indices=np.asarray(patch_idx, dtype=np.int64),
                gene_indices=np.asarray(gene_idx, dtype=np.int64),
                test_mask=test_mask.astype(bool),
            )
            (out_dir / "reference_rows.json").write_text(json.dumps(payload["reference_rows"], indent=2, sort_keys=True) + "\n", encoding="utf-8")
            payload["reference_rows_path"] = str(out_dir / "reference_rows.json")
            payload["predictions_path"] = str(out_dir / "point_predictions.npz")
            (out_dir / "baseline_metrics.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload
    finally:
        _close(adata)


def default_config() -> RealSmokeConfig:
    return RealSmokeConfig(expression_path=default_expression_path(), patches_path=default_patches_path())
