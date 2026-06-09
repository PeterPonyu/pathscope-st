from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .contracts import ClaimGateEvidence, ClaimStatus, evaluate_claim_gate
from .data_paths import processed_data_path
from .metrics import interval_coverage, pearson


def default_expression_path() -> Path:
    return processed_data_path("pathscope_coad_visium_he", "expression.h5ad")


def default_patches_path() -> Path:
    return processed_data_path("pathscope_coad_visium_he", "patches.h5")


@dataclass(frozen=True)
class RealSmokeConfig:
    expression_path: Path = field(default_factory=default_expression_path)
    patches_path: Path = field(default_factory=default_patches_path)
    max_spots: int = 768
    max_genes: int = 32
    max_gene_scan: int = 4096
    test_stride: int = 5
    ridge_alpha: float = 1.0
    seed: int = 23

    def validate(self) -> None:
        if self.max_spots < 20:
            raise ValueError("max_spots must be at least 20")
        if self.max_genes < 2:
            raise ValueError("max_genes must be at least two")
        if self.max_gene_scan < self.max_genes:
            raise ValueError("max_gene_scan must be at least max_genes")
        if self.test_stride < 2:
            raise ValueError("test_stride must be at least two")
        if self.ridge_alpha < 0.0:
            raise ValueError("ridge_alpha must be non-negative")


@dataclass(frozen=True)
class RealSmokeResult:
    metrics: dict[str, float]
    evidence: ClaimGateEvidence

    @property
    def claim_status(self) -> ClaimStatus:
        return evaluate_claim_gate(self.evidence)

    def to_jsonable(self) -> dict[str, object]:
        return {
            "metrics": self.metrics,
            "claim_status": self.claim_status.value,
            "missing_claim_evidence": list(self.evidence.missing()),
        }


def _resolve_h5ad(path: str | Path) -> Path:
    resolved = Path(path).expanduser()
    if not resolved.exists():
        raise FileNotFoundError(f"AnnData file does not exist: {resolved}")
    if resolved.suffix != ".h5ad":
        raise ValueError(f"expected a .h5ad file, got: {resolved}")
    return resolved


def _resolve_h5(path: str | Path) -> Path:
    resolved = Path(path).expanduser()
    if not resolved.exists():
        raise FileNotFoundError(f"patch file does not exist: {resolved}")
    if resolved.suffix not in {".h5", ".hdf5"}:
        raise ValueError(f"expected an HDF5 patch file, got: {resolved}")
    return resolved


def _read_h5ad(path: str | Path) -> Any:
    try:
        import anndata as ad  # type: ignore[import-untyped]
    except ImportError as exc:  # pragma: no cover - depends on optional real-data env
        raise RuntimeError("real-data smoke requires the optional anndata package") from exc
    return ad.read_h5ad(_resolve_h5ad(path), backed="r")


def _close(adata: Any) -> None:
    close = getattr(getattr(adata, "file", None), "close", None)
    if close is not None:
        close()


def _as_numpy(matrix: Any) -> Any:
    import numpy as np

    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32)


def _sample_pairs(obs_names: Any, patches_path: Path, max_spots: int, seed: int) -> tuple[list[int], list[int]]:
    import h5py  # type: ignore[import-untyped]
    import numpy as np

    with h5py.File(patches_path, "r") as h:
        if "barcode" not in h or "img" not in h:
            raise KeyError("patch HDF5 must contain 'barcode' and 'img' datasets")
        raw_barcodes = h["barcode"][:]
        valid = h["valid"][:] if "valid" in h else np.ones(len(raw_barcodes), dtype=bool)
    barcodes = [b.decode() if isinstance(b, bytes) else str(b) for b in raw_barcodes]
    patch_by_barcode = {barcode: idx for idx, barcode in enumerate(barcodes) if bool(valid[idx])}
    pairs = [(idx, patch_by_barcode[str(name)]) for idx, name in enumerate(obs_names) if str(name) in patch_by_barcode]
    if len(pairs) < 20:
        raise ValueError("not enough expression/patch pairs for real-data smoke")
    rng = np.random.default_rng(seed)
    if len(pairs) > max_spots:
        selected = rng.choice(len(pairs), size=max_spots, replace=False)
        pairs = [pairs[int(i)] for i in selected]
    pairs.sort(key=lambda item: item[0])
    return [p[0] for p in pairs], [p[1] for p in pairs]


def _read_patch_images(patches_path: Path, patch_idx: list[int]) -> Any:
    import h5py  # type: ignore[import-untyped]
    import numpy as np

    with h5py.File(patches_path, "r") as h:
        images = np.stack([h["img"][int(i)] for i in patch_idx], axis=0)
    return images


def _patch_features(images: Any) -> Any:
    import numpy as np

    arr = np.asarray(images, dtype=np.float32) / 255.0
    means = arr.mean(axis=(1, 2))
    stds = arr.std(axis=(1, 2))
    gray = arr.mean(axis=3)
    darkness = 1.0 - gray.mean(axis=(1, 2), keepdims=False)[:, None]
    contrast = gray.std(axis=(1, 2), keepdims=False)[:, None]
    return np.concatenate([means, stds, darkness, contrast], axis=1)


def _standardize_train_test(train: Any, test: Any) -> tuple[Any, Any]:
    import numpy as np

    mean = train.mean(axis=0, keepdims=True)
    scale = train.std(axis=0, keepdims=True)
    scale = np.where(scale > 0.0, scale, 1.0)
    return (train - mean) / scale, (test - mean) / scale


def _with_intercept(matrix: Any) -> Any:
    import numpy as np

    return np.concatenate([np.ones((matrix.shape[0], 1), dtype=np.float32), matrix.astype(np.float32)], axis=1)


def _fit_ridge(train_x: Any, train_y: Any, test_x: Any, alpha: float) -> tuple[Any, Any, Any]:
    import numpy as np

    train_x, test_x = _standardize_train_test(train_x, test_x)
    train_x = _with_intercept(train_x)
    test_x = _with_intercept(test_x)
    penalty = np.eye(train_x.shape[1], dtype=np.float64) * float(alpha)
    penalty[0, 0] = 0.0
    beta = np.linalg.solve(train_x.T @ train_x + penalty, train_x.T @ train_y)
    train_pred = train_x @ beta
    test_pred = test_x @ beta
    residual = train_y - train_pred
    sigma = residual.std(axis=0, ddof=min(train_x.shape[1], max(1, train_x.shape[0] - 1)))
    sigma = np.where(sigma > 1e-6, sigma, 1e-6)
    return test_pred, test_pred - 1.96 * sigma, test_pred + 1.96 * sigma


def _select_genes(adata: Any, obs_idx: list[int], max_genes: int, max_scan: int) -> list[int]:
    import numpy as np

    scan_count = min(int(adata.n_vars), max_scan)
    scan_idx = np.unique(np.linspace(0, int(adata.n_vars) - 1, scan_count, dtype=int))
    matrix = np.log1p(_as_numpy(adata[obs_idx, scan_idx].X))
    variance = matrix.var(axis=0)
    mean = matrix.mean(axis=0)
    usable = np.where(mean > 0.0, variance, -1.0)
    ranked = np.lexsort((scan_idx, -usable))
    selected = [int(scan_idx[i]) for i in ranked[:max_genes] if usable[i] >= 0.0]
    if len(selected) < 2:
        raise ValueError("not enough expressed genes in scanned expression matrix")
    return selected


def _mean_pearson_rows(truth: Any, pred: Any) -> float:
    import numpy as np

    vals = []
    for row_truth, row_pred in zip(truth, pred):
        vals.append(pearson(list(map(float, row_truth)), list(map(float, row_pred))))
    return float(np.mean(vals)) if vals else 0.0


def _mean_pearson_cols(truth: Any, pred: Any) -> float:
    import numpy as np

    vals = []
    for j in range(truth.shape[1]):
        vals.append(pearson(list(map(float, truth[:, j])), list(map(float, pred[:, j]))))
    return float(np.mean(vals)) if vals else 0.0


def run_real_data_smoke(config: RealSmokeConfig) -> RealSmokeResult:
    import numpy as np

    config.validate()
    patches_path = _resolve_h5(config.patches_path)
    adata = _read_h5ad(config.expression_path)
    try:
        obs_idx, patch_idx = _sample_pairs(adata.obs_names, patches_path, config.max_spots, config.seed)
        gene_idx = _select_genes(adata, obs_idx, config.max_genes, config.max_gene_scan)
        y = np.log1p(_as_numpy(adata[obs_idx, gene_idx].X))
        x = _patch_features(_read_patch_images(patches_path, patch_idx))
        positions = np.arange(len(obs_idx))
        test_mask = positions % config.test_stride == 0
        if int(test_mask.sum()) < 2 or int((~test_mask).sum()) < 2:
            raise ValueError("held-out split is too small")
        pred, lower, upper = _fit_ridge(x[~test_mask], y[~test_mask], x[test_mask], config.ridge_alpha)
        truth = y[test_mask]
        coverage = interval_coverage(list(map(float, truth.reshape(-1))), list(map(float, lower.reshape(-1))), list(map(float, upper.reshape(-1))))
        metrics = {
            "paired_spot_count": float(len(obs_idx)),
            "heldout_spot_count": float(int(test_mask.sum())),
            "gene_count": float(len(gene_idx)),
            "heldout_gene_pearson_mean": round(_mean_pearson_cols(truth, pred), 6),
            "heldout_spot_pearson_mean": round(_mean_pearson_rows(truth, pred), 6),
            "interval_coverage": round(float(coverage or 0.0), 6),
            "nominal_interval_coverage": 0.95,
            "mean_interval_width": round(float(np.mean(upper - lower)), 6),
            "mean_absolute_error": round(float(np.mean(np.abs(truth - pred))), 6),
        }
        return RealSmokeResult(metrics, ClaimGateEvidence(public_data_smoke=True))
    finally:
        _close(adata)
