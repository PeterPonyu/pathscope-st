from __future__ import annotations

import argparse
import json
from pathlib import Path

from .claim_status import graduation_claim_status
from .contracts import ClaimGateEvidence, evaluate_claim_gate
from .real_smoke import RealSmokeConfig, default_expression_path, default_patches_path
from .smoke import run_synthetic_smoke


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pathscope-st")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("smoke-synthetic")
    real = sub.add_parser("smoke-real")
    real.add_argument("--expression-path", "--st-path", dest="expression_path", type=Path, default=None)
    real.add_argument("--patches-path", type=Path, default=None)
    real.add_argument("--max-spots", type=int, default=768)
    real.add_argument("--max-genes", type=int, default=32)
    real.add_argument("--max-gene-scan", type=int, default=4096)
    real.add_argument("--test-stride", type=int, default=5)
    real.add_argument("--ridge-alpha", type=float, default=1.0)
    real.add_argument("--seed", type=int, default=23)
    real.add_argument("--results-dir", type=Path, default=None)
    parity = sub.add_parser("gate2-parity")
    parity.add_argument("--expression-path", "--st-path", dest="expression_path", type=Path, default=None)
    parity.add_argument("--patches-path", type=Path, default=None)
    parity.add_argument("--max-spots", type=int, default=768)
    parity.add_argument("--max-genes", type=int, default=32)
    parity.add_argument("--max-gene-scan", type=int, default=4096)
    parity.add_argument("--test-stride", type=int, default=5)
    parity.add_argument("--ridge-alpha", type=float, default=1.0)
    parity.add_argument("--seed", type=int, default=23)
    parity.add_argument("--reference-path", type=Path, default=None)
    parity.add_argument("--out-path", type=Path, default=None)
    parity.add_argument("--grid-size", type=int, default=6)
    parity.add_argument("--baseline-ridge-alpha", type=float, default=10.0)
    parity.add_argument("--knn-neighbors", type=int, default=8)
    parity.add_argument("--results-dir", type=Path, default=None)
    gate3 = sub.add_parser("gate3-analysis")
    gate3.add_argument("--expression-path", "--st-path", dest="expression_path", type=Path, default=None)
    gate3.add_argument("--patches-path", type=Path, default=None)
    gate3.add_argument("--max-spots", type=int, default=768)
    gate3.add_argument("--max-genes", type=int, default=32)
    gate3.add_argument("--max-gene-scan", type=int, default=4096)
    gate3.add_argument("--test-stride", type=int, default=5)
    gate3.add_argument("--ridge-alpha", type=float, default=1.0)
    gate3.add_argument("--seed", type=int, default=23)
    gate3.add_argument("--out-path", type=Path, default=None)
    gate3.add_argument("--results-dir", type=Path, default=None)
    sub.add_parser("claim-status")
    args = parser.parse_args(argv)
    if args.command == "smoke-synthetic":
        report = run_synthetic_smoke()
        print(json.dumps({"metrics": report.metrics, "claim_status": report.claim_status.value}, sort_keys=True))
        return 0
    if args.command == "smoke-real":
        from .real_smoke import run_real_data_smoke
        from .result_wiring import emit_real_smoke_results

        config = _real_config_from_args(args)
        report = run_real_data_smoke(config)
        payload = report.to_jsonable()
        payload["contract_results"] = _paths_payload(emit_real_smoke_results(report, config, results_dir=args.results_dir))
        print(json.dumps(payload, sort_keys=True))
        return 0
    if args.command == "gate2-parity":
        from .parity import build_parity_table
        from .real_smoke import run_real_data_smoke
        from .result_wiring import emit_gate2_results

        config = _real_config_from_args(args)
        report = run_real_data_smoke(config)
        if args.reference_path is not None:
            reference_rows = json.loads(args.reference_path.read_text(encoding="utf-8"))
        else:
            from .same_split_baseline import default_reference_context_path, run_same_split_baseline

            baseline = run_same_split_baseline(
                config=config,
                reference_context_path=default_reference_context_path(),
                grid_size=args.grid_size,
                baseline_ridge_alpha=args.baseline_ridge_alpha,
                knn_neighbors=args.knn_neighbors,
                local_smoke_metrics=report.metrics,
            )
            reference_rows = list(baseline["reference_rows"])
        table = build_parity_table(report.metrics, reference_rows, gate3_complete=True)
        payload = {"real_smoke": report.to_jsonable(), "parity": table}
        outputs = {}
        if args.out_path is not None:
            outputs["parity_table"] = args.out_path
        payload["contract_results"] = _paths_payload(
            emit_gate2_results(report, table, config, results_dir=args.results_dir, outputs=outputs)
        )
        _write_json_if_requested(args.out_path, payload)
        print(json.dumps(payload, sort_keys=True))
        return 0
    if args.command == "gate3-analysis":
        from .gate3 import run_gate3_analysis
        from .result_wiring import emit_gate3_results

        config = _real_config_from_args(args)
        payload = run_gate3_analysis(
            expression_path=config.expression_path,
            patches_path=config.patches_path,
            max_spots=config.max_spots,
            max_genes=config.max_genes,
            max_gene_scan=config.max_gene_scan,
            test_stride=config.test_stride,
            ridge_alpha=config.ridge_alpha,
            seed=config.seed,
        )
        outputs = {}
        if args.out_path is not None:
            outputs["gate3_table"] = args.out_path
        payload["contract_results"] = _paths_payload(
            emit_gate3_results(payload, config, results_dir=args.results_dir, outputs=outputs)
        )
        _write_json_if_requested(args.out_path, payload)
        print(json.dumps(payload, sort_keys=True))
        return 0
    if args.command == "claim-status":
        print(graduation_claim_status().value)
        return 0
    print(evaluate_claim_gate(ClaimGateEvidence()).value)
    return 0


def _real_config_from_args(args: argparse.Namespace) -> RealSmokeConfig:
    return RealSmokeConfig(
        expression_path=args.expression_path or default_expression_path(),
        patches_path=args.patches_path or default_patches_path(),
        max_spots=args.max_spots,
        max_genes=args.max_genes,
        max_gene_scan=args.max_gene_scan,
        test_stride=args.test_stride,
        ridge_alpha=args.ridge_alpha,
        seed=args.seed,
    )


def _write_json_if_requested(path: Path | None, payload: dict[str, object]) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _paths_payload(paths: dict[str, Path]) -> dict[str, str]:
    return {key: str(path) for key, path in paths.items()}


if __name__ == "__main__":
    raise SystemExit(main())
