# PathScope-ST — references (with code) & datasets

Consolidated reference + dataset index. Paper DOIs verified via Crossref and code
repositories via the GitHub API on 2026-06-09. See `BASELINE_REFERENCES.md` for the
full provenance and audit boundary.

## Reference papers & method baselines (with public code)

| Role | Method | Venue / year | DOI | Code |
|------|--------|--------------|-----|------|
| Primary | Path2Space — AI-predicted ST unlocks breast-cancer biomarkers from pathology | Cell 2026 | `10.1016/j.cell.2026.04.023` | Zenodo `10.5281/zenodo.14729337` |
| Primary (impl) | STimage — histology→ST prediction implementation baseline | — | — | https://github.com/BiomedicalMachineLearning/STimage |
| Baseline | HEST / HEST-1k — histology+ST benchmark infra | — | — | https://github.com/mahmoodlab/hest |
| Baseline | HisToSGE — high-res ST prediction from H&E | — | — | https://github.com/wenwenmin/HisToSGE |
| Baseline | Hist2ST — histology→ST deep baseline | — | — | https://github.com/biomed-AI/Hist2ST |
| Baseline | STAMP — spatially-aware pathology/ST pretraining | — | — | https://github.com/Hanminghao/STAMP |

## Datasets

Histology + spatial-transcriptomics paired data (user-supplied; no shipped catalog):
- HEST / HEST-1k benchmark (H&E ↔ ST pairs).
- Path2Space / STimage breast-cancer cohorts (Pearson / interval-coverage diagnostics).

> Verification: Path2Space DOI confirmed in Crossref (Cell, paywalled PDF); STimage / HEST /
> HisToSGE / Hist2ST / STAMP repos confirmed live via GitHub API (2026-06-09).
