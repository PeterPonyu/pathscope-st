# PathScope-ST baseline references

Verification date: 2026-06-07

## Baseline decision summary

| Role | Baseline | Decision |
|---|---|---|
| Primary | Path2Space codebase + STimage implementation baseline | Use as the first open-code/public-artifact starting point because it directly matches this track's input-output problem. |
| Secondary | See table below | Use only for comparison, adapter design, and ablation inspiration; do not copy implementation. |

## Primary baseline

- Paper title: AI-predicted spatial transcriptomics unlocks breast cancer biomarkers from pathology
- Venue/date: Cell, online 2026-05-08, DOI 10.1016/j.cell.2026.04.023
- Article URL: https://www.sciencedirect.com/science/article/pii/S0092867426004587
- Code/artifact URL: https://zenodo.org/records/14729337 and https://github.com/BiomedicalMachineLearning/STimage
- Verification date: 2026-06-07
- Default branch/artifact: Zenodo artifact; STimage master
- Observed HEAD SHA or DOI: Path2Space DOI 10.5281/zenodo.14729337; STimage 5a9696c3568a4f2728c96ed978bd60780eb437ec
- Local audit checkout/artifact: `baselines/Path2Space-codebase-original; baselines/STimage-original`
- License note: Path2Space Zenodo CC-BY-NC-4.0; STimage BSD-3-Clause license file observed locally
- Local use: Path2Space anchors the latest clinical virtual-ST paper; STimage provides a public GitHub implementation baseline for robust/interpretable gene and cell-type prediction from histology.
- Fallback: If this public code/artifact becomes unavailable, mark this track `deferred-unverified` until a comparable open-code baseline is found.
- Verification command/evidence:
  - `git ls-remote --symref <repo> HEAD` for GitHub/Git/Bioconductor repositories.
  - `zenodo.org/api/records/<record>` plus local file checks for Zenodo artifacts.

## Secondary verified references

| Baseline | Role | Code/artifact URL | Local audit checkout | Branch/artifact | Observed SHA/DOI | License note |
|---|---|---|---|---|---:|---|
| HEST | Histology + spatial transcriptomics dataset/benchmark infrastructure | https://github.com/mahmoodlab/hest | `baselines/HEST-original` | `main` | `3ddb5eaf5bd2` | CC BY-NC-SA 4.0 license file observed locally |
| HisToSGE | High-resolution ST prediction from histology images | https://github.com/wenwenmin/HisToSGE | `not cloned yet` | `main` | `d68d11740aa4` | license requires re-check |
| Hist2ST | Histology-to-ST deep baseline | https://github.com/biomed-AI/Hist2ST | `not cloned yet` | `main` | `7480e5d1f771` | license requires re-check |
| STAMP | Spatially-aware pathology/ST pretraining reference | https://github.com/Hanminghao/STAMP | `not cloned yet` | `main` | `78c4fd2f4e77` | license requires re-check |

## Brand independence note

Reference names in this file are provenance labels only. Local package names, CLI commands, figure labels, and manuscript novelty claims must use `PathScope-ST` terminology and the independent refinements in `README.md`, not upstream branding.

## Gate-2 baseline_comparison provenance (2026-06-08)

| Baseline | Gate-2 status | Evidence |
|---|---|---|
| Path2Space family reference | `REFERENCE_REPORTED` | Local frozen Zenodo audit copy `../baselines/Path2Space-codebase-original @ fefcc0a` lacks the demo `input_data/` and `output_data/` directories needed to replay predictions on the COAD card. Reference-only range recorded from AACR Cancer Immunology Research 2026 abstract C060: external validation Pearson `>0.4` for 546 genes and selected biomarkers `>0.8`. |
| STimage | `REFERENCE_REPORTED` | Local frozen clone `../baselines/STimage-original @ 5a9696c`; isolated editable-install dry-run blocked by package metadata/direct dependency handling and no predictions exist for the COAD held-out split. Reference-only values recorded from Nature Communications 2026: top-100 gene PCC range `0.2–0.8`, TCGA bulk-matched average PCC `0.48`, drug-response cohort average PCC `0.47`. |

Gate-2 parity artifact: `experiments/gate2_baseline_comparison/parity_table.json`.
