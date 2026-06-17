# 02_data

Inputs for `../01_code/`. Outputs belong in `../03_analyses/`.

| Item | Description | Used by |
|------|-------------|---------|
| `measurements/20211026_T0/` | Raw dissolved-O2 traces (one `.xlsx` per run), baseline | both scripts (T0) |
| `measurements/20211029_T3/` | Raw dissolved-O2 traces, post-stress | both scripts (T3) |
| `constraints/constraints_<date>.xlsx` | Per-run keep/recalc flag and start/stop window, hand-edited after QC | both scripts |

`constraints/` is a curated intermediate: the script seeds it as a template, you edit
the time windows after reviewing the QC plots, and the script reads it back. It is not
raw data, but it lives here because the analysis consumes it as an input.
