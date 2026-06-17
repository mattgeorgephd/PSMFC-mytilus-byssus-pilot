# respirometry

Standard metabolic rate (SMR) of mussels measured by closed-system respirometry,
before (baseline, T0) and after a 3-day stress exposure (T3). Re-runnable from this
repository, with one manual curation step (see below).

## Layout

```
respirometry/
├── respirometry.Rproj
├── 01_code/
│   ├── metabolic_rate_analysis_baseline.R   T0 timepoint
│   └── metabolic_rate_analysis_T3.R         T3 timepoint (identical logic, different date)
├── 02_data/
│   ├── measurements/{20211026_T0,20211029_T3}/   raw O2 traces, one .xlsx per run
│   └── constraints/                              hand-edited run constraints (see below)
└── 03_analyses/
    └── metabolic-rate/
        ├── QC_plots/{20211026_T0,20211029_T3}/   per-run oxygen-vs-time loess QC
        └── output/                               SMR slope tables
```

Both timepoint scripts write into the single `metabolic-rate/` output area, keyed by
date, because they are one procedure run at two timepoints and share combined outputs
(`SMR_slope_final.xlsx`, `processed_summary.xlsx`).

## How to run (two-phase, because of constraints)

1. Open `respirometry.Rproj` so `here::here()` anchors here.
2. Run a script down through the QC-plot loop. It reads `02_data/measurements/<date>/`,
   writes a constraints template to `02_data/constraints/constraints_<date>.xlsx` (only
   if one does not already exist), and writes QC plots to `03_analyses/metabolic-rate/QC_plots/<date>/`.
3. **Manual step:** inspect the QC plots, then edit `constraints_<date>.xlsx`: set the
   keep/recalc flag (column 3) and the start/stop time window (columns 4 and 5) for each run.
4. Run the rest of the script. It reads the edited constraints and computes the SMR slope
   into `03_analyses/metabolic-rate/output/`.

## Notes

- **Non-destructive constraints write.** The original script overwrote the constraints
  file on every run, which would wipe your manual start/stop edits on a re-run. The
  rewrite only writes the template when the file is absent. Delete the file to regenerate.
- **Leftover output name.** The baseline script writes `SMR_slope_heat_only_T5.xlsx`,
  a name carried over from other code. Consider renaming it to match the timepoint
  (e.g. `SMR_slope_20211026_T0.xlsx`); the committed output already uses that name.
