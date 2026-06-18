# `2_assemble_thread_summary.Rmd` (assemble the initial thread summary)

## What it does

Joins the tensometer extraction outputs (script 1) with per-mussel metadata to produce the initial,
pre-curation thread summary, split by species. It excludes nothing: every extracted trace is carried
through (including retries and flukes, flagged by the `note` column) so you can review each against
its QC plot and remove bad runs by hand. The curated trossulus table is what you save as
`02_data/thread-summary-clean.xlsx` for the analysis script.

## Pipeline position and the renumber

Extract (1) -> **assemble (2, new)** -> analyze (3). The companion `apply_assemble_setup.sh` renames
`2_analyze_thread_strength.Rmd` to `3_analyze_thread_strength.Rmd` via `git mv` (history preserved).
Note: the downstream input patch we discussed earlier now targets `3_analyze_thread_strength.Rmd`,
and since you renamed `summary.xlsx` to `thread-summary-OG.xlsx`, the analyze script will read
`thread-summary-clean.xlsx` once that patch is applied (after you create the clean file).

## Sources and join keys

| output column | source | key |
|---|---|---|
| max_force, integral, note | extraction outputs | the trace itself |
| group (before/after) | which folder the trace came from | n/a |
| species, mussel_trt | morphometrics `key` tab | mussel tag |
| treatment | "control" if before, else mussel_trt | derived |
| pad_area, failure, tech_max_force, tech_notes | tech `threads.xlsx` | mussel tag + thread + before/after |
| adhesion_kpa | max_force / pad_area * 1000 | derived |
| force_diff | max_force - tech_max_force | derived |
| rna_sequenced | presence in `rna_isolation_info` | mussel tag |

All tags are normalized to a 3-digit form (T37 -> T037) so the 2-digit control folder, the 3-digit
treatment folder, and the mixed-format metadata sheets join on the same key.

## Validation against `thread-summary-OG.xlsx`

The pipeline was checked row by row against your hand-built OG:

- Every trace receives a species and treatment; nothing is unmatched.
- The RNA flag reproduces the OG exactly (45/45 sequenced, 13/13 non-sequenced).
- Treatment matches the OG everywhere except T126-T135, where the pipeline outputs `control` and the
  OG had `DO`. The master `key` tab and the DE files both say control, so the OG label is an error the
  regeneration corrects.
- The import fix recovers 68 traces the OG never had (the previously undetected extension-less files).
- Computed peak force matches the technician-recorded force with a median difference of 0.0 N.

## Known issues to watch during curation

These are real and are surfaced by the summary, not hidden:

- **Pad area, 14 rows.** The tech sheet differs from the OG by a clean integer (1, 2, or 3 mm^2) for
  14 threads, and several tech values are implausibly small (0.53, 1.06 mm^2). The OG values look more
  reasonable. This suggests the tech sheet was edited after the OG, or had a calibration slip. Since
  pad_area drives `adhesion_kpa`, verify these before trusting their adhesion values. Affected:
  T014 t3, T029 t3, T041 t1, T046 t1-3, T094 t2-3, T125 t3, T130 t2-3, T133 t2, T134 t5, T135 t2.
- **Force, 5 rows.** `force_diff` will be large for T031 t2, T037 t1, T042 t2, T090 t1, T134 t5. For
  T090 and T134 the OG used the retry/fluke value, not the clean run; for the others the OG value
  matches neither trace and was likely hand-entered. Decide per row which value to keep.
- **Retry traces, 19 rows.** Clean and retry runs of the same thread appear as separate rows
  distinguished by `note`. Keep one per thread during curation.
- **Tech retest duplicates.** A few mussel/thread/group keys (e.g. T029) have two tech rows; the
  script keeps the first and reports the rest.
- **Missing values.** 1 trace has no pad area (no adhesion), and 39 traces have no technician force
  (the technician wrote a note instead of a clean value, which is itself the bad-run signal).
- **Control-mussel after-threads.** The summary now includes after-timepoint threads from control
  mussels (group treatment, treatment control), which the OG lacked because of the import bug. Decide
  whether the analysis pools these into "control" or drops them.

## Outputs

Written to `03_analyses/assemble-thread-summary/`:

- `thread-summary-trossulus.xlsx` (sheet "data"), 307 rows, the manuscript species.
- `thread-summary-gallo.xlsx` (sheet "data"), 104 rows, the capstone species.
- `force_computed_vs_technician.jpg`, the computed-vs-recorded force scatter.

Gallo before/after pairing is left for the capstone, since gallo tags genuinely diverge between
folders (the old retag situation); it does not affect the per-thread trossulus summary.

## Workflow

Run script 1 (extraction), then this script. Review `thread-summary-trossulus.xlsx` alongside the QC
plots, using `note`, `failure`, and `force_diff` to find bad runs. Remove what should be excluded,
resolve the flagged pad-area and force rows, and save the result as `02_data/thread-summary-clean.xlsx`
(sheet "data"). Then apply the downstream patch so script 3 reads it.
