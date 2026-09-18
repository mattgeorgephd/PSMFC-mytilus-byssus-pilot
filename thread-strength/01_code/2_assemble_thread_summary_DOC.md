# `2_assemble_thread_summary.Rmd`

Joins the trace-level extraction to the hand-measured plaque areas, computes adhesion, and
writes a curation-ready candidate plus two worklists.

## Pipeline position

```
1_extract_tensometer_data.Rmd  ->  03_analyses/extract-tensometer-data/thread-summary-raw-output.xlsx
              +  02_data/pad_area_measurements.xlsx                            (hand-measured)
2_assemble_thread_summary.Rmd  ->  03_analyses/assemble-thread-summary/
                                     thread-summary-candidate.xlsx
                                     pad-area-worklist.xlsx
   (manual: drop bad runs against the QC plots)
                               ->  03_analyses/thread-summary.xlsx
3_analyze_thread_strength.Rmd  ->  03_analyses/analyze-thread-strength/
```

## Where each fact comes from

One fact, one source. Nothing is copied into two files where the copies could drift.

| fact | grain | source |
|---|---|---|
| `max_force`, `integral`, `max_displacement` | thread | the trace, via script 1 |
| `thread_trt`, `phase`, `day` | thread | the tensometer subfolder, via script 1 |
| `species`, `mussel_trt`, `rna_sequenced` | mussel | `mussel-treatment-key.csv`, via script 1 |
| `pad_area`, `failure` | thread | `02_data/pad_area_measurements.xlsx` |

### `pad_area_measurements.xlsx`

| column | meaning |
|---|---|
| `mussel` | canonical tag, prefix + zero-padded 3 digits |
| `thread` | integer |
| `thread_trt` | `lab_control`, `baseline`, `treatment_control`, `OA`, `OW`, `DO` |
| `pad_area` | plaque cross-sectional area, mm². Blank where not yet measured |
| `failure` | `cohesive`, `peeling`, `tearing`, `thread`. Blank where not yet measured |
| `notes` | free text, yours |

One row per extracted trace (380), sorted by folder then mussel then thread; every row is
measured. A new trace is added as a row with blank `pad_area`, and the worklist below finds
it. The script still accepts `mussel_ID` / `thread_num`, so an older copy of the file joins.

## The join key

`mussel` + `thread` + `thread_trt`. Verified unique on both sides. `thread_trt` is not
optional: 47 animals were pulled before and after exposure, so `mussel` + `thread` alone
matches two different traces and would silently duplicate rows.

## Picture matching

Section 3 indexes `02_data/pictures/`, which is split `control` / `treatment`, matching the
pre/post split. A single image sometimes covers several threads (`T020_01-02-03.png`); those
are expanded so each thread gets its own entry. `picture_status` is `available` or `missing`
for every trace.

**This is what makes the worklist actionable**, because a trace with no image cannot be
measured from what is in the repository.

## Outputs

Written to `03_analyses/assemble-thread-summary/`.

### `thread-summary-candidate.xlsx`, sheet `data`

380 rows, the shape script 3 expects, plus `pad_notes`, `picture_status`, `source_folder` and
`file` for provenance. Script 3 ignores extra columns, so it can be saved as
`thread-summary.xlsx` as-is once the bad runs are removed.

### `pad-area-worklist.xlsx`

- **`to_measure`**: one row per trace with no plaque measurement. Carries the path to the QC
  plot and to the microscope image, plus `picture_status`.
- **`pairing_gaps`**: one row per animal that ought to pair before/after but does not, with
  `gap_type`:
  - `pair_blocked_unmeasured` — traces exist on both sides, but a plaque measurement is
    missing on at least one side. **Fixable, if an image exists.**
  - `missing_pre_trace` — day-3 threads exist, no pre-exposure trace was recorded.
  - `missing_post_trace` — pre-exposure threads exist, no day-3 trace was recorded.

## Current state

**Plaque measurements: 380 of 380 traces.** `to_measure` is empty, and the candidate table
equals `03_analyses/thread-summary.xlsx` value for value (force, area, failure mode,
adhesion) for all 380 rows. 293 traces have a committed microscope image; the other 87 were
measured from images that are not in `02_data/pictures/`.

### Pairing

| arm | animals with traces | traces on both sides | pre only | post only |
|---|---|---|---|---|
| control | 12 | 11 | 1 (T137) | 0 |
| OA | 18 | 12 | 5 | 1 |
| OW | 24 | 12 | 2 | 10 |
| DO | 20 | 12 | 6 | 2 |

The 27 unpaired animals are all trace gaps (`missing_pre_trace`, `missing_post_trace`),
not measurement gaps: 14 have a pre-exposure pull only (13 day-1 animals never re-pulled
by design, plus T137) and 13 a day-3 pull only. The control arm's eleven pairs are what
give the design a within-subject control contrast.
