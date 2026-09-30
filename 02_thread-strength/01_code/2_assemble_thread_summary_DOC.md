# `2_assemble_thread_summary.Rmd`

Joins the trace-level extraction to the hand-measured plaque areas, computes adhesion, and
writes a curation-ready candidate plus two worklists.

## Pipeline position

```
1_extract_tensometer_data.Rmd  ->  03_analyses/01_extract-tensometer-data/thread-summary-raw-output.xlsx
              +  02_data/pad_area_measurements.xlsx                            (hand-measured)
2_assemble_thread_summary.Rmd  ->  03_analyses/02_assemble-thread-summary/
                                     thread-summary-candidate.xlsx
                                     pad-area-worklist.xlsx
   (manual: drop bad runs against the QC plots)
                               ->  03_analyses/thread-summary.xlsx
3_analyze_thread_strength.Rmd  ->  03_analyses/03_analyze-thread-strength/
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

One row per extracted trace, sorted by folder then mussel then thread. A new trace is added
as a row with blank `pad_area`, and the worklist below finds it. The script still accepts
`mussel_ID` / `thread_num`, so an older copy of the file joins.

## The join key

`mussel` + `thread` + `thread_trt`. Verified unique on both sides. `thread_trt` is not
optional: for an animal pulled before and after exposure, `mussel` + `thread` alone
matches two different traces and would silently duplicate rows.

## Picture matching

Section 3 indexes `02_data/pictures/`, which is split `control` / `treatment`, matching the
pre/post split. A single image sometimes covers several threads (`T020_01-02-03.png`); those
are expanded so each thread gets its own entry. `picture_status` is `available` or `missing`
for every trace.

**This is what makes the worklist actionable**, because a trace with no image cannot be
measured from what is in the repository.

## Outputs

Written to `03_analyses/02_assemble-thread-summary/`.

### `thread-summary-candidate.xlsx`, sheet `data`

One row per trace, the shape script 3 expects, plus `pad_notes`, `picture_status`,
`source_folder` and `file` for provenance. Script 3 ignores extra columns, so it can be
saved as `thread-summary.xlsx` as-is once the bad runs are removed.

### `pad-area-worklist.xlsx`

- **`to_measure`**: one row per trace with no plaque measurement. Carries the path to the QC
  plot and to the microscope image, plus `picture_status`.
- **`pairing_gaps`**: one row per animal that ought to pair before/after but does not, with
  `gap_type`:
  - `pair_blocked_unmeasured` — traces exist on both sides, but a plaque measurement is
    missing on at least one side. **Fixable, if an image exists.**
  - `missing_pre_trace` — day-3 threads exist, no pre-exposure trace was recorded.
  - `missing_post_trace` — pre-exposure threads exist, no day-3 trace was recorded.

## Known data gaps

Some traces were measured from microscope images that are not in `02_data/pictures/`; they
are measured but carry `picture_status = missing`.

Some pairing gaps are trace gaps that no measurement can fix: the thirteen day-1 animals
pulled at baseline and never again by design, and T137 (control), which has no day-3
trace. The script prints the pairing status and the gaps by arm. The control arm's pairs
are what give the design a within-subject control contrast.
