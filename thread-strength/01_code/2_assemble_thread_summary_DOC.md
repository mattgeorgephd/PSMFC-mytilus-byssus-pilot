# `2_assemble_thread_summary.Rmd`

Joins the trace-level extraction to the hand-measured plaque areas, computes adhesion, and
writes a curation-ready candidate plus two worklists.

## Pipeline position

```
1_extract_tensometer_data.Rmd  ->  03_analyses/thread-summary-raw-output.xlsx   (every trace)
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

### What changed in `pad_area_measurements.xlsx`

The committed version carried `species`, `mussel_ID`, `thread_num`, `group`, `mussel_trt`,
`thread_trt`, `day`, `pad_area`, `failure` for 299 rows. Four of those columns were second
copies of facts owned elsewhere, and `group` was the retired column whose value `control`
collided with two other meanings. They were checked against the key before removal:
`mussel_trt` and `species` agreed on all 299 rows, and `day` is determined by the folder.

The file is now:

| column | meaning |
|---|---|
| `mussel` | canonical tag, prefix + zero-padded 3 digits |
| `thread` | integer |
| `thread_trt` | `lab_control`, `baseline`, `treatment_control`, `OA`, `OW`, `DO` |
| `pad_area` | plaque cross-sectional area, mm². **Blank where not yet measured** |
| `failure` | `cohesive`, `peeling`, `tearing`, `thread`. Blank where not yet measured |
| `notes` | free text, yours |

It now holds **all 375 rows**, one per extracted trace, sorted by folder then mussel then
thread. The 299 measured values are unchanged, verified value-for-value. The 76 outstanding
rows are present with blank `pad_area`, so filling them in is typing into existing cells
rather than adding rows. Filter on a blank `pad_area` to find them.

The script still accepts `mussel_ID` / `thread_num`, so an older copy of the file joins.

## The join key

`mussel` + `thread` + `thread_trt`. Verified unique on both sides. `thread_trt` is not
optional: 45 animals were pulled before and after exposure, so `mussel` + `thread` alone
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

375 rows, the shape script 3 expects, plus `pad_notes`, `picture_status`, `source_folder` and
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

**Plaque measurements: 299 of 375 traces.**

| folder | still needed | image available | image missing |
|---|---|---|---|
| `00_baseline` | 41 | 1 | 40 |
| `00_laboratory_control` | 5 | 1 | 4 |
| `01_treatment_control` | 3 | 2 | 1 |
| `02_OA_treatment` | 7 | 0 | 7 |
| `03_OW_treatment` | 5 | 1 | 4 |
| `04_DO_treatment` | 15 | 0 | 15 |
| **total** | **76** | **5** | **71** |

**Only 5 of the 76 can be measured from the images in this repository.** Every one of the 299
already-measured rows has a committed image, so the 71 are not a matching failure: those
images are not in `02_data/pictures/`. `pictures/control/` covers T01–T58 only, with no T07;
`pictures/treatment/` covers T14–T58, T89–T94 and T118–T135.

### Pairing

| arm | mussels | traces on both sides | usable adhesion pairs | blocked by measurement |
|---|---|---|---|---|
| DO | 20 | 12 | 9 | 3 |
| OA | 18 | 11 | 9 | 2 |
| OW | 24 | 12 | 12 | 0 |
| control | 12 | 10 | **0** | **10** |

Fifteen animals have traces on both sides but no usable pair, and in **every one of them the
missing side is the pre-exposure one**: T110, T111, T112 (DO); T118, T119 (OA); T126–T135
(control). None of their pre-exposure images is in the repository.

The control arm's ten pairs are the ones that matter most: they are the only route to a
within-subject control contrast, and their absence is what makes the mixed model rank
deficient in script 3.

Separately, 15 animals have pre-exposure traces but no day-3 trace (including T136 and T137,
the two day-3 control animals whose post pulls are missing), and 14 have day-3 traces but no
pre-exposure trace.
