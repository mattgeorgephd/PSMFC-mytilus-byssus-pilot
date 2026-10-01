# `03_assemble_thread_summary.Rmd`

Joins the trace-level extraction to the hand-measured plaque areas, computes adhesion, and
writes the thread summary that scripts 04 and 05 and the `09_gene-mechanics-correlation`
scripts read.

## Pipeline position

```
02_extract_tensometer_data.Rmd  ->  03_analyses/02_extract-tensometer-data/thread-summary-raw-output.xlsx
              +  02_data/pad_area_measurements.xlsx                            (hand-measured)
03_assemble_thread_summary.Rmd  ->  03_analyses/03_assemble-thread-summary/thread-summary.xlsx
04_analyze_thread_strength.Rmd  ->  03_analyses/04_analyze-thread-strength/
05_decompose_adhesion.Rmd       ->  03_analyses/05_decompose-adhesion/
```

## Where each fact comes from

One fact, one source. Nothing is copied into two files where the copies could drift.

| fact | grain | source |
|---|---|---|
| `max_force`, `integral`, `max_displacement` | thread | the trace, via script 02 |
| `thread_trt`, `phase`, `day` | thread | the tensometer subfolder, via script 02 |
| `species`, `mussel_trt`, `rna_sequenced` | mussel | `mussel-treatment-key.csv`, via script 02 |
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

`pad_area` and `failure` are measured by hand from microscope images of each plaque; they
cannot be derived from a force trace.

One row per extracted trace, sorted by folder then mussel then thread. A new trace is added
as a row with blank `pad_area`; the script reports how many traces have no row at all, and
scripts 04 and 05 exclude (and count) traces with no `pad_area`. The script still accepts
`mussel_ID` / `thread_num`, so an older copy of the file joins.

## The join key

`mussel` + `thread` + `thread_trt`. Verified unique on both sides. `thread_trt` is not
optional: for an animal pulled before and after exposure, `mussel` + `thread` alone
matches two different traces and would silently duplicate rows.

## Output

Written to `03_analyses/03_assemble-thread-summary/`.

### `thread-summary.xlsx`, sheet `data`

One row per trace: `species`, `sort_ID`, `mussel_ID`, `thread_num`, `note`, `phase`, `day`,
`mussel_trt`, `thread_trt`, `max_force`, `integral`, `max_displacement`, `pad_area`,
`adhesion_kpa`, `failure`, `rna_sequenced`, plus `pad_notes`, `source_folder` and `file` for
provenance. Scripts 04 and 05 ignore the extra columns.

Re-running the script overwrites this file, so any hand edit made to it is lost on the next
run. The committed copy is identical, cell for cell, to what the script writes from the
committed inputs.

## Known data gaps

Pairing is listed in script 02's `pairing` sheet: 13 animals were pulled at baseline only (by
design) and 13 at day 3 only; the 48 animals with both are the set the per-animal ANCOVA in
scripts 04 and 05 uses. The control arm's pairs are what give the design a within-subject
control contrast.
