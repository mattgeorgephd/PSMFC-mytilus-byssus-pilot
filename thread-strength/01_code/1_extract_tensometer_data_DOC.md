# 1_extract_tensometer_data.Rmd

Reads every raw tensometer trace under `02_data/tensometer_output/`, joins the mussel key,
and writes one trace-level table to `03_analyses/thread-summary-raw-output.xlsx`.

Run it from inside `thread-strength.Rproj`. Run `0_build_mussel_key.Rmd` first if the
morphometrics workbook has changed.

---

## 1. Inputs

| Path | Role |
|---|---|
| `02_data/tensometer_output/<folder>/<mussel>_<thread>[note][.txt]` | raw traces |
| `02_data/mussel-treatment-key.csv` | mussel-level metadata (optional; warns if absent) |

### Trace file format

The instrument writes each trace **wide**: a header line `Time  Displacement  Force`
followed by three tab-separated rows, one per channel. All 375 current files are in this
form (4 lines, 70 to 773 samples per channel).

The script reads the lines directly rather than `read.delim()` + transpose. The old approach
made `read.delim()` infer a column count from a row that can hold 773 fields, which works
only because R happens to scan the first five lines; reading lines is both faster and has no
such dependency. A file that is not in the wide format now raises a clear error instead of
being silently misread.

### Filename grammar

`<prefix><digits>_<digits><optional note>`, with or without a `.txt` extension. The trailing
text is the technician's note (`try2`, `fluke`, `thirdTry`). It is kept as its own `note`
column so every run is imported and can be reviewed against its QC plot. Files that do not
match are skipped **with a warning**, never silently.

> `mussel` + `thread` is **not** a unique key. 45 animals were pulled both before and after
> exposure, so the same pair appears in two folders. The unique key is
> `mussel` + `thread` + `thread_trt` (equivalently, + `source_folder`).

---

## 2. The label model

Three facts are kept in separate columns because **a mussel and its threads did not
necessarily experience the same thing**.

| column | grain | source | meaning |
|---|---|---|---|
| `thread_trt` | thread | folder name | what the **thread** was built in |
| `phase` | thread | folder name | `lab` / `pre` / `post` |
| `mussel_trt` | mussel | `mussel-treatment-key.csv` | the arm the **animal** was assigned to |

`mussel_trt` is a **destiny** label, not an experienced-condition label. For a `pre` thread it
describes the animal's future, not its past: *a baseline thread from an OA animal is not an
OA thread.* The two grains diverge in exactly one place, the pre-exposure folder, where all
158 traces carry `thread_trt = "baseline"` while `mussel_trt` is OA, OW, DO or control.
Keeping them in separate columns from separate sources is what stops them being conflated.

`phase` has **three** levels, not two. The lab-reference animals (day 0, never in the
experimental system) are not the same thing as the paired pre-exposure baselines. A binary before/after would pool them into one "pre" group and
quietly contaminate any paired contrast.

### Folder map

```
00_lab_reference   lab_control        lab    day 0   expect mussel_trt = lab_control
01_pre_exposure    baseline           pre    day 1   expect mussel_trt = (any arm)
02_post_control    treatment_control  post   day 3   expect mussel_trt = control
03_post_OA         OA                 post   day 3   expect mussel_trt = OA
04_post_OW         OW                 post   day 3   expect mussel_trt = OW
05_post_DO         DO                 post   day 3   expect mussel_trt = DO
```

The table lives in the `folder_labels` chunk and is the script's **only** point of
interpretation. It also carries the legacy folder names (`00_baseline`, `02_OA_treatment`,
...) so the script works before and after `rename_tensometer_folders.sh` is applied. Old and
new names never coexist, so the extra rows are inert; delete them once the rename is pushed.

An unmapped folder is still processed, with labels guessed from its name, `day = NA`, and a
warning. It is never silently mislabelled.

### `group` is deliberately not reproduced

In `thread-summary.xlsx`, `group`, `day` and `thread_trt` were perfectly collinear: every
crosstab was a clean diagonal, so three columns carried one column's worth of information.
Worse, `group == "control"` meant *day 0 or 1*, while `mussel_trt == "control"` means *the
control tank arm* and `thread_trt == "treatment_control"` is a third thing. `phase` replaces
`group` and does not collide with anything.

---

## 3. Output

`03_analyses/thread-summary-raw-output.xlsx`, three sheets. `data` is sheet 1, so
`read_excel(path)` with no `sheet` argument gets the trace table.

### Sheet `data`, one row per trace

| column | unit | notes |
|---|---|---|
| `source_folder` | | folder the trace came from |
| `thread_trt` | | what the thread was built in |
| `phase` | | `lab` / `pre` / `post` |
| `day` | | experiment day the thread was pulled |
| `mussel` | | canonical tag, prefix + 3 digits |
| `thread` | | integer |
| `note` | | technician note, `""` if none |
| `species`, `mussel_trt`, `rna_sequenced` | | joined from the mussel key |
| `max_force` | N | peak of the force channel |
| `integral` | N·s | area under the force/time curve |
| `max_displacement` | mm | peak extension |
| `n_points` | | samples in the trace |
| `duration_s` | s | |
| `n_na_force` | | missing force samples that were zero-filled |
| `file` | | source filename; also the QC plot name |

### Sheets `coverage` and `pairing`

`coverage` gives traces, mussels, threads per mussel, retry traces and zero-filled samples
per folder. `pairing` gives, per arm, how many animals have both a `pre` and a `post` pull.
Unpaired animals are not an error, but an unbalanced arm is what makes a per-arm marginal
mean from a random-intercept model unreliable, so the count belongs next to the data.

### QC plots

`03_analyses/extract-tensometer-data/QC_plots/<source_folder>/<file>.jpg`, one per trace:
raw force against time in blue with a loess smoother in red. Named after the source file, so
retries and flukes that share a mussel and thread never overwrite each other. The graphics
device is closed in a `finally` block, so a plotting failure cannot leave it open and corrupt
every later plot.

---

## 4. Two computational changes from the previous version

### The integral is now a real trapezoid

The old code was:

```r
current_loess <- loess(force ~ time, data = current_df)
auc <- sum(diff(current_loess$x) *
           (approx(current_loess$x, current_loess$y, n = length(current_loess$x)))$y[-1])
```

`loess()$y` is the **response**, not `$fitted`. So despite the comments, this never
integrated the smoothed curve; it integrated the raw trace, and it did so on a grid from
`approx(n = length(x))` (evenly spaced between `min(x)` and `max(x)`) that does not line up
with the `diff(x)` it was multiplied by.

It is now `sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)` on the raw trace: the same
intended quantity, computed correctly. Measured across all 375 traces, the two agree to a
**median of 0.006%** and a **maximum of 2.1%**. `loess()` is still fitted, for the QC plot
only, wrapped in `tryCatch` so a fit failure degrades to a plot without a smoother rather
than killing the run.

### Missing samples are handled by position, not blanket-zeroed

46 of the 375 files carry a literal `NaN` as the **first** Force sample. This is normal
instrument behaviour at the start of a pull, and it is set to zero.

The old code applied `is.na(x) <- 0` to every channel unconditionally, which would also zero
a mid-trace dropout without a word. Now a missing sample away from index 1 is still zeroed,
but **warns with the filename and index**, and the per-trace count survives into the output
as `n_na_force`. Time is also checked for monotonicity before integration and sorted, with a
warning, if it is not.

---

## 5. Integrity checks the script performs

1. **Stray files.** Anything in a trace folder not matching the filename grammar is reported.
2. **Unmapped folders.** A new subfolder is processed but flagged.
3. **Orphan mussels.** A tag with traces but no row in the mussel key is named in a warning.
4. **Folder versus morphometrics.** For a `post` folder, `thread_trt` and `mussel_trt` must
   agree. `expect_mussel_trt` in `folder_labels` encodes that, and any disagreement is
   tabulated. This catches a trace filed under the wrong folder, which no downstream script
   can detect on its own. It is `NA` for `pre`, where the two legitimately differ.
5. **Outlier screen.** Peak force by thread treatment, labelled, plus a Rosner test pooled
   and per arm (skipped where n < 25, which is the smallest sample Rosner is meant for).
   Flags candidates for curation; removes nothing.

---

## 6. Known data issues, as of the 375-trace extraction

- **76 traces have no plaque measurement** in `thread-summary.xlsx`, so `adhesion_kpa`
  cannot be computed for them: 41 pre-exposure, 15 DO, 7 OA, 5 lab, 5 OW, 3 post-control.
- **Three curated `max_force` values disagree with a fresh extraction**, and each propagated
  into `adhesion_kpa`:

  | trace | curated | extracted | integral agrees? | reading |
  |---|---|---|---|---|
  | `T125_03` (DO) | 0.026 | 0.265 | yes | transposed digits in the curated file |
  | `T029_02` (OW) | 0.105 | 0.130 | yes | curated peak looks hand-adjusted |
  | `T090_01` (OW) | 0.255 | 0.355 | **no**, 0.937 vs 0.381 | a different trace entirely |

- **`02_post_control` is missing post traces for T136 and T137.** Morphometrics lists twelve
  day-3 control animals (T126 to T137); the folder holds ten.
- **The `desiccation` arm has no tensometer traces at all** (12 animals, 24 h), nor do five
  day-3 DO animals.
- 15 animals have a pre-exposure pull but no post-exposure pull. 13 are day-1 animals that
  were never re-pulled by design; T136 and T137 are the two above.
