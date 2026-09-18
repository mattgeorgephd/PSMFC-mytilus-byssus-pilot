# 1_extract_tensometer_data.Rmd

Reads every raw tensometer trace under `02_data/tensometer_output/`, joins the mussel key,
and writes one trace-level table to `03_analyses/extract-tensometer-data/thread-summary-raw-output.xlsx`.

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
followed by three tab-separated rows, one per channel. All 380 current files are in this
form (4 lines, 70 to 773 samples per channel).

The script reads the lines directly rather than `read.delim()` + transpose, so nothing
depends on how many lines R scans to infer a column count. A file that is not in the wide
format raises a clear error instead of being silently misread.

### Filename grammar

`<prefix><digits>_<digits><optional note>`, with or without a `.txt` extension. The trailing
text is the technician's note (`try2`, `fluke`, `thirdTry`). It is kept as its own `note`
column so every run is imported and can be reviewed against its QC plot. Files that do not
match are skipped **with a warning**, never silently.

> `mussel` + `thread` is **not** a unique key. 47 animals were pulled both before and after
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
161 traces carry `thread_trt = "baseline"` while `mussel_trt` is OA, OW, DO or control.
Keeping them in separate columns from separate sources is what stops them being conflated.

`phase` has **three** levels, not two. The lab-reference animals (day 0, never in the
experimental system) are not the same thing as the paired pre-exposure baselines. A binary before/after would pool them into one "pre" group and
quietly contaminate any paired contrast.

### Folder map

```
00_laboratory_control   lab_control        lab    day 0   expect mussel_trt = lab_control
00_baseline             baseline           pre    day 1   expect mussel_trt = (any arm)
01_treatment_control    treatment_control  post   day 3   expect mussel_trt = control
02_OA_treatment         OA                 post   day 3   expect mussel_trt = OA
03_OW_treatment         OW                 post   day 3   expect mussel_trt = OW
04_DO_treatment         DO                 post   day 3   expect mussel_trt = DO
```

The table lives in the `folder_labels` chunk and is the script's **only** point of
interpretation. It also accepts the phase-first spellings (`00_lab_reference`,
`01_pre_exposure`, `02_post_control`, `03_post_OA`, `04_post_OW`, `05_post_DO`), so the
folders can be renamed without touching the script; only folders that exist are used.

An unmapped folder is still processed, with labels guessed from its name, `day = NA`, and a
warning. It is never silently mislabelled.

### One `phase` column, not `group`

`phase` (`lab` / `pre` / `post`) is the before/after axis. There is no `group` column: it
would be collinear with `phase` and `thread_trt`, and a value like `control` would collide
with `mussel_trt == "control"` (the control tank arm) and `thread_trt == "treatment_control"`.

---

## 3. Output

`03_analyses/extract-tensometer-data/thread-summary-raw-output.xlsx`, three sheets. `data` is sheet 1, so
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

## 4. How the two derived quantities are computed

### The integral is a trapezoid on the raw trace

`sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)` on the raw force/time trace. `loess()` is
fitted for the QC plot only, wrapped in `tryCatch` so a fit failure degrades to a plot
without a smoother rather than killing the run.

### Missing samples are handled by position, not blanket-zeroed

A literal `NaN` as the **first** Force sample is normal instrument behaviour at the start of
a pull and is set to zero (47 of the 380 files). A missing sample away from index 1 is also
zeroed, but **warns with the filename and index**, and the per-trace count survives into the
output as `n_na_force`. Time is checked for monotonicity before integration and sorted, with
a warning, if it is not.

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

## 6. Data state, as of the 380-trace extraction

- 380 traces from 86 animals: 38 lab-reference (12 animals), 161 pre-exposure (61), 34
  day-3 control (11), 39 OA (13), 67 OW (22), 41 DO (14). Every trace has a plaque
  measurement and a failure mode in `thread-summary.xlsx`, and the curated `max_force`
  values equal a fresh extraction for all 380.
- 47 animals were pulled at both timepoints (control 11, OA 12, OW 12, DO 12); 14 have a
  pre-exposure pull only and 13 a day-3 pull only (`pairing` sheet;
  `pad-area-worklist.xlsx`, sheet `pairing_gaps`).
- `01_treatment_control` holds eleven of the twelve day-3 control animals listed in the
  morphometrics (T126 to T136); T137 has no day-3 trace.
- The `desiccation` arm (12 animals, 24 h) has no tensometer traces; 21 other animals in
  the mussel key (16 day-1, 5 day-3) have none either.
