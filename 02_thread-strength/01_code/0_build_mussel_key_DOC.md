# 0_build_mussel_key.Rmd

Builds `02_data/mussel-treatment-key.csv`, the repo-internal lookup from mussel tag to
experimental attributes: the mussel-level facts, in one file.

Run it whenever `01_mussel-measurements/mussel-size-measurements.xlsx` changes. It is
deterministic: unchanged inputs produce a byte-identical file.

## Input

`01_mussel-measurements/mussel-size-measurements.xlsx` at the repository root, sheet `data`
(107 rows). The script stops with a clear message if the workbook is missing or the sheet is
missing a column it needs, rather than writing a partial key.

## Output

`02_thread-strength/02_data/mussel-treatment-key.csv`, one row per *M. trossulus* tag.

| column | meaning |
|---|---|
| `mussel` | canonical tag, prefix + zero-padded 3 digits |
| `species` | |
| `mussel_trt` | the arm the animal was assigned to: `OA`, `OW`, `DO`, `control`, `lab_control` |
| `group` | cohort label carried through from the workbook, **not** a before/after axis |
| `days_in_trt` | days in the exposure system (`0` or `3`), kept as text so a non-numeric value cannot break the key |
| `rna_sequenced` | |

CSV rather than xlsx on purpose: the key is small, it is reviewed in `git diff`, and a binary
workbook is invisible there.

## Why `group` is carried but not used as a timepoint

`group` in the workbook is a **mussel cohort** label and reads as a timepoint when it is not:

- `group = control`, `days_in_trt = 0` (41 animals) — the twelve lab-reference animals,
  T001 to T012; thirteen animals pulled at baseline and never again; and sixteen with no
  tensometer trace
- `group = treatment`, `days_in_trt = 3` (66 animals) — everyone else, covering **both** their
  pre-exposure and their post-exposure threads

So the pre-exposure folder contains traces from animals in *both* groups. The before/after
axis comes from the folder (`phase` in script 1), never from this column.

## Validation the script performs

- Every mussel tag parses, and no tag is duplicated. Either failure is a hard stop.
- Every mussel with a tensometer trace has a key row; any that does not is named in a
  warning.
- The arm implied by the folder matches `mussel_trt`; any disagreement is tabulated. This is
  the check that would catch a trace filed in the wrong folder.
- A coverage table of who is in the key but has no traces (21 animals in the current
  workbook).

## Known gap

`pad_area`, `failure` and the technician-recorded `maximum_force` were **thread-level**
columns of the project Google Sheet. No `01_mussel-measurements` sheet contains them, so this key
cannot supply them; they live in `02_data/pad_area_measurements.xlsx` (one row per trace),
which `2_assemble_thread_summary.Rmd` joins to the extraction.
