# 0_build_mussel_key.Rmd

Builds `02_data/mussel-treatment-key.csv`, the repo-internal lookup from mussel tag to
experimental attributes: the mussel-level facts, in one file.

Run it whenever `morphometrics/02_data/morphometrics - tross.xlsx` changes. It is
deterministic: unchanged inputs produce a byte-identical file.

## Input

`morphometrics/02_data/morphometrics - tross.xlsx`, sheet `mussels-all` (119 rows).
Note the spaces and hyphen in the filename. The script stops with a clear message if the
sheet is missing a column it needs, rather than writing a partial key.

## Output

`thread-strength/02_data/mussel-treatment-key.csv`, one row per *M. trossulus* tag.

| column | meaning |
|---|---|
| `mussel` | canonical tag, prefix + zero-padded 3 digits |
| `species` | |
| `mussel_trt` | the arm the animal was assigned to: `OA`, `OW`, `DO`, `control`, `lab_control`, `desiccation` |
| `group` | cohort label carried through from morphometrics, **not** a before/after axis |
| `days_in_trt` | text, because values are not all numeric (the desiccation arm is `24h`) |
| `date_sampled` | ISO date, empty where unrecorded |
| `rna_sequenced` | |

CSV rather than xlsx on purpose: the key is small, it is reviewed in `git diff`, and a binary
workbook is invisible there.

## Why `group` is carried but not used as a timepoint

`group` in morphometrics is a **mussel cohort** label and reads as a timepoint when it is not:

- `group = control`, `days_in_trt = 0` — the twelve lab-reference animals, T001 to T012
- `group = control`, `days_in_trt = 1` — thirteen animals pulled at baseline and never again
- `group = treatment`, `days_in_trt = 3` — everyone else, covering **both** their
  pre-exposure and their post-exposure threads

So the pre-exposure folder contains traces from animals in *both* groups. The before/after
axis comes from the folder (`phase` in script 1), never from this column.

## Validation the script performs

- Every mussel tag parses, and no tag is duplicated. Either failure is a hard stop.
- Every mussel with a tensometer trace has a key row. Currently **86 of 86**.
- The arm implied by the folder matches `mussel_trt`. Currently **380 of 380 traces agree**,
  with no exceptions. This is the check that would catch a trace filed in the wrong folder.
- A coverage table of who is in the key but has no traces. Currently 33 animals: the twelve
  `desiccation` animals, five day-3 animals and sixteen day-1 animals.

## Known gap

`pad_area`, `failure` and the technician-recorded `maximum_force` were **thread-level**
columns of the project Google Sheet. No morphometrics sheet contains them, so this key
cannot supply them; they live in `02_data/pad_area_measurements.xlsx` (one row per trace,
all 380 measured), which `2_assemble_thread_summary.Rmd` joins to the extraction.
