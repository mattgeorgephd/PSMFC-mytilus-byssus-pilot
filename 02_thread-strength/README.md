# 02_thread-strength

Plaque adhesion strength and thread mechanics from tensometer pull tests, before and after a
3-day stress exposure (control, OA, OW, DO), plus a lab-reference group.

Self-contained and reproducible from this repository apart from one manual input, the
plaque measurements, noted below. It reads one file from outside this folder
(`../01_mussel-measurements/mussel-size-measurements.xlsx`, script 0) and sources the shared
`../tools/pipeline_checks.R` (scripts 3 and 4, for `RUN_provenance.txt`).

## Layout

```
02_thread-strength/
├── thread-strength.Rproj              open this first; it anchors here::here()
├── 01_code/
│   ├── 0_build_mussel_key.Rmd         mussel workbook -> 02_data/mussel-treatment-key.csv
│   ├── 1_extract_tensometer_data.Rmd  raw traces     -> 03_analyses/01_extract-tensometer-data/
│   ├── 2_assemble_thread_summary.Rmd  raw output + plaque areas -> thread summary
│   ├── 3_analyze_thread_strength.Rmd  thread summary -> adhesion plots + ANCOVA
│   ├── 4_decompose_adhesion.Rmd       thread summary -> force / area / extension ANCOVA
│   ├── _ancova.R                      the per-animal ANCOVA both scripts source
│   └── *_DOC.md                       companion documentation, one per script
├── 02_data/
│   ├── tensometer_output/<phase folders>/   raw force/displacement .txt traces
│   ├── pad_area_measurements.xlsx           hand-measured plaque area + failure mode, per trace
│   └── mussel-treatment-key.csv             mussel tag -> arm, species, rna flag (script 0)
└── 03_analyses/
    ├── 01_extract-tensometer-data/          script 1: thread-summary-raw-output.xlsx
    │   └── QC_plots/<source_folder>/        per-trace loess QC jpgs
    ├── 02_assemble-thread-summary/          script 2: thread-summary.xlsx, input to 3, 4 and 09
    ├── 03_analyze-thread-strength/          script 3: figures + STATS_ancova_*.csv
    └── 04_decompose-adhesion/               script 4: force / area / extension ANCOVA
```

## Tensometer folder layout

```
02_data/tensometer_output/
├── 00_laboratory_control/   day 0  T001-T012, never entered the experimental system
├── 00_baseline/             day 1  shared holding system, threads built before exposure
├── 01_treatment_control/    day 3  common-garden control tank
├── 02_OA_treatment/         day 3
├── 03_OW_treatment/         day 3
└── 04_DO_treatment/         day 3
```

Script 1 maps each folder to `thread_trt`, `phase` and `day` through its `folder_labels`
table, which also accepts the phase-first spellings (`00_lab_reference`, `01_pre_exposure`,
`02_post_control`, `03_post_OA`, `04_post_OW`, `05_post_DO`) should the folders ever be
renamed. Only folders that exist are used.

## The label model

A mussel and its threads did not necessarily experience the same thing. Three columns keep
that straight:

| column | grain | source | meaning |
|---|---|---|---|
| `thread_trt` | thread | folder name | what the **thread** was built in |
| `phase` | thread | folder name | `lab` / `pre` / `post` |
| `mussel_trt` | mussel | `mussel-treatment-key.csv` | the arm the **animal** was assigned to |

`mussel_trt` is a **destiny** label. For a `pre` thread it describes the animal's future, not
its past: a baseline thread from an OA animal is not an OA thread. Thread-level facts come
from the folder path, mussel-level facts from the key, joined on the tag. Never put the
animal's arm in the folder path; that creates a second copy of the key that can drift.

`phase` has three levels because the lab-reference animals are not pre-exposure baselines. A
binary before/after would pool them and contaminate every paired contrast.


## How to run

1. Open `thread-strength.Rproj` in RStudio. Not the repository-root `.Rproj`; `here::here()`
   must resolve to this folder.
2. `01_code/0_build_mussel_key.Rmd` — only needed when
   `../01_mussel-measurements/mussel-size-measurements.xlsx` changes.
3. `01_code/1_extract_tensometer_data.Rmd` — reads every trace, writes
   `03_analyses/01_extract-tensometer-data/thread-summary-raw-output.xlsx` and the QC plots.
4. `01_code/2_assemble_thread_summary.Rmd` — joins `02_data/pad_area_measurements.xlsx` and
   writes `03_analyses/02_assemble-thread-summary/thread-summary.xlsx` (sheet `data`). It
   overwrites that file on every run, so do not edit it by hand.
5. `01_code/3_analyze_thread_strength.Rmd` — adhesion (kPa) figures and the per-animal
   ANCOVA (day-3 level adjusted for the animal's baseline; each arm vs control).
6. `01_code/4_decompose_adhesion.Rmd` — the same ANCOVA on peak force, plaque area and
   extension separately. Read this alongside script 3; run it after script 3.

`09_gene-mechanics-correlation` reads this folder's outputs (the thread summary, script 4's
`mussel_response_classification.csv`, both `DATA_ancova_animals.csv` files and script 1's
extraction), so re-run its driver after any change here.

## The manual step

`pad_area` (plaque cross-sectional area, mm²) and `failure` are measured by hand from
microscope images of each plaque and recorded in `02_data/pad_area_measurements.xlsx`, one
row per trace. They cannot be derived from a force trace. Script 2 joins the measurements to
the traces; a trace with a blank `pad_area` is carried through and then excluded, with a
count, by scripts 3 and 4.

`adhesion_kpa = max_force / pad_area * 1000` is recomputed by script 3 rather than trusted
from a cached spreadsheet formula.
