# thread-strength

Plaque adhesion strength and thread mechanics from tensometer pull tests, before and after a
3-day stress exposure (control, OA, OW, DO), plus a lab-reference group.

Self-contained and reproducible from this repository apart from one manual step, noted below.

## Layout

```
thread-strength/
├── thread-strength.Rproj              open this first; it anchors here::here()
├── 01_code/
│   ├── 0_build_mussel_key.Rmd         morphometrics -> 02_data/mussel-treatment-key.csv
│   ├── 1_extract_tensometer_data.Rmd  raw traces    -> 03_analyses/extract-tensometer-data/
│   ├── 2_assemble_thread_summary.Rmd  raw output    -> curation-ready candidate
│   ├── 3_analyze_thread_strength.Rmd  curated table -> adhesion plots + stats
│   ├── 4_decompose_adhesion.Rmd       curated table -> force / area / extension models
│   └── *_DOC.md                       companion documentation, one per script
├── 02_data/
│   ├── tensometer_output/<phase folders>/   raw force/displacement .txt traces
│   ├── pictures/{control,treatment}/        microscope images, source of pad_area
│   └── mussel-treatment-key.csv             mussel tag -> arm, species, rna flag
└── 03_analyses/
    ├── thread-summary.xlsx                  curated table, input to script 3
    ├── extract-tensometer-data/             output of script 1: thread-summary-raw-output.xlsx
    │   └── QC_plots/<source_folder>/        per-trace loess QC jpgs
    ├── assemble-thread-summary/             output of script 2: candidate + pad-area worklist
    ├── analyze-thread-strength/             output of script 3: figures + STATS_*.csv
    └── decompose-adhesion/                  output of script 4: force / area / extension models
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
2. `01_code/0_build_mussel_key.Rmd` — only needed when the morphometrics workbook changes.
3. `01_code/1_extract_tensometer_data.Rmd` — reads every trace, writes
   `03_analyses/extract-tensometer-data/thread-summary-raw-output.xlsx` and the QC plots.
4. `01_code/2_assemble_thread_summary.Rmd` — writes the curation candidate.
5. Curate by hand: add `pad_area` and `failure`, drop bad runs against the QC plots, save as
   `03_analyses/thread-summary.xlsx` (sheet `data`).
6. `01_code/3_analyze_thread_strength.Rmd` — adhesion (kPa) figures and models.
7. `01_code/4_decompose_adhesion.Rmd` — the same models on peak force, plaque area and
   extension separately. Read this alongside script 3: on this dataset the stressor effect
   is in the components and cancels in the ratio.

## The manual step

`pad_area` (plaque cross-sectional area, mm²) and `failure_mode` are measured from the
microscope images in `02_data/pictures/`. They cannot be derived from a force trace, so step
5 above is genuinely manual. Script 2 carries forward every measurement already present in
`thread-summary.xlsx` and reports exactly which traces still need one, so the manual work is
only ever on the new traces.

`adhesion_kpa = max_force / pad_area * 1000` is recomputed by script 3 rather than trusted
from a cached spreadsheet formula.
