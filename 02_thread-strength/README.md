# 02_thread-strength

Plaque adhesion strength and thread mechanics from tensometer pull tests, before and after a
3-day stress exposure (control, OA, OW, DO), plus a lab-reference group.

Self-contained and reproducible from this repository apart from one manual input, the
plaque measurements, noted below. It reads one file from outside this folder
(`../01_mussel-measurements/mussel-size-measurements.xlsx`, script 01) and sources the shared
`../tools/pipeline_checks.R` (scripts 04 and 05, for `RUN_provenance.txt`).

## Layout

```
02_thread-strength/
├── thread-strength.Rproj               open this first; it anchors here::here()
├── 01_code/
│   ├── 00_run_thread_strength.Rmd      batch runner: scripts 01-05 in order
│   ├── 01_build_mussel_key.Rmd         mussel workbook -> mussel key
│   ├── 02_extract_tensometer_data.Rmd  raw traces     -> trace-level extraction
│   ├── 03_assemble_thread_summary.Rmd  extraction + plaque areas -> thread summary
│   ├── 04_analyze_thread_strength.Rmd  thread summary -> adhesion plots + ANCOVA
│   ├── 05_decompose_adhesion.Rmd       thread summary -> force and area ANCOVA
│   ├── _ancova.R                       the per-animal ANCOVA scripts 04 and 05 source
│   └── *_DOC.md                        companion documentation, one per script
├── 02_data/
│   ├── tensometer_output/<phase folders>/   raw force/displacement .txt traces
│   └── pad_area_measurements.xlsx           hand-measured plaque area + failure mode, per trace
└── 03_analyses/
    ├── 01_build-mussel-key/                 script 01: mussel-treatment-key.csv (tag -> arm, species, rna flag)
    ├── 02_extract-tensometer-data/          script 02: thread-summary-raw-output.xlsx
    │   └── QC_plots/<source_folder>/        per-trace loess QC jpgs
    ├── 03_assemble-thread-summary/          script 03: thread-summary.xlsx, input to 04, 05 and 09
    ├── 04_analyze-thread-strength/          script 04: figures + STATS_ancova_*.csv
    ├── 05_decompose-adhesion/               script 05: force and area ANCOVA
    └── knit_html/                           runner reports and logs (git-ignored)
```

Each script writes only to its own numbered folder in `03_analyses/`, so the number on a
folder tells you which script made it.

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

Script 02 maps each folder to `thread_trt`, `phase` and `day` through its `folder_labels`
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

Open `thread-strength.Rproj` in RStudio (not the repository-root `.Rproj`; `here::here()`
must resolve to this folder) and knit `01_code/00_run_thread_strength.Rmd`. It runs scripts
01 to 05 in order, each in a fresh R process, and leaves an HTML report and a log per script
plus `run_log.csv` in `03_analyses/knit_html/`. Its `steps` parameter runs a subset (for
example `[4, 5]`). The repository-level `00_run_pipeline.Rmd` calls it as the first stage.

The scripts, in run order:

1. `01_build_mussel_key.Rmd`: reads
   `../01_mussel-measurements/mussel-size-measurements.xlsx` and writes
   `03_analyses/01_build-mussel-key/mussel-treatment-key.csv`.
2. `02_extract_tensometer_data.Rmd`: reads every trace, writes
   `03_analyses/02_extract-tensometer-data/thread-summary-raw-output.xlsx` and the QC plots.
3. `03_assemble_thread_summary.Rmd`: joins `02_data/pad_area_measurements.xlsx` and
   writes `03_analyses/03_assemble-thread-summary/thread-summary.xlsx` (sheet `data`). It
   overwrites that file on every run, so do not edit it by hand.
4. `04_analyze_thread_strength.Rmd`: adhesion (kPa) figures and the per-animal
   ANCOVA (day-3 level adjusted for the animal's baseline; each arm vs control).
5. `05_decompose_adhesion.Rmd`: the same ANCOVA on mean peak force, maximum peak force and
   plaque area separately. Read this alongside script 04; it reads script 04's output.

`09_gene-mechanics-correlation` reads this folder's outputs (the thread summary, script 05's
`mussel_response_classification.csv`, both `DATA_ancova_animals.csv` files and script 02's
extraction), so re-run its runner after any change here.

## The manual step

`pad_area` (plaque cross-sectional area, mm²) and `failure` are measured by hand from
microscope images of each plaque and recorded in `02_data/pad_area_measurements.xlsx`, one
row per trace. They cannot be derived from a force trace. Script 03 joins the measurements to
the traces; a trace with a blank `pad_area` is carried through and then excluded, with a
count, by scripts 04 and 05.

`adhesion_kpa = max_force / pad_area * 1000` is recomputed per thread by script 04 rather
than trusted from a cached spreadsheet formula.

## Force metrics and what is not measured

Each trace gives one thread's peak force (`max_force`, N, the column name at thread level).
Per animal and timepoint the threads are summarised two ways:

| metric | per-animal value | in the ANCOVA |
|---|---|---|
| `mean_force` | the animal's typical thread: the mean of its threads' peak forces | mean of the log peak forces (the geometric mean), the scale the model uses |
| `max_force` | the animal's strongest thread: the largest peak force | log of the largest peak force |

Descriptive tables (`DESC_*.csv`, script 04) give the arithmetic mean for `mean_force`; the
model uses the geometric mean, as it has since the per-animal ANCOVA was introduced. With
three threads per animal at day 3 (one animal has four) the two means differ little, but
they are not the same number. `max_force` rests on one thread per animal, so it is noisier
and carries the largest-of-n bias of a maximum: an animal with more threads tends to have a
larger maximum. Baselines had one to three threads (16 of 48 animals had one or two), so
`max_force` is exploratory in `09_gene-mechanics-correlation`; `mean_force` is primary.

**Extension and the area under the curve are not analysed.** The pulls ran at a constant
rate, so displacement could be turned into extension, but each thread was cut near the
junction of the plaque and the distal region, so the length of distal thread under test
varied from pull to pull. Extension therefore cannot be compared between threads, nor can the
area under the force-time curve (`integral`, N·s), which depends on the same length. Scripts
02 to 05 extract and test neither (the raw traces still hold the displacement channel);
peak force is the only quantity taken from a trace.
