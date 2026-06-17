# thread-strength

Plaque adhesion strength and thread mechanics from tensometer pull tests, before
and after a 3-day stress exposure (control, OA, OW, DO).

This analysis is **self-contained and reproducible from this repository**: all
raw inputs are present in `02_data/`, with one manual step noted below.

## Layout

```
thread-strength/
├── thread-strength.Rproj      open this first; it sets the working directory
├── 01_code/
│   ├── 1_extract_tensometer_data.Rmd   raw .txt traces -> max force + AUC tables
│   ├── 2_analyze_thread_strength.Rmd   summary.xlsx     -> adhesion plots + stats
│   └── thread-code/                    (legacy note)
├── 02_data/
│   ├── tensometer_output/{control,treatment}/   raw force/displacement .txt files
│   ├── pictures/{control,treatment}/            microscope images for plaque area
│   ├── summary.xlsx                             curated input for script 2 (see note)
│   └── instrument-reference/                    force-gauge VI, manual, wiring notes
└── 03_analyses/
    ├── extract-tensometer-data/        outputs of script 1
    │   ├── QC_plots/{control,treatment}/   per-thread loess QC jpgs
    │   └── summarized_data/                max_force_*.xlsx, integral_*.xlsx
    └── analyze-thread-strength/        outputs of script 2 (boxplots, line graphs)
```

Output subfolders under `03_analyses/` are named for the script that generates
them, so a file's provenance is always obvious.

## How to run

1. Open `thread-strength.Rproj` in RStudio (this anchors `here::here()` to this
   folder).
2. Knit or run `01_code/1_extract_tensometer_data.Rmd`. It reads the raw traces
   from `02_data/tensometer_output/` and writes force/AUC tables and QC plots into
   `03_analyses/extract-tensometer-data/`.
3. Knit or run `01_code/2_analyze_thread_strength.Rmd`. It reads `02_data/summary.xlsx`
   and writes adhesion boxplots and line graphs into
   `03_analyses/analyze-thread-strength/`.

All paths use `here::here(...)`, so the scripts run regardless of where this
folder sits, as long as they are run inside `thread-strength.Rproj`.

## Reproducibility note: summary.xlsx

`summary.xlsx` is **assembled by hand**, not produced by script 1. It combines
the force/AUC values from script 1 with plaque-area measurements (`pad_area`,
derived from the `02_data/pictures/` microscopy) and treatment labels. The master
copy of these data lives in the project Google Sheet
(<https://docs.google.com/spreadsheets/d/1GxLnNJjjjZ8xhBzz8nD-eUdpOwg6UY7yicG7ER5YIOQ/edit>).
Because this step is manual, re-running script 1 alone does not regenerate
`summary.xlsx`; update it from the sheet if the underlying data change.
