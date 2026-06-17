# summary-plots

Cross-cutting presentation figures that combine results from several analysis folders:
standard metabolic rate (SMR), thread production, and physiological condition (CI and GI),
plotted by treatment (OA, OW, DO).

## Layout

```
summary-plots/
├── summary-plots.Rproj
├── 01_code/
│   └── plot_analysis_mussel.R
└── 03_analyses/
    └── plot-analysis-mussel/   generated plots (PNG) + manually assembled PPTX
```

There is no `02_data/`. Inputs come from other analysis folders through a `repo_root`
pointer, `normalizePath(file.path(here::here(), ".."))`:

| Read | From |
|------|------|
| `processed_summary.xlsx` (SMR) | `respirometry/03_analyses/metabolic-rate/output/` |
| `morphometrics.xlsx` (thread_count, condition) | `morphometrics/02_data/` |

Note: thread *production* (count) is stored in `morphometrics.xlsx`, not in
`thread-strength`. The `thread-strength` folder holds plaque *adhesion* strength, which
this script does not plot.

## How to run

Open `summary-plots.Rproj` so `here::here()` anchors here, then run the script. It writes
PNGs into `03_analyses/plot-analysis-mussel/` (the PPTX decks there are assembled by hand).

## Known issue: trailing feeding-rate / stats section

The final section of the script is leftover from the oyster heatwave project. It reads
`feeding_rate/FR.xlsx`, which is not in this repo, and uses oyster variables (`ploidy`
D/T, `species == "gallo"`). It produces no plots. It is kept intact per the no-delete
policy but will error if that section is run top to bottom. Either add `feeding_rate/`
data or remove the section.
