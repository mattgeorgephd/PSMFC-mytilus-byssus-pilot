# `3_analyze_thread_strength.Rmd` documentation

Analyzes byssal-thread plaque adhesion (kPa) for *M. trossulus* before and after a 3-day
stress exposure, relative to a pre-exposure baseline and a day-3 control arm, with the
day-0 lab-reference animals as a descriptive comparison group.

Run it from inside `thread-strength.Rproj`, after script 2.
Every table it prints is also written to `03_analyses/03_analyze-thread-strength/`.

---

## 1. Input

`03_analyses/02_assemble-thread-summary/thread-summary.xlsx`, sheet `data`: the thread
summary written by script 2, which already carries `pad_area` and `failure`.

### Schema normalization

The load chunk accepts the thread summary under either column spelling:

| accepted | used internally |
|---|---|
| `mussel_ID` or `mussel` | `mussel` |
| `thread_num` or `thread` | `thread` |
| `thread_trt` or `treatment` | `thread_trt` |

`phase` is derived from `thread_trt` if the column is absent, since `thread_trt` determines
it completely. A missing required column is a hard stop that names what is missing and what
is present. Rows with no `pad_area` are dropped, and the script prints the per-arm count of
dropped rows.
`adhesion_kpa = max_force / pad_area * 1000` is recomputed, never read from the sheet.

---

## 2. Data model

Three columns, three grains. Each is named for what it describes.

| column | grain | values |
|---|---|---|
| `thread_trt` | thread | `lab_control`, `baseline`, `treatment_control`, `OA`, `OW`, `DO` |
| `phase` | thread | `lab`, `pre`, `post` |
| `mussel_trt` | mussel | `lab_control`, `control`, `OA`, `OW`, `DO` |

`mussel_trt` is a **destiny** label. For a `pre` thread it describes the animal's future, not
its past: a baseline thread from an OA animal is not an OA thread.

`timepoint` is the two-level before/after axis (`baseline` / `post`) derived from `phase`. It is
deliberately `NA` for `phase == "lab"`: those animals were never in the experimental system,
so they have no before/after position.

### The two controls

- **Treatment control** (`thread_trt == "treatment_control"`, `mussel_trt == "control"`):
  day-3 animals held three days under ambient conditions in the same system as the stressor
  arms. The reference for a stressor effect: the reference arm of the ANCOVA and of the
  DESeq2 TC contrasts.
- **Lab control** (`phase == "lab"`, `mussel_trt == "lab_control"`): day-0 animals that never
  entered the experimental system. They measure time in the system plus handling and are
  confounded with both, so they are a descriptive comparison group only
  (`DESC_lab_reference.csv`), never a baseline and never in the model.

---

## 3. Configuration

`INCLUDE_LAB_REFERENCE <- FALSE`. Whether the day-0 lab-reference animals count as part of
the baseline reference pool drawn as the left cluster of the line+box panels. FALSE is the
setting of record; TRUE pools 12 animals that never entered the system into that cluster
and confounds "stressor" with "being in the system at all", so the script warns loudly when
it is set. The ANCOVA never uses the pool: lab animals have no day-3 pull.

`SHOW_LAB_REFERENCE_IN_PANELS <- TRUE`. Draws the lab-reference animals in every line+box
panel as their own grey cluster at the far left, labelled "lab (day 0)", never merged into
the baseline cluster. Ignored when `INCLUDE_LAB_REFERENCE` is TRUE.

`RUN_provenance.txt` records both settings, the n at each filtering step and the package
versions behind the files in the output folder.

---

## 4. Figures

### Distribution panels `p1`–`p4`

x is `thread_trt`, so each bar is one condition a thread was actually built in; `baseline`
and `treatment_control` are separate bars. `p1` / `p3` show all threads / per-mussel means;
`p2` / `p4` restrict to the repeated-measures animals (threads at both timepoints).
Palette: `lab_control` grey, `baseline` blue, `treatment_control` lighter blue, `OA` green,
`OW` orange, `DO` purple.

### Line + box panels, one per day-3 arm

Left cluster = the baseline reference pool (the same pre-exposure animals on every panel);
right cluster = the animals with day-3 threads in that arm; a connecting line only for
animals present on both sides; a box beside each cluster; the lab-reference cluster at
the far left. `make_line_box()` returns `NULL` and warns if an arm has no day-3 threads. The
pair count is printed for every panel.

---

## 5. Statistics

One analysis, the per-animal ANCOVA in `01_code/_ancova.R` (shared with script 4). Arms were
assigned at random, so the day-3 level is modelled with the animal's own baseline level as
a covariate:

    log(day-3 adhesion) ~ arm + log(baseline adhesion)        lm, one row per animal

Each animal's value at a timepoint is the mean of its log thread adhesions (the log of the
geometric mean). Only animals with threads at both timepoints enter; the day-3 treatment
control is the reference arm. Adhesion is positive and right-skewed, with spread that rises
with its level, so it is modelled on the log scale and arm effects read as ratios.

- `STATS_ancova_tests.csv`: F for arm (adjusted for baseline) and for baseline (adjusted for
  arm), from `drop1()`.
- `STATS_ancova_vs_control.csv`: OA, OW and DO each vs control, as a ratio of geometric
  means at a common baseline (`estimate`, `pct_change`), with the 95% CI and p unadjusted
  (`ci_low`, `ci_high`, `p_unadjusted`) and Dunnett-adjusted over the three comparisons
  (`ci_low_dunnett`, `ci_high_dunnett`, `p_dunnett`; emmeans `adjust = "mvt"`, the exact
  multivariate-t method, seeded so it reproduces).
- `STATS_ancova_adjusted_means.csv`: arm geometric means at the mean baseline
  (`baseline_at`).
- `STATS_ancova_coefficients.csv`: the lm coefficients, including the baseline slope.
- `DATA_ancova_animals.csv`: the per-animal rows the model was fitted to.
- `DIAG_ancova_residuals.png`: Q-Q and residuals vs fitted.

No other test is run. The baseline table by future arm (`DESC_baseline_by_arm.csv`) and the
lab-reference comparison (`DESC_lab_reference.csv`) are descriptive: n, mean, SD and median
of the per-mussel means for adhesion, force, area and extension. Under random assignment a
baseline difference between arms is chance by construction, so it is not tested; the
ANCOVA adjusts for it. Sample sizes are small; results are descriptive of this pilot.

---

## 6. Outputs

Written to `03_analyses/03_analyze-thread-strength/`.

| file | contents |
|---|---|
| `BP_*.png` | four distribution panels |
| `LineBox_{CONTROL,OA,OW,DO}_*.png` | before/after panels, one per day-3 arm, with the lab-reference cluster |
| `STATS_ancova_tests.csv` | arm and baseline F-tests |
| `STATS_ancova_vs_control.csv` | each stressor arm vs control, unadjusted and Dunnett |
| `STATS_ancova_adjusted_means.csv` | arm means at the mean baseline |
| `STATS_ancova_coefficients.csv` | lm coefficients with 95% CI |
| `DATA_ancova_animals.csv` | per-animal model data |
| `DIAG_ancova_residuals.png` | residual diagnostics |
| `DESC_baseline_by_arm.csv` | baseline by future arm, descriptive |
| `DESC_lab_reference.csv` | lab (day 0), pre-exposure (day 1) and control (day 3), descriptive |
| `RUN_provenance.txt` | timestamp, settings, n, package versions |

---

## 7. Known cosmetic issue

`p2` uses the `aes(fill = x, fill = after_scale(colorspace::lighten(fill, .5)))` idiom to
lighten the violin fill. ggplot2 3.5 emits `Duplicated aesthetics after name standardisation:
fill` for this. The panel still renders.
