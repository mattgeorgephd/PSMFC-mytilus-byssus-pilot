# `3_analyze_thread_strength.Rmd` documentation

Analyzes byssal-thread plaque adhesion (kPa) for *M. trossulus* before and after a 3-day
stress exposure, relative to a pre-exposure baseline and a day-3 control arm, with the
day-0 lab-reference animals as a descriptive comparison group.

Run it from inside `thread-strength.Rproj`, after curating `03_analyses/thread-summary.xlsx`.
Every table it prints is also written to `03_analyses/analyze-thread-strength/`.

---

## 1. Input

`03_analyses/thread-summary.xlsx`, sheet `data`: the hand-curated table (380 traces from 86
animals; every trace has a `pad_area` and a `failure` mode). It is produced by reviewing
`03_analyses/assemble-thread-summary/thread-summary-candidate.xlsx` (script 2), adding
`pad_area` and `failure`, and dropping bad runs.

### Schema normalization

The load chunk accepts the curated file under either column spelling:

| accepted | used internally |
|---|---|
| `mussel_ID` or `mussel` | `mussel` |
| `thread_num` or `thread` | `thread` |
| `thread_trt` or `treatment` | `thread_trt` |

`phase` is derived from `thread_trt` if the column is absent, since `thread_trt` determines
it completely. A missing required column is a hard stop that names what is missing and what
is present. Rows with no `pad_area` are dropped with a per-arm count (currently none).
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

`timepoint` is the two-level modelling axis (`baseline` / `post`) derived from `phase`. It is
deliberately `NA` for `phase == "lab"`: those animals were never in the experimental system,
so they have no before/after position.

Current composition: 38 lab-reference threads (12 animals, day 0); 161 pre-exposure
threads (61 animals, day 1); 181 day-3 threads (11 control, 13 OA, 22 OW, 14 DO animals).
47 animals have threads at both timepoints (control 11, OA 12, OW 12, DO 12).

### The two controls

- **Treatment control** (`thread_trt == "treatment_control"`, `mussel_trt == "control"`):
  day-3 animals held three days under ambient conditions in the same system as the stressor
  arms. The reference for a stressor effect: the mixed-model interaction, the `vs_ctrl_*`
  columns in 1b, and the DESeq2 TC contrasts all use it.
- **Lab control** (`phase == "lab"`, `mussel_trt == "lab_control"`): day-0 animals that never
  entered the experimental system. They measure time in the system plus handling and are
  confounded with both, so they are a descriptive comparison group only (Statistics 3a),
  never a baseline and never in a model.

---

## 3. Configuration

`INCLUDE_LAB_REFERENCE <- FALSE`. Whether the day-0 lab-reference animals count as part of
the baseline reference pool used by the line+box left cluster and by the two-sample tests
in 1b. FALSE is the analysis of record; TRUE pools 12 animals that never entered the system
into the reference and confounds "stressor" with "being in the system at all", so the
script warns loudly when it is set. The paired tests and the mixed model never use the
pool, whatever the setting.

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
`p2` / `p4` restrict to the 47 repeated-measures animals. Palette: `lab_control` grey,
`baseline` blue, `treatment_control` lighter blue, `OA` green, `OW` orange, `DO` purple.

### Line + box panels, one per day-3 arm

Left cluster = the baseline reference pool (the same 61 pre-exposure animals on every
panel); right cluster = the animals with day-3 threads in that arm; a connecting line only
for animals present on both sides; a box beside each cluster; the lab-reference cluster at
the far left. `make_line_box()` returns `NULL` and warns if an arm has no day-3 threads. The
pair count is printed for every panel.

---

## 5. Statistics

**1a, between-arm day-3 comparison** on per-mussel means: one-way ANOVA (p = 0.016), Tukey
HSD (OW vs control p = 0.024, OA vs control 0.11, DO vs control 0.85) and a Kruskal-Wallis
backup.

**1b, within-arm before/after** on per-mussel means, three questions per arm:

| arm | pairs | paired t | Wilcoxon | vs pooled baseline (Welch) | vs day-3 control (Welch) |
|---|---|---|---|---|---|
| control | 11 | 0.28 | 0.28 | 0.76 | |
| OA | 12 | 0.0031 | 0.0068 | 0.0009 | 0.0064 |
| OW | 12 | 0.022 | 0.052 | 0.0003 | 0.0028 |
| DO | 12 | 0.67 | 0.57 | 0.50 | 0.44 |

The paired columns use only the paired animals; `welch_p` / `mannwhit_p` compare the arm's
day-3 animals with the pooled pre-exposure baseline of all arms (the "before" of the
before/after picture); `vs_ctrl_*` compare them with the day-3 treatment control (the
reference the mixed model and the expression contrasts use). The two references answer
different questions and the file names which is which (`baseline_pool`,
`vs_ctrl_reference` columns).

**2, mixed model** `lmer(adhesion_kpa ~ mussel_trt * timepoint + (1 | mussel))` on the 342
pre/post threads from 74 animals, plus the same model on log adhesion. Type III ANOVA:
`timepoint` p < 0.0001, `mussel_trt:timepoint` p = 0.013 (raw) and 0.0057 (log).
Within-arm change from `emmeans` (baseline minus post, kPa): control 9.0 (p = 0.26), OA 27.5
(p = 0.0001), OW 26.0 (p = 0.0002), DO 0.1 (p = 0.99). Residual Q-Q and residuals-vs-fitted
are saved (`DIAG_lmer_residuals.png`); adhesion is right-skewed and the log model is the
robustness check.

**3, the two controls, baseline balance and repeatability** (per-mussel means; adhesion,
force, area, extension):

- *3a lab-reference check* (`STATS_lab_reference_check.csv`): lab (day 0), pre-exposure
  (day 1) and day-3 control animals agree on adhesion (77.0, 77.3, 79.5 kPa; Welch
  p = 0.97 lab vs pre, 0.76 control vs pre), force (p = 0.37 / 0.49) and area (0.27 /
  0.49). Extension differs: the lab animals' threads extended further than the
  pre-exposure (0.19 vs 0.14 mm, p = 0.02) and the day-3 control threads (p = 0.03).
- *3b baseline balance by future arm* (`STATS_baseline_balance.csv`): arms were assigned
  after the baseline pull, so the per-mussel baseline should not differ by future arm.
  Adhesion (control 89, OA 78, OW 74, DO 71 kPa; ANOVA p = 0.19) and force (p = 0.67) are
  balanced. **Plaque area is not** (control 2.8, OA 3.0, OW 3.7, DO 3.7 mm2; ANOVA
  p = 0.006, Kruskal 0.006), nor is extension (p = 0.02): the animals that went on to the
  warming and hypoxia arms had larger baseline plaques, so part of their within-arm fall in
  area can be regression to the mean. Lead with the arm x timepoint interaction
  (`4_decompose_adhesion.Rmd`) for area and quote the paired % change beside its baseline.
- *3c repeatability* (`STATS_repeatability.csv`): ICC of the per-mussel random intercept
  from the arm x timepoint model is 0.21 (adhesion), 0.23 (log adhesion), 0.18 (log force),
  0.38 (log area), 0.13 (extension); the per-arm baseline-to-day-3 correlations of the
  paired animals are near zero for adhesion and force (control 0.00, OA -0.02, OW 0.29,
  DO 0.41 for adhesion). Pairing removes little between-animal variance on this dataset;
  the paired tests get their power from the per-mussel means, and the design argument for
  pairing is balance, not repeatability.

Parametric and non-parametric versions are reported side by side rather than relying on an
automatic transform. Sample sizes are small; results are descriptive of this pilot.

---

## 6. Outputs

Written to `03_analyses/analyze-thread-strength/`.

| file | contents |
|---|---|
| `BP_*.png` | four distribution panels |
| `LineBox_{CONTROL,OA,OW,DO}_*.png` | before/after panels, one per day-3 arm, with the lab-reference cluster |
| `DIAG_lmer_residuals.png` | Q-Q and residuals-vs-fitted for the adhesion model |
| `STATS_between_arm_day3_{ANOVA,TukeyHSD,KruskalWallis}.csv` | 1a |
| `STATS_within_arm_beforeafter.csv` | 1b: paired, pooled-baseline and day-3-control comparisons per arm, with `baseline_pool` and `vs_ctrl_reference` columns |
| `STATS_lmer_typeIII_ANOVA.csv` | raw and log adhesion models |
| `STATS_lmer_fixed_effects.csv`, `..._log.csv` | fixed-effect estimates |
| `STATS_lmer_marginal_means.csv` | arm means at post and at baseline |
| `STATS_lmer_within_arm_change.csv`, `STATS_lmer_arm_contrasts.csv` | `emmeans` contrasts |
| `STATS_lab_reference_check.csv` | 3a |
| `STATS_baseline_balance.csv` | 3b |
| `STATS_repeatability.csv` | 3c |
| `RUN_provenance.txt` | timestamp, settings, n at each filtering step, package versions |

---

## 7. Known cosmetic issue

`p2` uses the `aes(fill = x, fill = after_scale(colorspace::lighten(fill, .5)))` idiom to
lighten the violin fill. ggplot2 3.5 emits `Duplicated aesthetics after name standardisation:
fill` for this. The panel still renders.
