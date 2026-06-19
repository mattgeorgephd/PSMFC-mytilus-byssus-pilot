# `3_analyze_thread_strength.Rmd` documentation

Analysis of byssal-thread plaque adhesion (kPa) in *Mytilus trossulus* before and after a 3-day stress exposure (ocean acidification OA, ocean warming OW, hypoxia DO), against a day-0 baseline and a day-3 control arm.

This file rewrites the prior `3_analyze_thread_strength.Rmd`: it corrects the comments throughout, rebuilds the per-arm line figures with flanking boxplots, and replaces the invalid statistics with correct models plus a new mixed-model section.

---

## 1. Input

- `02_data/thread-summary-trossulus-clean.xlsx`, sheet `data` (299 thread-level rows, 80 mussels, *M. trossulus* only).
- Output folder created at `03_analyses/analyze-thread-strength/` (idempotent).

Adhesion is recomputed in R as `max_force / pad_area * 1000`, so the script does not depend on the spreadsheet's cached `adhesion_kpa` formula column.

## 2. Data model (the key to the filtering)

Three columns are easy to confuse, and the prior code keyed everything off the wrong one:

| Column | Meaning |
|---|---|
| `group` | `control` = day-0 **baseline** (pre-stress); `treatment` = day-3 (**post-stress**). This is the true before/after axis. |
| `mussel_trt` | The arm each mussel was **assigned** to (its fate): control, OA, OW, DO. Each mussel maps to exactly one fate (verified in code with a `stopifnot`). |
| `treatment` | The condition under which a single **thread** was made. **Every baseline thread is labeled `control`** regardless of fate, and the day-3 control-arm threads are **also** labeled `control`. |

Consequence: `treatment == "control"` (180 threads, 64 mussels) is a mixture of (i) day-0 baselines from all four arms and (ii) day-3 control-arm threads (mussels T126-T135). All filtering in this script therefore uses `group` + `mussel_trt`, and a clean `timepoint` factor (`baseline`/`post`) is derived from `group`.

A second structural fact worth knowing: the **control arm has no within-subject pairs**. The baseline control cohort is T001-T012 and the day-3 control cohort is T126-T135, different animals. So there is no before/after pairing for the control condition itself.

## 3. What changed from the prior version

**Comments / labels.** Removed the stale "ChatGPT prompt" headers and the copy-paste artifacts from the oyster/ploidy study (`timepoint:trt:ploidy`, "individual oysters"). Fixed the wrong color comment (`#5DADE2` is blue, not gray). Corrected the YAML title (`2_` → `3_`).

**ggplot2 modernization.** `size` is deprecated for lines/borders since ggplot2 3.4.0; the theme and geom borders now use `linewidth`. Point `size` is unchanged (points still use `size`).

**Refactor.** The three near-identical line-graph chunks are collapsed into one function, `make_line_box()`, called once per arm. This removes the copy-paste and the duplicated-bug risk.

## 4. Figures

### 4a. Distribution panels (`p1`-`p4`)
Kept and cleaned. `p1` (all threads) and `p3` (per-mussel means) are the two manuscript distribution panels. `p2`/`p4` restrict to mussels measured at both timepoints. These use the `treatment` column on the x-axis, so their "control" bar pools day-0 baselines with day-3 control threads; they are descriptive only. `p2` now uses the treatment palette (lightened) rather than ggplot's default hues, for consistency.

### 4b. Line + box panels (`make_line_box()`, one per arm)
This is the figure you asked to rebuild. Per the chosen convention (interpretation **B**):

- **Left cluster ("control")** = **all** day-0 baselines (`group == "control"`). This is the same 54 mussels on every panel, a fixed population-baseline reference.
- **Right cluster (the arm)** = only mussels that produced threads in that arm at day 3 (13 OA, 22 OW, 11 DO).
- **Connecting line** = drawn only for mussels present on both sides (true repeated measures: 9 OA, 12 OW, 9 DO). Unpaired mussels appear as a single point with no line. These line counts were checked against the arm-fate-with-both-timepoints sets and match exactly.
- **Flanking boxplots** summarize the **per-mussel means** (matching the plotted points). The control box sits at x = 0.60 and the arm box at x = 2.40 (width 0.22), while the points and connecting lines are confined to a narrow band around x = 1 and x = 2 (± 0.08). The boxes are therefore always clear of the points and never cross the lines, as requested.

**Geometry note.** The x-axis is continuous (control = 1, the arm = 2) so the boxes can be offset. Each mussel gets one fixed horizontal offset (computed once over all mussels) so the control cluster looks identical across the three panels.

**Color.** Control points and box in control blue (`#5DADE2`); arm points and box in the arm color (OA green, OW orange, DO purple); connecting lines neutral grey so they read as links rather than adding a 4th color family. The prior per-mussel rainbow ramp is dropped; if you prefer it back, it is a one-line change (map `color` to `mussel` and restore the `colorRampPalette` call). Note that under interpretation B the rainbow would be ~54 near-identical shades on the control side, which is mostly noise, hence the switch.

**Error bars.** Per-mussel SE (`sd / sqrt(n)`). Singletons (1 thread) have SE = NA and draw no bar (`na.rm = TRUE`); there are 2 such mussels on the control side and a few per arm.

**Axis.** y from 0 to 210, breaks every 25 (matches the distribution panels). Verified no point or error-bar top is clipped (max plotted point 134.7, max mean+SE 151.9).

## 5. Statistics

The prior tests were invalid and were replaced. The "all treatments" model (`aov(response ~ trt)`) had no error term, so it treated every thread as independent (pseudoreplication) while its comment claimed otherwise. The per-arm models (`aov(response ~ trt + Error(ID/thread))`) used `thread` (values 1,2,3) as a nesting stratum, which it is not, fed in all arms' baselines as the "control" level, and used `aov(Error())`, which requires balanced complete cells that this partially paired design does not have.

Replacement is in two parts.

### 5a/5b. Classic tests on per-mussel means
Collapsing threads to one mean per mussel per timepoint makes the **mussel** the unit of replication and removes pseudoreplication without a mixed model.

- **Between-arm (5a):** one-way ANOVA + Tukey HSD on day-3 per-mussel means across control/OA/OW/DO (each mussel once, so groups are independent), with Kruskal-Wallis as a non-parametric backup.
- **Within-arm before/after (5b):** for each stress arm, a **paired** test on mussels measured at both timepoints (matches the connecting lines) and a **two-sample** test of the arm's day-3 means against the full baseline pool (matches the boxes). Parametric and non-parametric versions are reported side by side. The control arm has no paired test (no within-subject pairs).

### 5c. Mixed model (lmer), new section
`adhesion_kpa ~ mussel_trt * timepoint + (1 | mussel)` on the full thread-level data. The random intercept absorbs multiple threads per mussel×timepoint and the within-mussel correlation across timepoints; the **arm×timepoint interaction** is the scientific target (does the baseline→day-3 change differ by arm?). Reported: Type III ANOVA (Satterthwaite df via `lmerTest`), fixed-effect estimates, residual diagnostics, a log-adhesion robustness fit, and `emmeans` contrasts (within-arm change; arm contrasts at day 3; baseline arm means as a randomization check).

**Caveat (documented in the code):** because the control arm is not paired within-subject, its `timepoint` effect in this model is a between-cohort contrast, not a within-subject one. The clean within-subject evidence is the paired tests in 5b; the lmer adds the full-data, all-animals view.

### Key results (validated in Python; confirm against the R output)
- **Between-arm at day 3:** no difference (one-way ANOVA F = 0.92, p = 0.44). Between-individual variance is large.
- **Within-arm, paired:** **OA shows a significant adhesion decrease** (mean Δ = -19.4 kPa, paired t p = 0.003, Wilcoxon p = 0.008, n = 9). OW (Δ = -13.5, p = 0.23, n = 12) and DO (Δ = +0.2, p = 0.99, n = 9) are not significant.

The contrast between these two is the whole point of getting the model right: the naive between-arm test at day 3 misses the OA effect that the paired/within-subject analysis detects, because pairing removes the between-individual variance. The old broken `aov` would have reported the wrong story.

## 6. Outputs (written to `03_analyses/analyze-thread-strength/`)

Figures: `BP_AllMussels_AllThreads.png`, `BP_OnlyRepeatedMussels_AllThreads.png`, `BP_AllMussels_MeanThreads.png`, `BP_OnlyRepeatedMussels_MeanThreads.png`, `LineBox_OA_BaselineVsTreatment.png`, `LineBox_OW_BaselineVsTreatment.png`, `LineBox_DO_BaselineVsTreatment.png`.

Stats tables: `STATS_between_arm_day3_ANOVA.csv`, `STATS_between_arm_day3_TukeyHSD.csv`, `STATS_within_arm_beforeafter.csv`, `STATS_lmer_fixed_effects.csv`, `STATS_lmer_within_arm_change.csv`, `STATS_lmer_arm_contrasts.csv`.

## 7. Dependencies and how to run

Open in RStudio and knit, or run chunk by chunk. The package chunk auto-installs anything missing. New dependencies beyond the prior version: `emmeans`, `broom`, `broom.mixed`, and `colorspace` (the last was already used implicitly by `p2` but was never in the install list).

## 8. Open decisions / flags for review on your diff

1. **Control cluster = all baselines (interpretation B), per your choice.** This puts the same 54-mussel baseline on every panel, so each panel's control box is identical and most control points on the OA panel are not OA animals. If you later want the cleaner within-arm before/after (interpretation C: only that arm's baselines on the left), it is a one-line filter change in `make_line_box()` (`filter(timepoint == "baseline", mussel_trt == trt_code)`).
2. **Classic tests collapse to per-mussel means.** This is the correct fix for the pseudoreplication; it is a different modeling choice than the old thread-level `aov`. The thread-level structure is handled correctly (and using all data) by the lmer section.
3. **Two-sample test added alongside the paired test** in 5b, because your figure shows both the boxes (population-level) and the lines (within-subject). Drop it if you only want the within-subject test.
4. **Colors changed** to side-based (blue control / arm-color treatment, grey lines) from the per-mussel rainbow. Reversible (see section 4b).
5. **Output filenames changed** to `LineBox_*` to reflect the new content (they are no longer "OnlyRepeatedMussels").
6. **lmer reported on raw and log adhesion.** If the residual diagnostics show clear right-skew (likely), prefer the log fit's inference; the interaction term is reported for both so you can compare directly.
