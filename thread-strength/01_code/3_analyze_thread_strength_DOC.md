# `3_analyze_thread_strength.Rmd` documentation

Analyzes byssal-thread plaque adhesion (kPa) for *M. trossulus* before and after a 3-day
stress exposure, relative to a pre-exposure baseline and a day-3 control arm.

Run it from inside `thread-strength.Rproj`, after curating
`03_analyses/thread-summary.xlsx`.

---

## 1. Input

`03_analyses/thread-summary.xlsx`, sheet `data`: the hand-curated table. It is produced by
reviewing `03_analyses/assemble-thread-summary/thread-summary-candidate.xlsx` (script 2),
adding `pad_area` and `failure`, and dropping bad runs.

The previous version read `02_data/thread-summary-trossulus-clean.xlsx`, which was deleted
in commit `91dea50`. The script could not run at all.

### Schema normalization

The load chunk accepts the curated file under either the old or the new column names, so the
script keeps working across the rename:

| accepted | used internally |
|---|---|
| `mussel_ID` or `mussel` | `mussel` |
| `thread_num` or `thread` | `thread` |
| `thread_trt` or `treatment` | `thread_trt` |

`phase` is derived from `thread_trt` if the column is absent, since `thread_trt` determines
it completely. A missing required column is a hard stop that names what is missing and what
is present, rather than an obscure failure two chunks later.

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

The retired `group` column is not used. It was fully recoverable from `thread_trt`, and its
value `control` collided with both `mussel_trt == "control"` and
`thread_trt == "treatment_control"`.

---

## 3. `INCLUDE_LAB_REFERENCE`, the one configuration decision

At the top of the script, defaulting to `FALSE`.

The previous version built its baseline reference pool as `group == "control"`, described in
its own comment as "ALL day-0 baselines ... the SAME 54 mussels". That pool was actually
**43 pre-exposure animals plus 11 lab-reference animals** (T001–T011) that never entered the
experimental system.

| | pool | mean adhesion |
|---|---|---|
| pre-exposure baselines | 43 | 72.42 kPa |
| lab-reference animals | 11 | 79.40 kPa |
| pooled (old behaviour) | 54 | 73.84 kPa |

The lab animals sit 6.98 kPa above the pre-exposure animals. On their own that gap is not
significant (Welch p = 0.44), but pooling them both raises the reference mean and adds 11
observations to it, and both effects push the two-sample p-values down:

| arm | `welch_p`, pre only (default) | `welch_p`, pre + lab (old) |
|---|---|---|
| `treatment_control` | 0.993 | 0.856 |
| OW | 0.0203 | 0.0103 |
| OA | 0.0263 | 0.0150 |
| DO | 0.469 | 0.400 |

**So the default setting roughly doubles the OA and OW two-sample p-values relative to your
existing figures.** Set the toggle to `TRUE` for one run if you want the old numbers side by
side.

Scope of the toggle:

- **Affects:** the line+box left cluster, and `welch_p` / `mannwhit_p` in section 1b.
- **Does not affect:** the paired tests (lab animals have no day-3 pull, so they never pair),
  or the mixed model, which is always restricted to `phase %in% c("pre", "post")`.

---

## 4. The control arm is no longer hardcoded as unpairable

The previous version asserted, in a comment and in its code path:

> The control arm has no within-subject pairs (baseline cohort T001-T012 differs from the
> day-3 control cohort T126-T135), so it has no paired test.

**That is false.** It was an artifact of pooling the lab-reference animals into the baseline.
All ten day-3 control animals, T126–T135, have baseline traces on disk in the pre-exposure
folder. They were never curated into `thread-summary.xlsx` because none of their 25 baseline
traces has a plaque measurement.

Pairing is now computed from the data, per arm, and reported by every code path that depends
on it. Current state:

| arm | day-3 mussels | paired to a baseline |
|---|---|---|
| `treatment_control` | 10 | **0** |
| OW | 22 | 12 |
| OA | 13 | 9 |
| DO | 11 | 9 |

Measuring those 25 traces takes the control arm to 10 pairs and gives the design a proper
within-subject control contrast. It also removes the rank deficiency described in §6.

A line+box panel is now built for the control arm as well, saved as
`LineBox_CONTROL_BaselineVsDay3.png`.

---

## 5. Figures

### Distribution panels `p1`–`p4`

x is now `thread_trt`, so each bar is one condition a thread was actually built in. The
previous version plotted a `treatment` column whose `control` bar **pooled** day-0/day-1
baselines with day-3 control-arm threads. `baseline` and `treatment_control` are now separate
bars, so these panels are not directly comparable to the old ones.

The palette gained `lab_control` (grey) and `treatment_control` (lighter blue); `baseline`
keeps the blue the old pooled `control` level used.

### Line + box panels

Unchanged in geometry. The left cluster is now the baseline reference pool as defined by
`INCLUDE_LAB_REFERENCE`, the x label reads `baseline` rather than `control`, and a panel is
built for every day-3 arm present in the data rather than a hardcoded three.

`make_line_box()` returns `NULL` and warns if an arm has no day-3 threads, instead of
producing an empty plot.

---

## 6. Statistics

**1a, between-arm day-3 comparison** on per-mussel means: one-way ANOVA, Tukey HSD, and a
Kruskal-Wallis backup. Skipped with a message if fewer than two arms have data.

**1b, within-arm before/after** on per-mussel means. Two questions per arm: a paired test on
animals measured at both timepoints, and a two-sample test of the arm's day-3 distribution
against the baseline pool. Run for **every** day-3 arm including control. Paired columns are
`NA` when there are fewer than two pairs, and the arms in that state are named in a message
pointing at the likely cause.

**2, mixed model** `lmer(adhesion_kpa ~ mussel_trt * timepoint + (1 | mussel))` on all
threads with `phase` in `pre`/`post`.

### The rank deficiency you will currently see

With the control arm at zero baseline threads, its interaction cell is empty. `lme4` drops a
coefficient, `lmerTest` warns about missing cells, and `emmeans` returns `nonEst` for the
control arm's within-arm change. **This is correct behaviour, not a model failure.** The
script now runs an explicit design-balance check before fitting and names the offending arm,
so the cause is stated rather than left to be inferred from an `lme4` message.

It resolves itself when the T126–T135 baselines are curated.

### Current results

Type III ANOVA, adhesion: `timepoint` F = 6.11, p = 0.014; `mussel_trt:timepoint`
F = 3.13, p = 0.045. On log adhesion: `timepoint` p = 0.0018, interaction p = 0.067.

Paired tests, unchanged by anything in this revision: OA p = 0.0030, OW p = 0.2314,
DO p = 0.9857.

The paired-versus-mixed disagreement documented in the open-items register still stands, and
the recommendation is unchanged: report the overall timepoint decline, which both models
agree on, and the paired within-subject tests as the primary stressor-specific evidence. Do
not cite mixed-model per-arm contrasts as stressor-specific claims while baseline pairing is
this unbalanced.

---

## 7. Outputs

Written to `03_analyses/analyze-thread-strength/`.

| file | contents |
|---|---|
| `BP_*.png` | four distribution panels |
| `LineBox_{CONTROL,OA,OW,DO}_*.png` | before/after panels, one per day-3 arm; since 17 September 2026 with the day-0 lab-reference animals as their own grey cluster at the far left (`SHOW_LAB_REFERENCE_IN_PANELS`), never pooled into the baseline |
| `DIAG_lmer_residuals.png` | **new.** Q-Q and residuals-vs-fitted for the adhesion model |
| `STATS_between_arm_day3_{ANOVA,TukeyHSD,KruskalWallis}.csv` | Kruskal-Wallis is **new** to disk |
| `STATS_within_arm_beforeafter.csv` | carries a `baseline_pool` column recording the toggle; since 17 September 2026 also `n_ctrl_post`, `delta_vs_ctrl`, `vs_ctrl_welch_p`, `vs_ctrl_mannwhit_p`: each stressor arm's day-3 animals against the **day-3 treatment control** (the reference the mixed model and the DESeq2 contrasts use), beside the pooled-baseline columns `welch_p` / `mannwhit_p`, which answer a different question |
| `STATS_lab_reference_check.csv` | **17 September 2026.** Per metric (adhesion, force, area, extension), per-mussel means: lab (day 0) vs pre-exposure (day 1), day-3 control vs pre-exposure, lab vs day-3 control; Welch and Mann-Whitney |
| `STATS_baseline_balance.csv` | **17 September 2026.** Per-mussel baseline by FUTURE arm (n, mean, sd, median), one-way ANOVA and Kruskal-Wallis per metric: the randomisation check behind the paired design |
| `STATS_repeatability.csv` | **17 September 2026.** ICC from the arm x timepoint mixed model's variance components (adhesion raw and log, log force, log area, extension) and the per-arm Pearson / Spearman correlation of the paired animals' baseline and day-3 per-mussel means |
| `STATS_lmer_typeIII_ANOVA.csv` | **new.** Both the raw and the log model (register items 5 and 6) |
| `STATS_lmer_fixed_effects.csv`, `..._log.csv` | the log version is **new** |
| `STATS_lmer_marginal_means.csv` | **new.** Arm means at post and at baseline |
| `STATS_lmer_within_arm_change.csv`, `STATS_lmer_arm_contrasts.csv` | unchanged |
| `RUN_provenance.txt` | Timestamp, toggle states, n at each filtering step, package versions |

Everything that was print-only is now on disk. That closes register items §1.5 (log-adhesion
model) and §1.6 (Type III ANOVA and emmeans marginal means).

---

## 8. Other changes

- **Rows with no `pad_area` are dropped loudly**, with a per-arm count, instead of passing
  through as `NA` and being silently discarded by each model in turn. The n behind every
  contrast is now visible.
- **An unrecognised `mussel_trt` value is a hard stop.** Previously `factor(levels = ...)`
  silently converted `lab_control` to `NA`, and those rows then vanished from the `lmer`
  without a word.
- **Five declared packages were dropped**: `bestNormalize`, `agricolae`, `nlme`, `multcomp`
  and `rstatix` are never called anywhere in the script. A fresh machine no longer installs
  them.
- The fate check is now a `stop()` that names the offending mussels rather than a bare
  `stopifnot`.

## 9. Known cosmetic issue

`p2` uses the `aes(fill = x, fill = after_scale(colorspace::lighten(fill, .5)))` idiom to
lighten the violin fill. ggplot2 3.5 emits `Duplicated aesthetics after name standardisation:
fill` for this. The panel still renders. This is pre-existing and was not changed, because
rewriting it would alter the appearance of a figure you have already been working from.

---

## 9. 17 September 2026: the two controls as outputs

Three checks that justify the design were by-products of the figures or lived only in
this document; they are now tables (Statistics 3 in the script, files above):

- **Lab-reference check.** The day-0 lab-reference animals, the day-1 pre-exposure
  animals and the day-3 treatment-control animals agree on adhesion (77.0, 77.3, 79.5 kPa;
  Welch p = 0.97 / 0.76), force (p = 0.37 / 0.49) and area (p = 0.27 / 0.49). Extension is
  the exception: the lab animals' threads extended further than the pre-exposure (0.19 vs
  0.14 mm, p = 0.02) and the day-3 control threads (p = 0.03).
- **Baseline balance by future arm.** Adhesion (control 89, OA 78, OW 74, DO 71 kPa;
  ANOVA p = 0.19) and force (p = 0.67) are balanced. **Plaque area is not** (control 2.8,
  OA 3.0, OW 3.7, DO 3.7 mm2; ANOVA p = 0.006, Kruskal p = 0.006), nor is extension
  (p = 0.02): the animals that went on to the warming and hypoxia arms had larger baseline
  plaques. The paired change in area (control +10 %, OA +7 %, OW -18 %, DO -37 %) therefore
  starts from unequal baselines, and part of the within-arm change in area can be regression
  to the mean. The between-arm comparison at day 3 (`STATS_decomp_typeIII_ANOVA.csv`,
  arm x timepoint interaction) is the estimate to lead with for area; quote the paired %
  change beside its baseline.
- **Repeatability.** ICC of the per-mussel random intercept is low: 0.21 (adhesion),
  0.23 (log adhesion), 0.18 (log force), 0.38 (log area), 0.13 (extension); the per-arm
  baseline-to-day-3 correlations of the paired animals are near zero for adhesion and force
  (control 0.00, OA -0.02, OW 0.29, DO 0.41 for adhesion). Pairing removes little
  between-animal variance on this dataset; the paired tests gain their power from the
  per-mussel means, not from the pairing itself, and the mixed model's random intercept is
  doing correspondingly little work.

`INCLUDE_LAB_REFERENCE = TRUE` now warns loudly; the day-0 animals appear in every line+box
panel as their own labelled cluster and enter no pool and no model.
