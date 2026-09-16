# `4_decompose_adhesion.Rmd`

Models peak force, plaque area and extension at break **separately**, with the same
arm × timepoint mixed model script 3 fits to adhesion. Exists because adhesion is a ratio
and, on this dataset, the ratio hides the stressor effect.

## Why it exists

Adhesion (kPa) = peak force / plaque area. On the log scale that is exactly
`log(adhesion) = log(force) − log(area)`, so a fitted change in log-adhesion is the
difference of the fitted changes in log-force and log-area.

Script 3, on the full 375-thread dataset with the control arm paired, finds:

- `timepoint` main effect on adhesion: F = 18.7, p < 0.0001. Everything declines.
- `mussel_trt:timepoint` interaction on adhesion: F = 0.71, **p = 0.55**. No arm declines
  differently from any other, including control.

Read on its own that says "no stressor effect". Decomposed, the same data say something
else entirely:

| response (log scale) | arm × timepoint | control | OA | OW | DO |
|---|---|---|---|---|---|
| adhesion | p = 0.83 | −22% | −21% | −31% | −23% |
| **peak force** | **p = 0.0005** | −12% | −18% | **−44%** | **−51%** |
| **plaque area** | **p < 0.0001** | +12% | +7% | **−19%** | **−37%** |
| extension at break | p = 0.0026 | −0.034 mm | −0.018 | +0.017 | +0.001 |

Percentages are mixed-model post/baseline ratios within each arm. Bold entries are
p < 0.001.

Warming and hypoxia animals build **smaller plaques that hold with far less force**.
Control and acidification animals build **slightly larger plaques that hold with slightly
less force**. In every arm the two changes divide out to a similar ~20–30% fall in kPa,
which is why the ratio shows no interaction. The stressor signal is in the components.

This is consistent with the earlier top-line finding that warming was the most mechanically
damaging stressor and that stress acts on plaque size, and it now has the control arm to
anchor it.

## Input

`03_analyses/thread-summary.xlsx`, sheet `data`. Same schema normalisation as script 3;
same exclusion of the lab-reference animals from the models (no before/after position).
Rows without `pad_area` are excluded (currently none).

## Responses

| response | scale | why |
|---|---|---|
| `adhesion_kpa` | log | reference; reproduces script 3's log model |
| `max_force` | log | positive, right-skewed; changes read as ratios |
| `pad_area` | log | as above; and makes the decomposition additive |
| `max_displacement` | raw | extension at break, mm |

`emmeans` cannot see a `log()` wrapped inside a programmatically built formula, so the
transformation is declared explicitly (`tran = "log"`). Without that, `type = "response"`
silently returns the log scale. The first draft of this script had exactly that bug.

## Outputs, `03_analyses/decompose-adhesion/`

| file | contents |
|---|---|
| `STATS_decomp_typeIII_ANOVA.csv` | arm, timepoint, interaction; one block per response |
| `STATS_decomp_within_arm_change.csv` | post vs baseline per arm per response. For logged responses: `ratio`, its 95% CI, and `pct_change` |
| `STATS_decomp_paired_by_arm.csv` | per-animal paired change (mean of threads per timepoint), per arm per response, with 95% CI and paired t / Wilcoxon |
| `STATS_failure_mode.csv` | failure-mode counts by thread treatment; chi-square and Fisher for the day-3 arms |
| `FIG_paired_change_by_metric.png` | paired % change with 95% CI, one panel per response |
| `FIG_failure_mode_by_arm.png` | failure-mode composition by thread treatment |
| `RUN_provenance.txt` | timestamp and n at each step |

## Failure mode

Among day-3 threads, failure-mode composition differs by arm: chi-square X² = 16.6,
df = 9, p = 0.056; Fisher's exact (simulated) p = 0.044. Peeling threads hold at roughly
55 kPa against 80 kPa for cohesive failures, so a shift toward peeling is a shift toward
weaker attachment. DO has the most peeling (49%) and almost no tearing; the day-3 control
arm has the most cohesive failures (61%).

## Palette note

`arm_colors` reuses the project palette so these panels sit beside script 3's figures.
That palette does not pass a colour-vision-deficiency check (the OW orange is too light,
the DO purple reads as grey, OA/OW are borderline under protanopia). In these figures the
arm is always labelled on the axis, so colour is never the only encoding. Worth revisiting
before the manuscript figures are finalised.
