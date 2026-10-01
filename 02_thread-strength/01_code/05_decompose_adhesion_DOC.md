# `05_decompose_adhesion.Rmd`

Fits the per-animal ANCOVA script 04 fits to adhesion (`01_code/_ancova.R`) to mean peak
force, maximum peak force and plaque area **separately**. Exists because adhesion is a ratio,
and modelling its parts separately shows which component drives an arm difference in it.

## Why it exists

Adhesion (kPa) = peak force / plaque area, so on the log scale
`log(adhesion) = log(force) − log(area)`. Two arms with the same difference in adhesion can
differ in how force and area moved, so each component gets its own model.
`FIG_ancova_vs_control.png` draws adhesion (from script 04's output) beside the two force
metrics and area. Because each response is adjusted for its own baseline, the force and area
estimates are not constrained to combine exactly into the adhesion estimate.

Extension at break is not analysed: each thread was cut near the junction of the plaque and
the distal region, so the length of distal thread under test differed between pulls and
extension cannot be measured (see the folder README).

## Model

    log(day-3 level) ~ arm + log(baseline level)     mean_force, max_force, pad_area

One row per animal with threads at both timepoints. The animal's value is the mean of its log
thread values (the log of the geometric mean) for `mean_force` and `pad_area`, and the log of
its largest thread peak force for `max_force`. `max_force` rests on one thread per animal
and timepoint, and a maximum grows with the number of threads (baselines had one to three),
so it is the noisier, secondary force metric. Arms were assigned at random; the day-3
treatment control is the reference arm.
OA, OW and DO are each compared with control, with 95% CI and
p unadjusted and Dunnett-adjusted over the three comparisons (emmeans `adjust = "mvt"`,
seeded). No other test of the arms is run.

The response is logged before fitting, so `emmeans` is told the transformation explicitly
(`tran = "log"`); without it `type = "response"` would silently stay on the log scale.

## Input

`03_analyses/03_assemble-thread-summary/thread-summary.xlsx`, sheet `data`, with the same schema normalisation as
script 04. Rows without `pad_area` are excluded. The figure also reads
`03_analyses/04_analyze-thread-strength/STATS_ancova_vs_control.csv`, so run script 04 first;
script 05 stops if that file is missing.

## Outputs, `03_analyses/05_decompose-adhesion/`

| file | contents |
|---|---|
| `STATS_ancova_tests.csv` | F for arm and for baseline, per response |
| `STATS_ancova_vs_control.csv` | each arm vs control: ratio, 95% CI and p, unadjusted and Dunnett, per response |
| `STATS_ancova_adjusted_means.csv` | arm means at the mean baseline |
| `STATS_ancova_coefficients.csv` | lm coefficients with 95% CI |
| `DATA_ancova_animals.csv` | per-animal model data |
| `DIAG_ancova_residuals.png` | Q-Q and residuals vs fitted, per response |
| `FIG_ancova_vs_control.png` | each arm vs control for adhesion, mean peak force, maximum peak force and plaque area |
| `mussel_response_classification.csv` | per-animal description, see below |
| `FIG_animal_response_force_vs_area.png` | every paired animal: its own change in force against change in area |
| `STATS_failure_mode.csv` | failure-mode counts by thread treatment; chi-square (with the minimum expected count) and a `fisher_exact_day3` row (Monte Carlo, B = 10,000, fixed seed) for the day-3 arms |
| `FIG_failure_mode_by_arm.png` | failure-mode composition by thread treatment |
| `RUN_provenance.txt` | timestamp, model, n per response |

## Per-animal response classification (handoff to gene-mechanics)

`mussel_response_classification.csv`, one row per animal with day-3 threads, change
columns populated for those with baselines (`paired = TRUE`). Per metric (`adhesion_kpa`,
`mean_force`, `max_force`, `pad_area`): `_pre`, `_post` (arithmetic mean of the threads; the
largest thread for `max_force`), `_change` (log-ratio), `_pct_change`, `_direction` (`decreased` / `increased`), `_vs_control` (the animal's %
change minus the control arm's mean % change) and `_rel_direction`
(`below_control` / `above_control`).

Two reference frames because they answer different questions. The raw direction says
whether this animal's threads got weaker. The control arm's own change is the trajectory
in the system without a stressor, so the control-referenced frame asks whether an animal
fell *more than the control trajectory*, which is the stress-specific question.

Composite: `response_class` is `weaker` if both mean peak force and adhesion decreased,
`stronger` if both increased, else `mixed`; `response_score` is the mean standardised
log-ratio across mean peak force, area and adhesion (higher = held up better). `max_force`
is reported but kept out of the composite, which already carries force once.

The script prints the response class and the number of animals whose force, area and
adhesion decreased, by arm. `FIG_animal_response_force_vs_area.png` shows
every paired animal on the two axes.

Script 01 in `09_gene-mechanics-correlation/` joins the class and the score into its paired
manifest for inspection; they enter no model there.

## Failure mode

Among the day-3 threads, the script tests whether failure-mode composition differs by arm:
chi-square, and Fisher's exact (Monte Carlo, B = 10,000, fixed seed), which is the one to
quote when expected counts are below 5 (the script flags this). It also prints mean and
median adhesion by failure mode, for interpreting any difference in composition.

## Palette note

`arm_colors` is `TREATMENT_COLORS` from `../../tools/plot_style.R`, the colours every figure in
the repository uses (control grey, OA green, OW orange, DO purple). They keep the hues of the
earlier palette at steps that pass a colour-vision-deficiency check (worst all-pairs deltaE 9.2
deutan, 16.7 normal vision); the earlier values failed it. The arm is also labelled on the
axis, so colour is never the only encoding.
