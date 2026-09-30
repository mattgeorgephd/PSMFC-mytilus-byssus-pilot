# `4_decompose_adhesion.Rmd`

Fits the per-animal ANCOVA script 3 fits to adhesion (`01_code/_ancova.R`) to peak force,
plaque area and extension at break **separately**. Exists because adhesion is a ratio, and
modelling its parts separately shows which component drives an arm difference in it.

## Why it exists

Adhesion (kPa) = peak force / plaque area, so on the log scale
`log(adhesion) = log(force) − log(area)`. Two arms with the same difference in adhesion can
differ in how force and area moved, so each component gets its own model.
`FIG_ancova_vs_control.png` draws adhesion (from script 3's output) beside force, area and
extension. Because each response is adjusted for its own baseline, the force and area
estimates are not constrained to combine exactly into the adhesion estimate.

## Model

    log(day-3 level) ~ arm + log(baseline level)     force, area
    day-3 extension  ~ arm + baseline extension      extension (mm)

One row per animal with threads at both timepoints; each animal's value is the mean of its
log thread values (the log of the geometric mean) for force and area, the plain mean for
extension. Arms were assigned at random; the day-3 treatment control is the reference arm.
OA, OW and DO are each compared with control, with 95% CI and
p unadjusted and Dunnett-adjusted over the three comparisons (emmeans `adjust = "mvt"`,
seeded). No other test of the arms is run.

The response is logged before fitting, so `emmeans` is told the transformation explicitly
(`tran = "log"`); without it `type = "response"` would silently stay on the log scale.

## Input

`03_analyses/thread-summary.xlsx`, sheet `data`, with the same schema normalisation as
script 3. Rows without `pad_area` are excluded. The figure also reads
`03_analyses/03_analyze-thread-strength/STATS_ancova_vs_control.csv`, so run script 3 first;
script 4 stops if that file is missing.

## Outputs, `03_analyses/04_decompose-adhesion/`

| file | contents |
|---|---|
| `STATS_ancova_tests.csv` | F for arm and for baseline, per response |
| `STATS_ancova_vs_control.csv` | each arm vs control: ratio (force, area) or difference in mm (extension), 95% CI and p, unadjusted and Dunnett |
| `STATS_ancova_adjusted_means.csv` | arm means at the mean baseline |
| `STATS_ancova_coefficients.csv` | lm coefficients with 95% CI |
| `DATA_ancova_animals.csv` | per-animal model data |
| `DIAG_ancova_residuals.png` | Q-Q and residuals vs fitted, per response |
| `FIG_ancova_vs_control.png` | each arm vs control for adhesion, force, area and extension |
| `mussel_response_classification.csv` | per-animal description, see below |
| `FIG_animal_response_force_vs_area.png` | every paired animal: its own change in force against change in area |
| `STATS_failure_mode.csv` | failure-mode counts by thread treatment; chi-square (with the minimum expected count) and a `fisher_exact_day3` row (Monte Carlo, B = 10,000, fixed seed) for the day-3 arms |
| `FIG_failure_mode_by_arm.png` | failure-mode composition by thread treatment |
| `RUN_provenance.txt` | timestamp, model, n per response |

## Per-animal response classification (handoff to gene-mechanics)

`mussel_response_classification.csv`, one row per animal with day-3 threads, change
columns populated for those with baselines (`paired = TRUE`). Per metric (adhesion, force,
area, extension): `_pre`, `_post`, `_change` (log-ratio; plain difference for extension),
`_pct_change`, `_direction` (`decreased` / `increased`), `_vs_control` (the animal's %
change minus the control arm's mean % change) and `_rel_direction`
(`below_control` / `above_control`).

Two reference frames because they answer different questions. The raw direction says
whether this animal's threads got weaker. The control arm's own change is the trajectory
in the system without a stressor, so the control-referenced frame asks whether an animal
fell *more than the control trajectory*, which is the stress-specific question.

Composite: `response_class` is `weaker` if both force and adhesion decreased, `stronger` if
both increased, else `mixed`; `response_score` is the mean standardised log-ratio across
force, area and adhesion (higher = held up better). Extension is reported but excluded from
the composite because its direction is not a strength direction.

The script prints the response class and the number of animals whose force, area,
adhesion and extension decreased, by arm. `FIG_animal_response_force_vs_area.png` shows
every paired animal on the two axes.

Script 20 in `gene-mechanics-correlation/` joins the class and the score into its paired
manifest for inspection; they enter no model there.

## Failure mode

Among the day-3 threads, the script tests whether failure-mode composition differs by arm:
chi-square, and Fisher's exact (Monte Carlo, B = 10,000, fixed seed), which is the one to
quote when expected counts are below 5 (the script flags this). It also prints mean and
median adhesion by failure mode, for interpreting any difference in composition.

## Palette note

`arm_colors` reuses the project palette so these panels sit beside script 3's figures.
That palette does not pass a colour-vision-deficiency check (the OW orange is too light,
the DO purple reads as grey, OA/OW are borderline under protanopia). In these figures the
arm is always labelled on the axis, so colour is never the only encoding. Worth revisiting
before the manuscript figures are finalised.
