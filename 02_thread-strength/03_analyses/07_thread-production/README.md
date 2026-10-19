# 07_thread-production

Written by `../../01_code/07_thread_production.Rmd` from
`../../../01_mussel-measurements/mussel-thread-production.xlsx` (sheet `data`, raw workbook, read
only) and script 01's mussel key (`../01_build-mussel-key/mussel-treatment-key.csv`, a check that
every mussel's arm, group and days agree).

Both counts are threads made in 3 days: `threads_baseline` in the 3 days on a mica plate in the
shared holding system before the exposure (the workbook's `thread_count_before`), and
`threads_exposure` in the 3 days of the exposure, in the arm (`thread_count_after`). The analysis
of record is the per-animal ANCOVA of the thread mechanics, run by the same code (`_ancova.R`,
`fit_ancova()`): log(threads during the exposure) ~ arm + log(baseline threads), each stressor arm
against the day-3 treatment control, unadjusted and Dunnett-adjusted over the three comparisons.
A negative binomial GLM of the day-3 count with the same terms (`MASS::glm.nb`, log link) is the
sensitivity analysis. Only the 66 exposed mussels have both counts; the 12 day-0 lab-reference
mussels and the 29 mussels assigned to an arm they never entered have a baseline count only and
are described, not tested.

| File | Contents |
|---|---|
| `DATA_thread_counts.csv` | every mussel of the workbook (107): ID, species, `group`, `mussel_trt` (the arm), `days_in_trt`, `set` (lab reference, day 0; assigned, not exposed; exposed, day 3), the two counts (blank where there is none), `in_ancova` |
| `DESC_thread_counts.csv` | n, mean, SD, median, quartiles, range and geometric mean of each count, by set and arm; no tests |
| `STATS_ancova_tests.csv` | F for arm and for the baseline count, each adjusted for the other (`drop1`), residual df, n |
| `STATS_ancova_coefficients.csv` | lm coefficients with 95% CI |
| `STATS_ancova_adjusted_means.csv` | geometric mean count per arm at the mean log baseline (`baseline_at`, 9.96 threads) |
| `STATS_ancova_vs_control.csv` | each stressor arm vs control: ratio of geometric means, % change, 95% CI, p unadjusted, and Dunnett CI and p |
| `DATA_ancova_animals.csv` | the 66 mussels the model was fitted to: `y` and `baseline` (log counts) and the two counts |
| `DIAG_ancova_residuals.png` | Q-Q and residuals vs fitted of the ANCOVA |
| `STATS_negbin_tests.csv` | likelihood-ratio tests of arm and of the baseline count in the negative binomial GLM (each a `glm.nb` fit with and without the term, with its own theta), theta and its SE, and the dispersion of the Poisson GLM with the same terms (Pearson chi-square / df) |
| `STATS_negbin_vs_control.csv` | each stressor arm vs control in the negative binomial GLM: ratio of means, % change, 95% CI, p unadjusted, and Dunnett CI and p |
| `FIG_thread_production_boxes.png` | threads made at baseline (blue) and during the exposure (arm colour), one panel per arm: box plots (median, quartiles, whiskers to 1.5 IQR) and every mussel (spread sideways only); each stressor panel gives the ANCOVA's % difference from the control and its Dunnett p |
| `FIG_thread_production_paired.png` | the same counts and labels, a line per mussel |
| `RUN_provenance.txt` | settings, mussels per group, whether the two models agree, `checks` (arm, group and days are script 01's; every exposed mussel has both counts and no other has a day-3 count; every mussel has a baseline count; `fit_ancova()` fits the direct lm of the log counts; the figures draw the model's mussels), commit, input MD5s, R and package versions |

## Results (run of 9c46f77, 2026-10-19)

- **Mussels.** Control 12, OA 13, OW 22, DO 19. Five DO mussels (T113 to T117) have counts but no
  tensometer pull, so the DO arm is larger here than in scripts 04 to 06 (14). Baseline counts
  were alike across arms (geometric means 8.6 to 11.7 threads, medians 9.5 to 12).
- **Hypoxia.** DO mussels made 45% fewer threads during the exposure than the controls, adjusted
  for their baseline count: ratio 0.55 (95% CI 0.35 to 0.84; p = 0.0071, Dunnett p = 0.019);
  adjusted geometric means 4.6 and 8.4 threads. The overall test of arm: F(3, 61) = 9.00,
  p = 5.1e-5.
- **OA and OW.** Not different from the control: OA +49% (ratio 1.49, 0.94 to 2.37; p = 0.090,
  Dunnett p = 0.20), OW +14% (1.14, 0.75 to 1.73; p = 0.52, Dunnett p = 0.84).
- **Baseline.** The baseline count predicts the count during the exposure only weakly
  (F(1, 61) = 3.05, p = 0.086; slope 0.24 on the log scale).
- **Sensitivity.** The negative binomial GLM agrees in direction and significance for every arm:
  DO 0.55 (0.36 to 0.82; p = 0.0034, Dunnett p = 0.0090), OA 1.41 (Dunnett p = 0.21), OW 1.13
  (Dunnett p = 0.83); arm LR chi-square 23.0 on 3 df, p = 4.0e-5; theta 6.0, and the counts vary
  3.1 times as much as a Poisson's.
- **Described only.** The controls made about as many threads during the exposure as before it
  (geometric means 8.1 and 8.6); DO mussels made 4.7 against 11.7. The day-0 lab reference made
  9.4 at baseline (geometric mean), and the mussels assigned to an arm they never entered 8.6 to
  10.4.

An earlier analysis (the manuscript draft's thread production figures) compared each arm's counts
before and during the exposure without the day-3 control, which cannot separate a stressor from
3 more days in the laboratory; the analysis of record compares each arm with the day-3 control,
adjusted for the baseline count (`AGENTS.md`, statistics of record).
