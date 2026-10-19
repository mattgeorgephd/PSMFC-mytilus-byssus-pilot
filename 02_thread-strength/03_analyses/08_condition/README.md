# 08_condition

Written by `../../01_code/08_condition.Rmd` from
`../../../01_mussel-measurements/mussel-condition.xlsx` (sheet `data`, raw workbook, read only)
and script 01's mussel key (`../01_build-mussel-key/mussel-treatment-key.csv`, a check that none
of these mussels is a thread or gene-expression mussel).

Sixty mussels sampled for condition only: the day-0 lab reference T074 to T088 (15) and, at day
3, the treatment control T150 to T164 (15), OA T165 to T174, OW T175 to T184 and DO T185 to
T194 (10 each). Each was measured (shell length, height and width, mm; volume) and dissected,
and its gonad and remaining tissue weighed (g; dry masses according to the manuscript's Methods,
which the workbook does not label). Two indices, recomputed from these columns:

- condition index CI = 1000 x tissue mass (g) / shell length (cm)^3, in mg cm^-3
- gonad index GI = 100 x gonad mass / tissue mass, in % of the tissue mass (tissue mass is the
  gonad and the rest together)

They agree with the workbook's `condition_index` (2 decimals) and `gonad_index` (whole percent).
Each stressor arm is compared with the day-3 control in a one-way linear model of the index
(the model of record), unadjusted and Dunnett-adjusted over the three comparisons, with three
sensitivity analyses: the log index (ratios), an allometric ANCOVA that does not assume the
index's scaling (log tissue mass ~ arm + log shell length; log gonad mass ~ arm + log somatic
mass; Packard and Boardman 1988, Physiol Zool 61:1-9), and a Kruskal-Wallis test of the four
day-3 arms. The lab reference is described and drawn, not tested.

| File | Contents |
|---|---|
| `DATA_condition.csv` | every mussel: ID, `group` (lab, control, OA, OW, DO), `mussel_trt`, `day`, shell length, height, width and volume, tissue, gonad and somatic mass, CI, GI (unrounded) |
| `DESC_condition.csv` | n, mean, SD, median, quartiles and range of shell length, tissue mass, gonad mass, CI and GI, by group; no tests |
| `STATS_condition_tests.csv` | per index and model (`primary` marks the model of record): F for arm (and for the covariate, with its slope and SE), residual df, n, Shapiro-Wilk p of the residuals; Kruskal-Wallis chi-square and p |
| `STATS_condition_vs_control.csv` | per index and model, each stressor arm vs control: difference (index) or ratio (log models), % change for ratios, 95% CI, p unadjusted, and Dunnett CI and p |
| `DIAG_condition_residuals.png` | Q-Q and residuals vs fitted of the two models of record (CI, GI) |
| `FIG_condition_index.png` | CI by group: box plots (median, quartiles, whiskers to 1.5 IQR) and every mussel; the lab reference left of the dashed line; the Dunnett p of each stressor arm vs the day-3 control |
| `FIG_gonad_index.png` | GI, drawn the same way |
| `RUN_provenance.txt` | settings, mussels per group, whether the models agree, `checks` (every mussel is the day-0 lab reference or in a day-3 arm; the recomputed indices are the workbook's; none of the mussels is in the thread and gene-expression key; the figures draw every mussel), commit, input MD5s, R and package versions |

## Results (run of 9c46f77, 2026-10-19)

No stressor arm differed from the day-3 control in either index, in any of the four analyses.

- **CI.** Control 7.93 mg cm^-3 (SD 1.56), OA 7.90, OW 7.47, DO 7.28; arm F(3, 41) = 0.40,
  p = 0.75. Against the control: OA -0.03 (95% CI -1.42 to 1.36; Dunnett p > 0.99), OW -0.46
  (-1.85 to 0.93; p = 0.86), DO -0.65 (-2.04 to 0.74; p = 0.69). The log index (DO -11%,
  Dunnett p = 0.56), the allometric ANCOVA (DO -10%, p = 0.62) and the Kruskal-Wallis test
  (p = 0.98) agree. Tissue mass scaled with shell length to the power 3.38 (SE 0.25), close to
  the cube the index assumes.
- **GI.** Control 29.2% (SD 8.4), OA 24.6, OW 30.1, DO 29.1; arm F(3, 41) = 1.00, p = 0.40.
  Against the control: OA -4.5 percentage points (-11.0 to 1.9; Dunnett p = 0.38), OW +0.9
  (-5.5 to 7.4; p = 0.98), DO -0.04 (-6.5 to 6.4; p > 0.99). The log index (OA -16%, Dunnett
  p = 0.37), the allometric ANCOVA (OA -21%, p = 0.36) and the Kruskal-Wallis test (p = 0.41)
  agree. Gonad mass scaled with somatic mass to the power 0.96 (SE 0.11), close to the
  proportionality the index assumes.
- **Lab reference (described).** CI 8.01 (SD 1.18) and GI 28.5% (SD 3.9), like the day-3
  controls.
- **Precision.** With 10 mussels per stressor arm and 15 controls, the 95% CIs of the
  differences span about 1.4 CI units (18% of the control mean) and 6.5 GI points either side
  of the estimate, so smaller effects are not excluded.
- **Residuals.** The index models' residuals are close to normal (Shapiro-Wilk p = 0.19 for CI,
  0.052 for GI), more so than the log models' (0.019, 0.024), the reason the index is the model
  of record; every model gives the same answer.
