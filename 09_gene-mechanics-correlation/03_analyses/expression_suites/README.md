# expression_suites

Written by `../../01_code/06_expression_suites.Rmd`, once per tissue (`_F` foot, `_G` gill),
from script 01's paired animals, VST and candidate genes (`../gene_mechanics/`) and script 05's
gene sets (`../go_mechanics/`). The question: do combinations of genes ("suites") go with
strong or weak attachment, across arms (which includes the treatment effect) and within arms
(arm and the animal's own baseline held constant)?

| File | Contents |
|---|---|
| `suite_axes_<T>.csv` | one row per axis: `axis_id`, `type` (DEG program, GO term, co-expression module, expression PC), `label`, `n_genes`, `var_expl` (variance its score explains among its genes; for a PC, among all genes), `eta2_arm` (share of the score's variance between arms), `mean_control` to `mean_DO` (arm means of the score, in SD units) |
| `suite_axis_tests_<T>.csv` | one row per axis and metric: within arms (`slope`, `se`, `p_lm`, `partial_r` with `partial_r_lo` and `partial_r_hi`; ANCOVA `level_day3 ~ score + treatment + level_baseline`), across arms (`across_r`, `p_across`; the same without treatment), between arms (`between_r`, the correlation of the four arm means of the score and of the baseline-adjusted level, descriptive), `p_interaction` (slopes differ between arms), and BH q within metric over all axes (`q_lm`), within metric and axis type (`q_type`) and within tier (`q_family`) |
| `suite_axis_members_<T>.csv` | the genes of each axis with protein names: every member of a program or module with its weight in the score (and `kME`, its correlation with a module's score); the 50 genes with the largest loading for each expression component |
| `suite_scores_<T>.csv` | one row per paired animal: sample, mussel, treatment, the day-3 and baseline metrics, `<metric>_pct_control` (baseline-adjusted day-3 level as % of the control mean) for each primary metric, and every axis score (SD units) |
| `suite_prediction_<T>.csv` | one row per prediction scenario (primary metric x predictor set x across or within): animals, predictors (`candidates`, `top_variable`, `axes`, `programs` (the six DEG programs) and `axes_with_treatment`, the power check: the treatment's indicator columns added to all the axes, across arms only), out-of-sample `q2`, the permutation null's median and 95% range, `n_perm`, `p_perm`, and `treatment_only_q2` (the treatment alone as predictor, across arms) |
| `suite_prediction_null_<T>.csv` | every permutation's Q2, by scenario |
| `suites_between_within_<T>.png` | the six DEG programs against the primary metrics: arm means (+/- 1 SE) joined in arm order, and the slopes within each arm, with the within-arm r, p and q |
| `suites_axes_map_<T>.png` | every axis: across-arm against within-arm partial correlation, by axis type and number of genes; red ring for within-arm q < 0.1 |
| `suites_state_space_<T>.png` | the animals in the space of the DEG programs (principal components of their scores, with the programs as arrows) and of the first two expression components, coloured by baseline-adjusted mean peak force |
| `suites_prediction_<T>.png` | the prediction test against its permutations and the treatment alone, with the power check; the within-arm tests of each axis type against chance (QQ plots) |
| `RUN_provenance_<T>.txt` | settings, axis counts, the lowest within-arm q, the prediction Q2 range, the DEG programs' and the power check's Q2, `checks` (script 05's sets reproduce its partial correlations; the prediction test is the same serially as on the workers; and others), commit, input MD5s, R and package versions |

## Results (run of 73f8c4e, 2026-10-06)

- **Axes.** Foot: 108 (6 DEG programs, 44 GO terms, 48 co-expression modules of 30 to 425
  genes, 10 expression components), 10 to 10,017 genes each, 45 animals. Gill: 89 (6, 15, 58
  modules of 30 to 552 genes, 10), 10 to 13,017 genes, 46 animals.
- **Within arms, no axis tracks thread mechanics.** For mean peak force and plaque area the
  lowest q is 0.85 in foot and 0.89 in gill (0.62 and 0.83 over all four metrics); 2 of 216 and
  3 of 178 tests have p < 0.05, fewer than the 11 and 9 that chance alone would give, and the
  p-values of every axis type lie on or below the chance line. The DEG programs' within-arm
  partial r is -0.19 to 0.11 in foot and -0.22 to 0.22 in gill.
- **Across arms, the warming and hypoxia programs follow the treatment.** Their scores move with
  the arms' force: "Foot OW, up" has across-arm r -0.55 and between-arm r -0.85, "Foot OW,
  down" 0.47 and 0.79, and the gill OW programs -0.52 and 0.59 (between-arm -0.92 and 0.82),
  while their slopes within arms are flat (`suites_between_within_<T>.png`).
- **Prediction.** The six DEG programs predict the force of held-out animals across arms
  (foot Q2 0.18, permutation p 0.02; gill 0.12, p 0.05), less than the treatment alone (0.38
  in both), and not within arms (-0.01 and -0.04). The larger sets (262 or 324 candidate genes,
  the 2,000 most variable genes, all the axes) give Q2 from -1.00 to 0.03. One within-arm
  scenario beats its permutations: the gill candidates with plaque area (Q2 0.03, p 0.0099,
  the smallest p 100 permutations give), in line with script 01's HSP70-family gene (q 0.036);
  it is one of 32 scenarios (power checks aside), so it does not survive a correction over them.
- **Power check.** Given the treatment's columns together with all the axes, the elastic net
  predicts force with Q2 0.05 in foot and 0.07 in gill, against 0.38 for the treatment alone.
  With 45 or 46 animals it cannot find a signal of the treatment's size among about a hundred
  predictors, so the large sets' Q2 near or below 0 does not rule out modest multi-gene
  signals; the six programs, a handful of predictors, are the informative prediction test.

Programs were chosen because their genes differ between arms, so arms separate along them by
construction; the within-arm tests are the informative ones. With 10 to 12 animals per arm a
within-arm association needs a partial r of about 0.41 for p < 0.05 in one test (80% power) and
about 0.6 after correction for 100 tests. Exploratory.
