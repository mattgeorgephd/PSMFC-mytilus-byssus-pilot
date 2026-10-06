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

Programs were chosen because their genes differ between arms, so arms separate along them by
construction; the within-arm tests are the informative ones. With 10 to 12 animals per arm a
within-arm association needs a partial r of about 0.41 for p < 0.05 in one test (80% power) and
about 0.6 after correction for 100 tests. Exploratory.
