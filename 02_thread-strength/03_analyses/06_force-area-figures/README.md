# 06_force-area-figures

Written by `../../01_code/06_force_area_figures.Rmd`, from script 05's per-animal table
(`../05_decompose-adhesion/mussel_response_classification.csv`), the thread summary
(`../03_assemble-thread-summary/thread-summary.xlsx`, for the day-0 lab reference) and script
04's `../04_analyze-thread-strength/DESC_lab_reference.csv` (a check). Descriptive: the arm
comparisons of record are the baseline-adjusted ANCOVAs of scripts 04 and 05.

Script 05's own figure, `../05_decompose-adhesion/FIG_animal_response_force_vs_area.png`, is kept
as it was drawn (linear % axes, no ellipses). The figures here add a 95% ellipse around each group,
from a multivariate t distribution fitted to its points (`ggplot2::stat_ellipse(type = "t")`,
`MASS::cov.trob`; Venables and Ripley 2002, *Modern Applied Statistics with S*, 4th ed.), as in
`05_differential-expression` step 14's PCA. Both are drawn on the log scale, the scale of the
ANCOVAs, so that the ellipses are fitted there.

| File | Contents |
|---|---|
| `FIG_force_vs_area_change_ellipses.png` | every animal with baseline and day-3 threads (12 per arm): its change in mean peak force against its change in plaque area (the ratio day 3 / baseline on a log scale, labelled as % change), coloured by arm, shaped by script 05's response class, labelled with the mussel ID, with a 95% ellipse around each arm. The day-0 lab reference has no baseline, so it has no change to show |
| `FIG_force_vs_area_day3_levels.png` | one point per animal: the mean peak force and the mean plaque area of its day-3 threads (log axes), every animal with day-3 threads (control 12, OA 13, OW 22, DO 14), each arm with its ellipse; left without, right with the 12 day-0 lab-reference animals (their day-0 threads; light grey, dashed ellipse). Both panels share their axes |
| `RUN_provenance.txt` | settings, animals per group, `checks` (the lab reference's and the day-3 controls' per-animal means equal script 04's `DESC_lab_reference.csv`; the change on the axes is script 05's; the ellipses are those of the log ratios), commit, input MD5s, R and package versions |

Each animal is summarised as scripts 04 and 05 summarise animals: the arithmetic mean of its
threads' peak forces and the mean of their plaque areas, the per-animal values of script 05's
table (the ANCOVAs use geometric means).
