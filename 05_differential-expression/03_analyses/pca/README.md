# pca

Written by `../../01_code/14_pca.Rmd`, from step 01's count matrix and sample table
(`../count_matrix/`, nuclear genes only, with each library's RNA isolation date from
`library_crosswalk.csv`), step 02's contrast samples (`../DEG_lists/contrast_samples.csv`), step
04's apeglm tables (`../DEG_lists/<Foot|Gill>/<T><X>_TC_apeglm.csv`, for the isolation-day check) and the
protein names of the annotation of record (`03_blast` step 05's `LOC_GO_list.txt`). Every sample
set the original analysis drew a PCA for, each with its own blind VST and the 500 genes with the
largest variance (as `DESeq2::plotPCA`), with and without the day-0 lab control (LC), which is
not a control anywhere (its feet were dissected differently). Ellipses are 95% for each group
(multivariate t, `ggplot2::stat_ellipse(type = "t")`); PERMANOVA (Anderson 2001) and the spread
test (PERMDISP, centroid form; Anderson 2006) compare the groups of each set. Axes are not
comparable between panels, except the two tissue panels of one row of `PCA_both_tissues.png`.

| File | Contents |
|---|---|
| `PCA_arms.png` | every day-3 arm (control, OA, OW, DO), foot and gill in separate panels; top row without, bottom row with the lab control |
| `PCA_pairwise_foot.png`, `PCA_pairwise_gill.png` | each stressor against the day-3 controls (the samples of the contrasts of record, step 02), without (top) and with (bottom) the lab control, and the lab control against the day-3 controls |
| `PCA_both_tissues.png` | one PCA of both tissues (tissue is PC1), each tissue drawn in its own panel zoomed to its range; top without, bottom with the lab control |
| `PCA_all_samples.png` | one PCA of every library, by sampled region (foot phenol gland to tip, rest of the foot at day 0, gill) |
| `PCA_comparison.png` | the share of the variation the groups explain (PERMANOVA R², with its permutation p) in every set, in the 500 most variable genes and in all genes, without and with the lab control |
| `PCA_comparison.csv` | one row per sample set: `set`, `tissue`, `comparison`, `lab_control`, `groups` (libraries per group), `libraries`, `ntop`, `pc1_pct` and `pc2_pct` (variance explained), `pc1_byssal_share` and `pc2_byssal_share` (share of the component's squared loadings on genes named as byssal proteins), `permanova_R2`, `permanova_F`, `permanova_p` (top genes), `permanova_R2_all_genes`, `permanova_p_all_genes`, `spread_F`, `spread_p` |
| `PCA_scores.csv` | PC1 and PC2 of every library in every set (`set`, `tissue`, `sample`, `mussel`, `group`, `isolation_date`), and of the PCA of every library (`set` "every library"); `09` step 07 reads the foot and gill "every arm" sets |
| `PCA_isolation_day.csv` | the RNA isolation day within arms, by tissue, in the 500 most variable genes and in all genes: `isolation_days` (the day or days of each arm), `arms_with_two_days`, `R2_day_within_arm` (the day's share of the variation after the arm), `p` (the days permuted among the libraries of each arm), `null_median`. In foot each arm was isolated on one day, so it cannot be tested there |
| `isolation_day_contrasts.csv` | the indirect check for foot: what the three stressor contrasts of each tissue share, since all three compare with the same controls (a difference between the controls' isolation day and the others' would be shared): `controls_isolated`, `controls_share_day_with`, `DEGs_OA` to `DEGs_DO`, `genes_tested_in_all_three`, `DE_in_all_three`, `same_direction`, and the correlations of the three contrasts' apeglm fold changes (`r_lfc_OA_OW`, `r_lfc_OA_DO`, `r_lfc_OW_DO`); gill, whose controls share their day with DO, is the comparison |
| `PCA_loadings.csv` | the 25 genes with the largest loadings on PC1 and PC2 of each set: `LOC_ID`, `loading`, `byssal_name` (named as a byssal protein), `protein_name` (best hit) |

## Results (run of 610f1dc, 2026-10-17; `isolation_day_contrasts.csv` from b275fc2)

PERMANOVA R² (share of the variation between the groups) and permutation p (999 permutations;
0.001 is the smallest possible), in the 500 most variable genes and in all genes:

| set | foot R² (p) | foot, all genes | gill R² (p) | gill, all genes |
|---|---|---|---|---|
| every day-3 arm | 0.077 (0.13) | 0.093 (0.001) | 0.102 (0.001) | 0.097 (0.001) |
| OA vs control | 0.046 (0.28) | 0.056 (0.004) | 0.094 (0.001) | 0.077 (0.001) |
| OW vs control | 0.048 (0.24) | 0.054 (0.017) | 0.051 (0.037) | 0.057 (0.001) |
| DO vs control | 0.055 (0.20) | 0.075 (0.001) | 0.088 (0.001) | 0.076 (0.001) |
| lab control vs control | 0.346 (0.001) | 0.212 (0.001) | 0.120 (0.001) | 0.166 (0.001) |
| every arm with the lab control | 0.271 (0.001) | 0.167 (0.001) | 0.137 (0.001) | 0.152 (0.001) |

- **Foot.** The arms do not separate in the 500 most variable genes, and only weakly in all
  genes (R² 0.05 to 0.09). Foot PC1 (21% of the variance) is the byssal axis: 56% of its squared
  loadings fall on genes named as byssal proteins (foot proteins 2, 10, 12 and 15, YGH-rich
  proteins, collagen-like proteins 5 to 7), and its scores follow the animals' mean mfp
  expression (Spearman 0.92, `09` step 07). The largest source of variation among day-3 feet is
  which animals express the byssal genes, not the treatment.
- **Gill.** The arms separate (R² 0.10, p 0.001), OA and DO from the controls more than OW
  (0.094 and 0.088 against 0.051). The spread differs between the arms (p 0.018), which also
  raises R². Gill PC1 (7%) is driven by five adjacent genes with no BLAST hit (LOC134690918 to
  LOC134690922).
- **The lab control** separates from the day-3 controls in both tissues, most in foot (R² 0.35,
  PC1 47%), as the different day-0 dissection predicts; with it, PC1 of every foot set is the
  day-0 against day-3 axis.
- **Both tissues together:** tissue is PC1 (79%; R² 0.80 for tissue and arm together).
- **RNA isolation day** (`PCA_isolation_day.csv`). In gill, where OA and OW each span two
  isolation days, the day within those arms explains 4.2% of the variation, exactly what
  shuffling the days gives (null median 4.2%; p 0.47 in the top 500 genes, 0.88 in all genes),
  and the controls and DO, isolated on the same day, separate as much as OA from the controls.
  In foot each arm was isolated on its own day (the controls on 2021-12-23, 17 to 19 days before the
  others), so no direct test is possible. The indirect check (`isolation_day_contrasts.csv`)
  finds no sign of a shared difference: a difference between the foot controls' isolation and
  the others' would show in all three foot contrasts, but only 7 genes are DE in all three (all
  in the same direction), against 21 in gill (20), and the foot contrasts' fold changes agree
  less than gill's (r 0.18, 0.50 and 0.23 for OA with OW, OA with DO and OW with DO, against
  0.29, 0.62 and 0.43). A modest effect of the isolation day on some genes cannot be excluded.
