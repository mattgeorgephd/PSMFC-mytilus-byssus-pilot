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
