# figures

| File | Produced by | Shows |
|---|---|---|
| `PCA_foot_day3.png`, `PCA_gill_day3.png`, `PCA_both_tissues_day3.png` | `03_deseq_contrasts` | variance-stabilised counts of the day-3 samples by treatment |
| `PCA_all_samples.png` | `03_deseq_contrasts` | every library, coloured by sampled region (phenol gland to tip, rest of the foot, gill) |
| `TC_venn_gill_stressors.png`, `TC_venn_foot_stressors.png`, `TC_venn_foot_vs_gill_<X>.png` | `08_deg_venn` | overlap of the TC DEG lists |
| `TC_volcano_foot.png`, `TC_volcano_gill.png` | `09_volcano_plots` | the manuscript volcano panels (TC): the DEGs only (padj < 0.05), apeglm log2 fold change against adjusted p, one panel per stressor with its counts, axes shared and set from the data (`volcano_TC.png` adds the genes that are not DEGs) |
| `TC_DEG_numbers.png` | `10_number_degs` | TC DEGs by direction |
| `DEG_counts_all_contrasts.png` | `11_deg_figures_all_contrasts` | DEGs up and down for the 7 contrasts (TC and foot vs gill) |
| `volcano_TC.png`, `volcano_FG.png` | `11_deg_figures_all_contrasts` | apeglm log2 fold change against adjusted p for each family (red up, blue down) |
| `MT_mitochondrial_expression.png` | `13_mitochondrial_expression` | **manuscript figure for the mitochondrial genes**: A, each of the 12 proteins (counted on the mitochondrial genome alone) and their sum, each stressor against the day-3 control, Wald log2 fold change with 95% interval, filled where BH-adjusted p < 0.05; B, the mitochondrial share of each day-3 library |
| `MT_haplotypes.png` | `13_mitochondrial_expression` | the mitochondrial haplotype groups: per library, the share of COX1 reads that align uniquely to the mitogenome's COX1 in the genome alignment, against the count on the mitogenome alone (log scale) |

The LC figures (`volcano_LC.png`, `TC_vs_LC_overlap.png`) are retired to
`../_superseded/LC_contrasts/figures/`; `MT_nuclear_copy_share.png` to `_superseded/` here
(README inside).

Colours come from `tools/plot_style.R` in every figure: red up, blue down; control grey, OA
green, OW orange, DO purple.
