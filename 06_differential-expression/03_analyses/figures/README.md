# figures

| File | Produced by | Shows |
|---|---|---|
| `PCA_foot_day3.png`, `PCA_gill_day3.png`, `PCA_both_tissues_day3.png` | `03_deseq_contrasts` | variance-stabilised counts of the day-3 samples by treatment |
| `PCA_all_samples.png` | `03_deseq_contrasts` | every library, coloured by sampled region (phenol gland to tip, rest of the foot, gill) |
| `TC_venn_gill_stressors.png`, `TC_venn_foot_stressors.png`, `TC_venn_foot_vs_gill_<X>.png` | `08_deg_venn` | overlap of the TC DEG lists |
| `TC_volcano_foot.png`, `TC_volcano_gill.png` | `09_volcano_plots` | the manuscript volcano panels (TC) |
| `TC_DEG_numbers.png` | `10_number_degs` | TC DEGs by direction |
| `DEG_counts_all_contrasts.png` | `11_deg_figures_all_contrasts` | DEGs up and down for all 16 contrasts |
| `volcano_TC.png`, `volcano_LC.png`, `volcano_FG.png` | `11_deg_figures_all_contrasts` | apeglm log2 fold change against adjusted p for each family (red up, blue down) |
| `TC_vs_LC_overlap.png` | `11_deg_figures_all_contrasts` | DEGs found against the day-3 control, the day-0 control, or both |

Colours come from `tools/plot_style.R` in every figure: red up, blue down; control grey, OA
green, OW orange, DO purple.
