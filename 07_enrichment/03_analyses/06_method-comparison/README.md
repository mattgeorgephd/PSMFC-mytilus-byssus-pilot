# 06_method-comparison

How far topGO `weight01`, topGO `classic`, goseq and clusterProfiler agree on the 12 TC runs
(biological process); written by `../../01_code/06_method_comparison.Rmd`.

| File | Contents |
|---|---|
| `method_counts_TC_BP.csv` | enriched terms per run and method |
| `method_agreement_TC_BP.csv` | per run and method pair: enriched in each, shared, Jaccard (when both enriched something), terms compared, Spearman correlation of p-values |
| `method_pair_summary_TC_BP.csv` | per method pair: runs where both, or only one, enriched terms; median Jaccard and Spearman |
| `consensus_terms_TC_BP.csv` | terms enriched by `weight01` and at least one other method: the most robust results |
| `method_comparison_TC_BP.png` | enriched terms per method and run, and the overlap of each pair |

`classic` and clusterProfiler run the same test on the same propagated annotation (Spearman
1.00); goseq differs from them only by the length weighting (0.94 to 0.996).
