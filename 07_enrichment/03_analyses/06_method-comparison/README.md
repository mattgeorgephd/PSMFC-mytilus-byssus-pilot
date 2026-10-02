# 06_method-comparison

How far topGO `weight01`, topGO `classic`, goseq and clusterProfiler agree on the 12 TC runs,
for each ontology; written by `../../01_code/06_method_comparison.Rmd`. Files end in the
ontology (`BP`, `MF`, `CC`).

| File | Contents |
|---|---|
| `method_counts_TC_<ont>.csv` | enriched terms per run and method |
| `method_agreement_TC_<ont>.csv` | per run and method pair: enriched in each, shared, Jaccard (when both enriched something), terms compared, Spearman correlation of p-values |
| `method_pair_summary_TC_<ont>.csv` | per method pair: runs where both, or only one, enriched terms; median Jaccard and Spearman |
| `consensus_terms_TC_<ont>.csv` | terms enriched by `weight01` and at least one method with FDR control (`classic`, goseq or clusterProfiler, BH < 0.05): the most robust results, and the `fdr_supported` flag of `09_gene-mechanics-correlation` script 05 |
| `method_comparison_TC_<ont>.png` | A, enriched terms per run, one small panel per method; B, the overlap of each method pair |

`classic` and clusterProfiler run the same test on the same propagated annotation (Spearman
1.00 in BP); goseq differs from them only by the length weighting (0.92 to 1.00).
