# 04_clusterprofiler

clusterProfiler `enricher` (hypergeometric, BH), one run per contrast and direction with its
own universe, merged with `merge_result()` into a compareCluster view; written by
`../../01_code/04_clusterprofiler.Rmd`.

| File | Contents |
|---|---|
| `clusterprofiler_enriched.csv` | terms with BH-adjusted p < 0.05 in any run: count, gene ratio, background ratio, p, adjusted p, genes |
| `clusterprofiler_all_terms_TC_<BP,MF,CC>.csv` | every term with at least one DEG in the TC runs, one file per ontology; read by step 06 |
| `clusterprofiler_run_summary.csv` | per run: genes, DEGs, terms reported and enriched |
| `clusterprofiler_<TC,FG>_<BP,MF,CC>_dotplot.png` | `enrichplot::dotplot()` of the merged result: five terms per run with any enriched term |
