# 03_analyses

One folder per script in `../01_code/`, each with its own README.

| Folder | Produced by | Contents |
|---|---|---|
| `01_go-inputs/` | `01_go_inputs` | `gene_annotation.tsv` (the annotation every method uses), `gene_sets_summary.csv` |
| `02_topgo/` | `02_topgo` | topGO `weight01` results of every run, the TC p-values of every tested term (BP, MF, CC), dot plots |
| `03_goseq/` | `03_goseq` | goseq results, the length-bias diagnostic, dot plots |
| `04_clusterprofiler/` | `04_clusterprofiler` | clusterProfiler results, compareCluster dot plots |
| `05_rrvgo/` | `05_rrvgo` | topGO terms (BP, MF, CC) grouped into clusters, parent-term plots |
| `06_method-comparison/` | `06_method_comparison` | counts, overlap and correlation between methods; consensus terms |
| `_superseded/` | the DAVID / REVIGO workflow | the submitted lists and what the web tools returned |
| `knit_html/` | the runner | HTML reports and logs (git-ignored) |
