# 03_analyses

| Folder | Produced by | Contents |
|---|---|---|
| `goslims/` | `01_go_slims` | GO slim (BP) tables per TC contrast, `goslim_summary_TC.csv`, `goslim_TC_heatmap.png`, `goslim_provenance.txt` |
| `Top_gene_summaries/` | `02_uniprot_summaries`, `03_ortholog_lists` | NCBI gene summaries (`<code>_topgene_summs.csv`) and bivalve orthologs (`<code>_topgene_summs_ortho.csv`, `ortho_species.tab.gz`) of the top-50 DEGs, refreshed 2026-10-01, with `RUN_provenance_*.txt`. The top-50 lists themselves are written by `06_differential-expression` (`03_analyses/top_DEGs/Top_50_genes/`). README inside |
| `_superseded/` | the earlier GO slim script | `goslims_genome/`, superseded tables |
| `knit_html/` | the runner | HTML reports and logs (git-ignored) |
