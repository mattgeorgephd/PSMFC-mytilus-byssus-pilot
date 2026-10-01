# 03_analyses

| Item | Produced by | Contents |
|------|-------------|----------|
| `gene_mechanics/` | scripts 01, 02 | paired manifest, handoffs, per-gene ANCOVA associations, module associations, detection-floor flags, animal reconciliation, leave-one-out influence and best-hits tables, figures; `_F` and `_G` |
| `expr_tables/` | script 03 | RNA × thread manifest, top-25 up/down expression tables, sample metadata; `_F` and `_G` |
| `byssus_genes/` | script 04 | byssus/foot gene list with expression, category scores; `_F` and `_G` |
| `knit_html/` | runner 00 | rendered HTML and a log of each script per tissue, plus `run_log.csv`; git-ignored, regenerable |
| `_superseded/foot_byss_gene_plot.pdf`, `gill_byss_gene_plot.pdf` | legacy script `01_code/_superseded/11-byssal_thread_by_sample.Rmd` | byssus-gene expression plots |

Everything here except `_superseded/` is regenerable by knitting
`01_code/00_run_gene_mechanics_by_tissue.Rmd`. `vst_paired_<T>.csv` (7 to 11 MB each) is a
handoff script 01 rewrites every run.
