# 03_analyses

| Item | Produced by | Contents |
|------|-------------|----------|
| `gene_mechanics/` | scripts 20, 21 | paired manifest, handoffs, per-gene and module associations, permutation, detection-floor flags and leave-one-out influence tables, figures; `_F` and `_G` |
| `expr_tables/` | script 22 | RNA × thread manifest, top-25 up/down expression tables, sample metadata; `_F` and `_G` |
| `byssus_genes/` | script 23 | byssus/foot gene list with expression, category scores; `_F` and `_G` |
| `knit_html/` | driver 24 | rendered HTML of each script per tissue plus `run_log.csv`; git-ignored, regenerable |
| `foot_byss_gene_plot.pdf`, `gill_byss_gene_plot.pdf` | legacy script 11 | byssus-gene expression plots |

Everything here is regenerable by knitting `01_code/24-run_gene_mechanics_by_tissue.Rmd`. `vst_paired_<T>.csv` (7 to 11 MB each) is a handoff script 20 rewrites every run.
