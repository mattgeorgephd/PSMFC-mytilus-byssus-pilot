# 03_analyses

Everything here is written by the scripts in `../01_code/` and rebuilt by
`00_run_differential_expression.Rmd`.

| Folder / file | Produced by | Contents |
|---|---|---|
| `count_matrix/` | `01_clean_count_matrix` | `gene_count_matrix_clean.csv`, `treatmentinfo_clean.csv` (with `region`), `library_crosswalk.csv`, `mitochondrial_loci.csv` (the rows left out of the genome fits) |
| `DEG_lists/contrasts.csv`, `contrast_samples.csv` | `02_define_contrasts` | the 7 contrasts (6 TC, `FG_TC`) and the samples in each |
| `DEG_lists/filter_summary.csv` | `03_deseq_contrasts` | samples, genes in the matrix and genes kept by the count filter, per contrast |
| `DEG_lists/Foot/`, `DEG_lists/Gill/`, `DEG_lists/Foot_vs_Gill/` | `04_shrinkage_filtration` | per contrast: `<code>_apeglm.csv` (every gene kept by the filter), `<code>_siggene.csv` (padj < 0.05), `<code>_filter_counts.csv`, `<code>_MA_plots.pdf` |
| `DEG_lists/DEG_counts.csv` | `04_shrinkage_filtration` | genes tested and DEGs (all, up, down) per contrast |
| `DEG_lists/sensitivity_fourlevel/` | `05_fourlevel_sensitivity` | one four-level model per tissue against the pairwise TC lists |
| `DEG_lists/GOterms_genome/` | `06_join_annotation` | TC DEG lists joined to the BLAST / UniProt / GO annotation (`_sigs_merged`, `_sigs_ID`, `_sigs_unID`); `clean_zenodo_files/` from `12_deg_list_cleanup` |
| `DEG_lists/DEG_join_summary.csv` | `06_join_annotation` | per TC contrast: DEGs, annotated and unannotated genes, mitochondrial DEGs (a check, 0). Report `n_DEG_genes`, not merged rows |
| `mitochondrial/` | `13_mitochondrial_expression` | the 12 mitochondrial proteins on their own: summed counts, DESeq2 per TC contrast, mitochondrial share per library (README inside) |
| `dds/` | `03_deseq_contrasts` | the fitted, filtered `DESeqDataSet` of each contrast (git-ignored) |
| `figures/` | `03`, `08`-`11`, `13` | PCA, DEG counts, volcano, Venn and mitochondrial figures |
| `top_DEGs/Top_50_genes/` | `07_top_degs` | per TC contrast, the 25 most up- and 25 most down-regulated annotated DEGs |
| `_superseded/` | the previous scripts | per-contrast count and sample tables; the retired LC contrasts (`LC_contrasts/`); kept as a record |
| `knit_html/` | the runner | HTML reports and logs (git-ignored) |

Read by: `07_enrichment` (contrasts and apeglm tables), `08_gene-annotation`
(`GOterms_genome/`, `top_DEGs/`, `mitochondrial_loci.csv`) and `09_gene-mechanics-correlation`
(count matrix, sample table, `mitochondrial_loci.csv`, TC DEG lists, `GOterms_genome/`,
`mitochondrial/mt_share_by_sample.csv`).
