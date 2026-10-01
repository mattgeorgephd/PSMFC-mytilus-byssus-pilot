# 03_analyses

| Subfolder | Produced by | Contents |
|-----------|-------------|----------|
| `hisat/` | `01_hisat_stringtie.Rmd` (HPC) | the reference StringTie tables (`t_data.ctab`, `e_data.ctab`, `i_data.ctab`, `e2t.ctab`, `i2t.ctab`) and MultiQC alignment-rate reports. The per-sample StringTie folders (GTF and ctabs) and BAMs are not committed (too large); `t_data.ctab` gives transcript lengths to `07_enrichment` and the transcript-to-gene map to `02_prepDE.Rmd` |
| `prepDE/` | `prepDE.py` (HPC) and `02_prepDE.Rmd` | `transcript_count_matrix.csv` and `gene_count_matrix.csv`, the input of `06_differential-expression` (README inside) |
| `fastqc/trimmed/`, `fastqc/untrimmed/` | FastQC | per-sample read quality reports (HTML + ZIP) |
| `_superseded/kallisto/` | `_superseded/07-kallisto.Rmd` | kallisto quant per sample; superseded by HISAT2 |
| `knit_html/` | the runner | HTML reports and logs (git-ignored) |
