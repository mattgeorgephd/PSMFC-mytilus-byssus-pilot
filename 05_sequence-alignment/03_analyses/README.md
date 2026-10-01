# 03_analyses

| Subfolder | Produced by | Contents |
|-----------|-------------|----------|
| `hisat/` | `01_hisat_stringtie.Rmd` (HPC) | the reference StringTie tables (`t_data.ctab`, `e_data.ctab`, `i_data.ctab`, `e2t.ctab`, `i2t.ctab`) and MultiQC alignment-rate reports. The per-sample StringTie folders (GTF and ctabs) and BAMs are not committed (too large); `t_data.ctab` gives transcript lengths to `07_enrichment` and the transcript-to-gene map to `02_prepDE.Rmd` |
| `prepDE/` | `prepDE.py` (HPC) and `02_prepDE.Rmd` | `transcript_count_matrix.csv` and `gene_count_matrix.csv`, the input of `06_differential-expression` (README inside) |
| `fastqc/trimmed/`, `fastqc/untrimmed/` | FastQC | per-sample read quality reports (HTML + ZIP) |
| `fastqc/multiqc_report_trimmed_merged_2022-08-09.html` | MultiQC 1.12 (FastQC module), August 2022; copied from gannet `panopea/PSMFC-mytilus-byssus-pilot/` | an earlier QC of the trimmed, lane-merged reads of the first 73 libraries (2.5 to 6.3 million reads, mean length 71 to 78 bp). The other 56 (T031 to T058, T110 to T131) were trimmed later (gannet copies dated 2022-11-30). The trimmed files were made upstream of this repository; who trimmed them and with which settings is not recorded (`byssus-exp-analysis/code/01-diff-exp-analysis.Rmd` downloads them ready-trimmed from a list, `download.txt`, that was not kept) |
| `_superseded/kallisto/` | `_superseded/07-kallisto.Rmd` | kallisto quant per sample; superseded by HISAT2 |
| `knit_html/` | the runner | HTML reports and logs (git-ignored) |
