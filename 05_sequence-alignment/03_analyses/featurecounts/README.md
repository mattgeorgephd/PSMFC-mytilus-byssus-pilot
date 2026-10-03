# featurecounts

The count matrix of record, written by `../../01_code/04_count_matrix_of_record.Rmd`.

| File | Contents |
|---|---|
| `gene_count_matrix.csv` | genes (47,806) x the 131 libraries (columns named by sample ID, `T001F`); counts of uniquely aligned reads on the sense strand, by featureCounts (Subread 2.1.1) on the RefSeq annotation with Iso-Seq-extended 3' ends (`04_iso-seq-transcriptome` steps 05 and 06). Rows are named as prepDE names them (`gene-LOC134721619|LOC134721619`; `tools/gene_ids.R`), so `gene_key()` joins them to the annotation. `06_differential-expression` step 01 reads it |
| `RUN_provenance.txt` | code commit and input MD5s (step 06's matrix and its provenance, `t_data.ctab`) |

The alignment (HISAT2 2.2.1, defaults and `--dta`, RefSeq splice sites given at alignment
time) and the counting are recorded in `04_iso-seq-transcriptome/03_analyses/06_genome-recount/`
(README and `RUN_provenance.txt`).
