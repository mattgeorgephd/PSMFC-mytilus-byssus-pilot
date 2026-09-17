# 02_data

DESeq2 inputs.

| Item | Description | Source |
|------|-------------|--------|
| `gene_count_matrix_clean.csv` | Gene-level StringTie count matrix | sequence-alignment (HISAT2 + StringTie) |
| `transcript_count_matrix.csv` | Transcript-level count matrix | sequence-alignment |
| `treatmentinfo_clean.csv` | Sample-to-treatment design table | curated |
| `psmfc_mussel_rna_summary.csv` | RNA sample summary | sequencing submission |
| `PSMFC-mytilus-byssus-pilot-RNA-tagseq_raw.csv` | Raw sample/RNA metadata | sequencing submission |
| `gene_count_matrix_clean/` | Companion directory for the cleaned matrix | derived |
| `secretion_state.csv` | Per foot sample: byssal plaque-protein module score (mean log2 CPM of mfp-2, mfp-4, FP10, FP12, FP15, tyrosinase-like 1), `secretion_state` (`on`/`off`), the rule and threshold, per-gene log2 CPM. Keyed by `sample` and `mussel`. Read by `02_7` and by gene-mechanics script 20 (covariate only when `USE_SECRETION_STATE = TRUE`) | `01_code/01_7-secretion_state.Rmd` |

The gene count matrix is generated upstream from the StringTie ctabs in
`sequence-alignment/03_analyses/hisat/` and placed here as the DE input.
