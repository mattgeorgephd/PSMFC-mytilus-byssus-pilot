# 02_data

DESeq2 inputs.

| Item | Description | Source |
|------|-------------|--------|
| `gene_count_matrix_clean.csv` | Gene-level StringTie count matrix | sequence-alignment (HISAT2 + StringTie) |
| `transcript_count_matrix.csv` | Transcript-level count matrix | sequence-alignment |
| `treatmentinfo_clean.csv` | Sample-to-treatment design table | curated |
| `psmfc_mussel_rna_summary.csv` | RNA sample summary | sequencing submission |
| `PSMFC-mytilus-byssus-pilot-RNA-tagseq_raw.csv` | Raw sample/RNA metadata | sequencing submission |
| `gene_count_matrix_clean` (no extension) | Space-delimited copy of the gene count matrix, read by `07_enrichment/01_code/09-GOenrichment_listcreation_REVIGO.Rmd` | derived |

The gene count matrix is generated upstream from the StringTie ctabs in
`05_sequence-alignment/03_analyses/hisat/` and placed here as the DE input. The raw StringTie
matrix that `01_code/01_5-gene_count_matrix.Rmd` cleans (`gene_count_matrix.csv`) is not in
the repository, so the cleaned matrix it produced is the committed input.
