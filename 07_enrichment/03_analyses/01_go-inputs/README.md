# 01_go-inputs

Written by `../../01_code/01_go_inputs.Rmd`; read by every later step.

| File | Contents |
|---|---|
| `gene_annotation.tsv` | one row per gene in any contrast's apeglm table (12,905): `gene` (the DESeq2 ID), `LOC_ID`, best BLAST hit (`protein_name`), median reference-transcript `length`, `mt_encoded` (best hit is an mtDNA-encoded protein) and the GO IDs by ontology (`GO_BP`, `GO_MF`, `GO_CC`, ";"-separated, direct annotation only) |
| `gene_sets_summary.csv` | per contrast and direction: universe size, genes with BP annotation and with a length, DEGs, annotated DEGs, mitochondrially encoded protein DEGs |
