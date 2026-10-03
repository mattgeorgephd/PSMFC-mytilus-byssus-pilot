# 01_go-inputs

Written by `../../01_code/01_go_inputs.Rmd`; `gene_annotation.tsv` is read by steps 02, 03 and 04.

| File | Contents |
|---|---|
| `gene_annotation.tsv` | one row per gene in any contrast's apeglm table (11,754): `gene` (the DESeq2 ID), `LOC_ID` (`gene_key()`, `tools/gene_ids.R`), best BLAST hit (`protein_name`, highest bitscore), median reference-transcript `length` (from the gene record where a reference gene has no transcript record; missing only for the StringTie-only genes), and the GO IDs of the best hit by ontology (`GO_BP`, `GO_MF`, `GO_CC`, ";"-separated, direct annotation only) |
| `gene_sets_summary.csv` | per contrast and direction: universe size, genes with BP annotation and with a length, DEGs, DEGs with BP annotation |
| `RUN_provenance.txt` | the GO release in `GO.db`, the thresholds, the code commit, R and package versions for every method in this folder, and an MD5 checksum of each input |

The mitochondrial loci are not in any universe (`05_differential-expression` leaves them out
of its fits); step 01 stops if one appears.
