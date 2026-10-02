# mitogenome

Written by `../../01_code/05_mitogenome_counts.Rmd` (`online: true`): every library's trimmed
reads aligned to the mitochondrial genome alone (NC_007687.1; `../../02_data/`) with HISAT2
2.2.1, no spliced alignment, and counted per feature with featureCounts (Subread 2.1.1; sense
strand, the few reads on two adjacent features or with two placements shared fractionally).
`06_differential-expression` step 13 tests these counts; why they replace the genome
alignment's for the mitochondrial genes is in the step's text.

| File | Contents |
|---|---|
| `mitogenome_gene_counts.csv` | the 38 features (12 protein genes, 2 rRNAs, 24 tRNAs; rows in genome order, as `mitogenome_genes.saf`) x the 131 libraries (columns named by sample ID), HISAT2 default scoring: **of record** |
| `mitogenome_gene_counts_permissive.csv` | the same with a permissive minimum alignment score (`--score-min L,0,-0.6`, about six mismatches in a 63-base read instead of two): step 13's sensitivity setting |
| `mapping_summary.csv` | per library: reads, mean read length, reads aligned to the mitogenome with each setting, the default's share of the library, reads on the 12 protein genes and on the 2 rRNAs |
| `RUN_provenance.txt` | tool versions, settings, code commit and input MD5s |
| `libraries/<library>/` | per-library counts, featureCounts summaries and HISAT2 logs (git-ignored) |
