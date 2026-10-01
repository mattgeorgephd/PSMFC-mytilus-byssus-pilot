# mitochondrial

The mitochondrial protein genes analysed on their own, written by
`../../01_code/13_mitochondrial_expression.Rmd`. The genome contrasts leave these loci out
(`../count_matrix/mitochondrial_loci.csv`). Each protein is counted once, as the sum of its
mitochondrial-genome gene and its copies on unplaced scaffolds, because reads of one
transcript are split between them in proportions that vary between libraries.

| File | Contents |
|---|---|
| `mt_protein_counts.csv` | per library, the summed count of each of the 12 proteins (ND1-6, ND4L, CYTB, COX1-3, ATP6) and of all 12 (`all_12`) |
| `mt_de_TC.csv` | per TC contrast and protein: `baseMean`; `log2FoldChange`, `lfcSE`, `ci_low`, `ci_high` (DESeq2 Wald estimate, unshrunk, with 95% interval); `lfc_apeglm`, `lfc_apeglm_sd` (apeglm shrunk estimate and posterior SD); `pvalue` (Wald); `padj_mt` (BH over the 12 proteins within the contrast; NA for `all_12`, which is read from its own p); `significant`; `passes_count_filter` (the genome filter, at least 10 counts in a third of the samples; every protein is reported regardless); `copy_share` (share of the protein's reads on the copies); `outlier_replaced` (DESeq2 replaced an outlier count and refit the protein: the Wald estimate and p come from the refit, the apeglm estimate from the original counts; Gill OA ND1 and ATP6, from library T025G) |
| `mt_share_by_sample.csv` | per library: tissue, treatment, day, region, mitochondrial protein reads as a percentage of the library's nuclear reads (`mt_reads_pct_of_nuclear`, read by `09_gene-mechanics-correlation` script 05) and the share held by the copies |

Fits use the contrast's samples, `~ treatment` with the day-3 control as reference, size
factors from the nuclear genes only, and dispersions shrunk toward the transcriptome-wide
trend. Figures: `../figures/MT_mitochondrial_expression.png` (manuscript) and
`../figures/MT_nuclear_copy_share.png`.
