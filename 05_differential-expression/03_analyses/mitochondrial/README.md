# mitochondrial

The mitochondrial protein genes analysed on their own, written by
`../../01_code/13_mitochondrial_expression.Rmd`. The genome contrasts leave the mitochondrial
loci out (`../count_matrix/mitochondrial_loci.csv`). Since 2026-10-02 each protein is counted on
the mitochondrial genome alone (`04_sequence-alignment/03_analyses/mitogenome/`), not as the sum
of its gene and its nuclear-scaffold copies in the genome count matrix; the reasons are in the
script's text.

| File | Contents |
|---|---|
| `mt_protein_counts.csv` | per library, the count of each of the 12 proteins (ND1-6, ND4L, CYTB, COX1-3, ATP6) on the mitogenome and of all 12 (`all_12`) |
| `mt_de_TC.csv` | per TC contrast and protein: `baseMean`; `log2FoldChange`, `lfcSE`, `ci_low`, `ci_high` (DESeq2 Wald estimate, unshrunk, with 95% interval); `lfc_apeglm`, `lfc_apeglm_sd` (apeglm shrunk estimate and posterior SD); `pvalue` (Wald); `padj_mt` (BH over the 12 proteins within the contrast; NA for `all_12`, which is read from its own p); `significant`; `passes_count_filter` (the genome filter, at least 10 counts in a third of the samples; every protein is reported regardless); `outlier_replaced` (DESeq2 replaced an outlier count and refit the protein: the Wald estimate and p come from the refit, the apeglm estimate from the original counts) |
| `mt_de_TC_sensitivity.csv` | the same test with the animal's haplotype group as a covariate (`~ haplotype + treatment`) and with the permissive alignment's counts: `setting`, `code`, `protein`, Wald `log2FoldChange`, `lfcSE`, `pvalue`, `padj_mt`, `significant` |
| `mt_haplotypes.csv` | per animal: libraries, the COX1 index (reads aligned uniquely to the mitogenome's COX1 in the genome alignment of record, over the count on the mitogenome alone; median over the animal's libraries), its spread across them, and the group (reference-like, group A, group B) |
| `mt_counts_vs_previous.csv` | per protein and haplotype group: libraries, and the median ratio of the count here to the previous count (the protein's mitogenome gene and copies summed in the StringTie + prepDE matrix) |
| `mt_share_by_sample.csv` | per library: tissue, treatment, day, region, haplotype group, and mitochondrial protein reads as a percentage of the library's nuclear reads (`mt_reads_pct_of_nuclear`, read by `09_gene-mechanics-correlation` script 05) |

Fits use the contrast's samples, `~ treatment` with the day-3 control as reference, size
factors from the nuclear genes of the count matrix of record only, and dispersions shrunk
toward the transcriptome-wide trend. Figures: `../figures/MT_mitochondrial_expression.png`
(manuscript) and `../figures/MT_haplotypes.png`.
