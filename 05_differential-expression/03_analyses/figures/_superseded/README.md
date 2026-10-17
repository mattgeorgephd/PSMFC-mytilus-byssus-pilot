# figures/_superseded

Retired figures, kept as records.

| File | Made by | Retired | Replaced by |
|---|---|---|---|
| `PCA_foot_day3.png`, `PCA_gill_day3.png`, `PCA_both_tissues_day3.png`, `PCA_all_samples.png` | `03_deseq_contrasts.Rmd` until 2026-10-17 | 2026-10-17 | Step 14 (`14_pca.Rmd`, `../../pca/`) draws every sample set of the original analysis, each tissue on its own, with 95% ellipses, without and with the lab control, and compares them (PERMANOVA, spread). The four figures here were the day-3 panels by treatment (no ellipses, no lab control), both tissues on one plot and every library by region |
| `MT_nuclear_copy_share.png` | `13_mitochondrial_expression.Rmd` until 2026-10-02 | 2026-10-02 | Step 13 no longer sums each mitochondrial protein's gene and nuclear-scaffold copies in the genome count matrix; it tests reads counted on the mitochondrial genome alone (`05_sequence-alignment` step 05), where the copies cannot take them. The figure showed the share of each protein's reads that the copies held in the previous (StringTie + prepDE) matrix, per library. `../MT_haplotypes.png` shows why the genome alignment does not serve |
