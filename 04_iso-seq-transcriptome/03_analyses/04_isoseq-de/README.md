# 04_isoseq-de

Written by `../../01_code/04_isoseq_de_comparison.Rmd`: the six TC contrasts of
`06_differential-expression` fitted on the Iso-Seq gene counts (`../03_salmon/`) exactly as 06
fits them on the genome counts (same samples, `~ treatment`, fit on all features, then at least
10 counts in a third of the samples, apeglm, DEG at padj < 0.05), with the mitochondrial
features left out, and compared gene by gene with the genome results.

| File | Contents |
|---|---|
| `<code>_isoseq_apeglm.csv` | every tested feature: `feature`, `baseMean`, apeglm `log2FoldChange`, `lfcSE`, Wald `pvalue`, `padj` |
| `isoseq_DEG_counts.csv` | per contrast: features tested and DEGs (up, down) on the Iso-Seq reference, with the genome numbers beside them |
| `reference_agreement.csv` | per contrast: DEGs on each reference, Iso-Seq DEGs on annotated genes, novel loci and unassigned isoforms, genes tested on both, DEGs on both, of those in the same direction, DEGs on one reference only, genome DEGs not tested on Iso-Seq or with no isoform, and the Spearman correlation of fold changes over genes tested on both |
| `counts_per_gene_both_references.csv` | every nuclear gene counted on both references: total count over the day-3 libraries on each and their ratio |
| `isoseq_only_DEGs.csv` | DEGs on the Iso-Seq reference only: whether on an annotated gene, whether the genome branch tested it, mean count, fold change and padj on each reference |
| `FIG_isoseq_vs_genome_lfc.png` | fold change on each reference for genes tested on both, per contrast, marked by where each is a DEG |
| `FIG_mapping_rates.png` | per library, salmon's mapping rate to the isoforms against HISAT2's alignment rate to the genome |
| `RUN_provenance.txt` | settings and input MD5s |

The step was checked by feeding it the genome count matrix in place of the Iso-Seq counts: it
then reproduces 06 exactly (same genes tested, same DEGs in all six contrasts).

## Result (2026-10-02)

| contrast | features tested (Iso-Seq / genome) | DEGs (Iso-Seq / genome) | genome DEGs also DEG on Iso-Seq | same direction | fold-change Spearman |
|---|---|---|---|---|---|
| Foot OA | 13,057 / 7,638 | 194 / 75 | 58 (77%) | 58 | 0.84 |
| Foot OW | 13,604 / 7,949 | 262 / 165 | 89 (54%) | 89 | 0.82 |
| Foot DO | 12,938 / 7,614 | 867 / 363 | 291 (80%) | 291 | 0.87 |
| Gill OA | 18,421 / 10,772 | 884 / 423 | 304 (72%) | 301 | 0.81 |
| Gill OW | 18,199 / 10,762 | 366 / 180 | 123 (68%) | 123 | 0.81 |
| Gill DO | 18,283 / 10,727 | 600 / 310 | 227 (73%) | 227 | 0.82 |

- **Agreement.** Over genes tested on both references, fold changes correlate at 0.81 to 0.87,
  and every genome DEG that is also an Iso-Seq DEG has the same sign but three (all in Gill OA).
  1 to 22 genome DEGs per contrast could not be tested on the Iso-Seq reference.
- **More DEGs on the Iso-Seq reference, mostly from more reads per gene.** On genes counted by
  both, the Iso-Seq count is a median 1.59 times the genome count (interquartile 1.17 to
  2.81); the genome counts are prepDE estimates, about the aligned read bases divided by 75,
  while the reads average about 63 bases, so equal read numbers would give a ratio of about
  1.19, and the Iso-Seq reference places about 1.3 times the reads on a typical gene (2.4 times or
  more on a quarter of genes). Of the Iso-Seq-only DEGs, 69 to 81% are on annotated
  genes; most of those were tested on the genome too, with lower counts there (median mean
  count 23 to 44 against 37 to 73 on Iso-Seq) and a near-miss padj (median 0.15 to 0.23); the
  rest the genome filter had dropped. Novel loci contribute 26 to 128 DEGs per contrast and
  unassigned isoforms 4 to 39.
- **Where the extra reads come from.** In `../03_salmon/read_classes_by_library.csv`, 18% of
  assigned reads are on isoforms that run past their gene's annotated 3' end and 14.5% on novel
  loci, which suggests the RefSeq gene models often end before the 3' ends that Tag-seq reads.
