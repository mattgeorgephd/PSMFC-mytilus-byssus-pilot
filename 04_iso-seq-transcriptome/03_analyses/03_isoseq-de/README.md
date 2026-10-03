# 03_isoseq-de

Written by `../../01_code/03_isoseq_de_comparison.Rmd`: the six TC contrasts of
`06_differential-expression` fitted on the Iso-Seq gene counts (`../02_salmon/`) exactly as 06
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

Against the genome analysis of record (featureCounts on the Iso-Seq-extended RefSeq
annotation, since 2026-10-02):

| contrast | features tested (Iso-Seq / genome) | DEGs (Iso-Seq / genome) | genome DEGs also DEG on Iso-Seq | same direction | fold-change Spearman |
|---|---|---|---|---|---|
| Foot OA | 13,057 / 7,908 | 194 / 161 | 93 (58%) | 93 | 0.87 |
| Foot OW | 13,604 / 8,138 | 262 / 211 | 124 (59%) | 124 | 0.85 |
| Foot DO | 12,938 / 7,797 | 867 / 678 | 483 (71%) | 483 | 0.89 |
| Gill OA | 18,421 / 10,773 | 884 / 584 | 381 (65%) | 380 | 0.84 |
| Gill OW | 18,199 / 10,678 | 366 / 227 | 157 (69%) | 157 | 0.84 |
| Gill DO | 18,283 / 10,673 | 600 / 404 | 284 (70%) | 284 | 0.85 |

- **Agreement.** Over genes tested on both references, fold changes correlate at 0.84 to 0.89,
  and every genome DEG that is also an Iso-Seq DEG has the same sign but one (Gill OA). 2 to 35
  genome DEGs per contrast could not be tested on the Iso-Seq reference. The genome DEGs that
  are not Iso-Seq DEGs are mostly near misses there (median Iso-Seq padj 0.10 to 0.13; 65 to
  77% below 0.2), and 658 of those 665 keep their direction.
- **More DEGs on the Iso-Seq reference, 1.2 to 1.6 times as many.** Both references now count
  reads, and on genes counted by both the Iso-Seq count is a median 1.25 times the genome count
  (interquartile 0.97 to 2.02). Of the Iso-Seq-only DEGs, 61 to 75% are on annotated genes;
  most of those (63 to 88%) were tested on the genome too, with lower counts there (median mean
  count 26 to 66 against 36 to 93 on Iso-Seq) and a near-miss padj (median 0.11 to 0.17); the
  rest the genome filter had dropped. Novel loci contribute 26 to 128 DEGs per contrast and
  unassigned isoforms 4 to 39.
- **Against the previous genome counts** (StringTie + prepDE, until 2026-10-02) the Iso-Seq
  reference had found about twice as many DEGs (194 to 884 against 75 to 423) and a median 1.59
  times the counts, of which a factor of about 1.19 was prepDE's estimate of read bases over 75
  for reads averaging 63 bases. Option B (`05_sequence-alignment` steps 05 and 06, and step 04 here) showed that the
  counter accounted for more of that gap than the 3' ends did (`../04_augmented-de/README.md`).
- **Where the extra reads come from.** In `../02_salmon/read_classes_by_library.csv`, 18% of
  assigned reads are on isoforms that run past their gene's annotated 3' end and 14.5% on novel
  loci. The genome counts of record already include the 3' extensions (`ext3`); the rest of the
  difference lies in the reads on novel loci and on the isoforms' other exons, which only the
  `full` annotation of `05_sequence-alignment` step 05 counts, and in reads with several placements, which salmon
  shares among isoforms and featureCounts leaves out.
