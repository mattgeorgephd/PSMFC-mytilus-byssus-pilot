# 07_augmented-de

Written by `../../01_code/07_augmented_de_comparison.Rmd` (option B, a parallel analysis): the
six treatment-control contrasts of `06` fitted on each of step 06's six count matrices (two
counters x three annotations) with `06`'s rules (same samples, `~ treatment`, fit on every
gene, at least 10 counts in a third of the samples, apeglm, DEG at padj < 0.05, the 310
mitochondrial loci left out), and compared with the record and with each counter's RefSeq
control. Given the record's own counts, the step reproduces `06` exactly (same genes tested,
same DEGs, fold changes equal to 1e-15; checked 2026-10-02).

| File | Contents |
|---|---|
| `<stringtie|featurecounts>_<refseq|ext3|full>_apeglm.csv.gz` | every tested gene of every contrast: `code`, `feature` (gene key, or `novel:...` in `full`), `baseMean`, apeglm `log2FoldChange`, `lfcSE`, Wald `pvalue`, `padj` |
| `deg_summary.csv` | per matrix and contrast: genes tested, DEGs (up, down), the record's DEGs, shared, shared in the same direction, record DEGs lost (and not tested here), new DEGs (and on novel loci), Spearman of fold changes over genes tested in both |
| `annotation_effect.csv` | per counter and contrast: `ext3` and `full` against that counter's RefSeq control (DEGs shared, gained, gained on novel loci, lost, fold-change Spearman) |
| `control_vs_record.csv` | the StringTie + prepDE RefSeq control against the record: per-library totals, mitochondrial share per library, day-3 gene totals |
| `mitochondrial_share.csv` | per matrix: median mitochondrial share of the libraries (%), the record's, and the Spearman correlation of the per-library shares with the record's |
| `counts_per_gene.csv` | per gene, the total count over the day-3 libraries in the record and on each matrix |
| `new_DEGs.csv` | DEGs of `ext3` and `full` (StringTie + prepDE) that are not DEGs in the record: whether on a novel locus, whether `ext3` extended the gene, whether the record tested it, both mean counts, fold changes and padj |
| `byssal_genes.csv` | the byssal structural genes (best BLAST hit named as in `09`) that are a DEG in any matrix or in the record, by contrast and source |
| `FIG_deg_counts.png` | DEGs per contrast: the record and the three annotations, per counter |
| `FIG_lfc_vs_record.png` | fold change on each annotation (StringTie + prepDE) against the record, genes tested in both |
| `RUN_provenance.txt` | settings, input MD5s |

## Result (2026-10-02)

DEGs per contrast (apeglm padj < 0.05):

| contrast | record (06) | StringTie + prepDE: RefSeq control / 3' extension / full | featureCounts: RefSeq / 3' extension / full |
|---|---|---|---|
| Foot OA | 75 | 83 / 93 / 118 | 163 / 163 / 170 |
| Foot OW | 165 | 180 / 202 / 208 | 208 / 211 / 208 |
| Foot DO | 363 | 377 / 416 / 489 | 636 / 676 / 705 |
| Gill OA | 423 | 410 / 431 / 528 | 588 / 583 / 596 |
| Gill OW | 180 | 175 / 198 / 228 | 227 / 227 / 222 |
| Gill DO | 310 | 299 / 306 / 340 | 398 / 404 / 441 |

- **The realignment alone (control against the record).** Per-library totals are 0.992 of
  the record's (0.988 to 0.996), the mitochondrial share per library ranks identically
  (Spearman 1.0), and day-3 gene totals correlate at 0.997. The control keeps 91 to 99% of the
  record's DEGs per contrast (fold-change Spearman 0.993 to 0.996) and adds 6 to 25. About 5%
  of the DEG calls therefore move with the alignment alone; differences smaller than that say
  nothing about the annotation.
- **3' extension (StringTie + prepDE, the record's counter).** Against its control: 12 to 49
  DEGs gained and 2 to 16 lost per contrast, fold-change Spearman 0.98 to 0.99. Against the
  record: 91 to 97% of the record's DEGs kept, 21 to 67 new, of which 15 to 56 on extended
  genes; the new ones were near misses in the record (median record padj 0.06 to 0.09; 76 to
  95% below 0.2) with similar counts there (median baseMean ratio 1.0 to 1.2). Gene keys, names
  and GO annotation are unchanged.
- **Full models (StringTie + prepDE).** Against its control: 50 to 178 gained, 15 to 62 lost
  (Spearman 0.94 to 0.95). Against the record: 68 to 84% of the record's DEGs kept, 57 to 183
  new, 7 to 37 of them on novel loci, which have no names or GO terms.
- **The counter matters more than the annotation.** On the same alignment, featureCounts
  finds 16 to 96% more DEGs than StringTie + prepDE on the RefSeq annotation (4 to 75% more on
  `ext3`; -3 to 44% on `full`). Its counts are about 1.16
  times higher on a typical gene (prepDE estimates read bases over 75, while the reads average
  about 63 bases) and its gene-wise dispersions about 9% lower (Foot DO: median 0.065 against
  0.072), so it has more power; it also counts on the sense strand only and leaves out reads
  aligned to several places. Within featureCounts, the annotation changes little: `ext3` gains
  9 to 51 and loses 9 to 37 DEGs, `full` gains 31 to 119 and loses 24 to 68.
- **Mitochondrial reads.** With StringTie + prepDE the mitochondrial share of each library
  matches the record (Spearman 1.0 on RefSeq, 0.999 on `ext3`, 0.992 on `full`). featureCounts
  keeps only about a third of the mitochondrial reads (median share 0.8% against 2.3%;
  Spearman 0.86 on every annotation), most likely because reads of the mitochondrial genes
  also align to their nuclear copies and featureCounts leaves out reads aligned to several
  places; as configured it would not serve the mitochondrial analysis (`06` step 13).
- **Byssal genes** (`09`'s definition): the three ACDC genes are DEGs in Gill OA and Gill DO
  in the record and in every matrix. featureCounts adds foot proteins 11 and 13 in Foot DO
  with fold changes near 0 and small padj, the pattern of genes expressed in a few animals
  (as on the Iso-Seq reference, step 04); they are not treatment effects.
- **Genes whose fold change moves.** 30 gene-contrast pairs of the 3' extension (StringTie)
  shift by more than one log2 unit from the record, 24 of them on extended genes (only one
  where the extension overlaps a gene on the other strand, which the unstranded StringTie
  count could pick up); some gained many reads (LOC134683242,
  Foot DO: mean count 147 to 762), others barely changed in count, so StringTie's
  redistribution of coverage among the lengthened transcripts may be part of it. They should
  be inspected before `ext3` is adopted.
