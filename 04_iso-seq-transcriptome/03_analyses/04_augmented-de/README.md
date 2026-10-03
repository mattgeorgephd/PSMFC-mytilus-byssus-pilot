# 04_augmented-de

Written by `../../01_code/04_augmented_de_comparison.Rmd` (option B). Option B was adopted on
2026-10-02: the count matrix of record is now featureCounts on the `ext3` annotation (taken
by `05_sequence-alignment` step 07). This step keeps the evidence for that choice: the six
treatment-control contrasts of `06` fitted on each of `05_sequence-alignment` step 06's six count matrices (two
counters x three annotations) and on the previous record (StringTie + prepDE on the original
HPC alignment, `05_sequence-alignment/03_analyses/prepDE/`), all with `06`'s rules (same
samples, `~ treatment`, fit on every gene, at least 10 counts in a third of the samples,
apeglm, DEG at padj < 0.05, the 331 mitochondrial loci left out).

Two checks run every time: the refit of the previous matrix reproduces the previous record
(same genes tested and same DEGs in all six contrasts as `06` before the change; fold changes
within 1e-6, checked 2026-10-02), and the fit of featureCounts `ext3` reproduces `06`'s DEG lists
(same genes, same DEGs, fold changes within 1e-8; the step stops otherwise).

| File | Contents |
|---|---|
| `<stringtie|featurecounts>_<refseq|ext3|full>_apeglm.csv.gz`, `previous_record_apeglm.csv.gz` | every tested gene of every contrast: `code`, `feature` (gene key, or `novel:...` in `full`), `baseMean`, apeglm `log2FoldChange`, `lfcSE`, Wald `pvalue`, `padj` |
| `deg_summary.csv` | per matrix and contrast: genes tested, DEGs (up, down), the previous record's DEGs, shared, shared in the same direction, previous DEGs lost (and not tested here), new DEGs (and on novel loci), Spearman of fold changes over genes tested in both |
| `annotation_effect.csv` | per counter and contrast: `ext3` and `full` against that counter's RefSeq control (DEGs shared, gained, gained on novel loci, lost, fold-change Spearman) |
| `record_change.csv` | every gene that is a DEG in the previous record or in the record, per contrast: `status` (DEG in both, previous record only, record only), whether `ext3` extended it, and both mean counts, fold changes and padj |
| `control_vs_previous.csv` | the StringTie + prepDE RefSeq control against the previous record: per-library totals, mitochondrial share per library, day-3 gene totals |
| `mitochondrial_share.csv` | per matrix: median mitochondrial share of the libraries (%), the previous record's, and the Spearman correlation of the per-library shares with the previous record's |
| `counts_per_gene.csv` | per gene, the total count over the day-3 libraries in the previous record and on each matrix |
| `new_DEGs.csv` | per matrix, the DEGs that are not DEGs in the previous record: whether on a novel locus, whether `ext3` extended the gene, whether the previous record tested it, both mean counts, fold changes and padj |
| `byssal_genes.csv` | the byssal structural genes (best BLAST hit named as in `09`) that are a DEG in any matrix or in the previous record, by contrast and source |
| `FIG_deg_counts.png` | DEGs per contrast: the previous record and the three annotations, per counter |
| `FIG_lfc_vs_previous.png` | fold change on each featureCounts annotation against the previous record, genes tested in both (the 3' extension row is the record) |
| `RUN_provenance.txt` | settings, input MD5s |

## Result (2026-10-02)

DEGs per contrast (apeglm padj < 0.05); the record is featureCounts on `ext3`:

| contrast | previous record | StringTie + prepDE: RefSeq control / 3' extension / full | featureCounts: RefSeq / **3' extension (record)** / full |
|---|---|---|---|
| Foot OA | 75 | 83 / 93 / 118 | 163 / **161** / 170 |
| Foot OW | 165 | 180 / 202 / 208 | 208 / **211** / 208 |
| Foot DO | 363 | 377 / 417 / 489 | 636 / **678** / 705 |
| Gill OA | 423 | 410 / 431 / 528 | 588 / **584** / 596 |
| Gill OW | 180 | 175 / 198 / 228 | 227 / **227** / 222 |
| Gill DO | 310 | 299 / 306 / 340 | 398 / **404** / 441 |

- **What the adoption changed** (`record_change.csv`). The record keeps 84, 70, 90, 84, 78 and
  84% of the previous record's DEGs (Foot OA, OW, DO, Gill OA, OW, DO: 63 of 75, 116 of 165,
  325 of 363, 356 of 423, 140 of 180, 259 of 310), with the same sign in all but two of the
  shared ones (both in Gill OA), and adds 98, 95, 353, 228, 87 and 145. The new DEGs were mostly
  near misses in the previous record (median previous padj 0.09 to 0.11; 73 to 83% below 0.2),
  and 54 to 61% of them are genes `ext3` extended; the lost ones are near misses
  in the record (median padj 0.07 to 0.14). Fold changes over the genes tested in both correlate
  at 0.92 to 0.94 (Spearman). 158 of the 50,149 gene-contrast pairs tested in both move by more
  than one log2 unit, and in 154 of them one of the two apeglm estimates is near 0 (|log2 fold
  change| < 0.2): apeglm's shrunk estimate of a gene driven by a few animals flips between near
  0 and the full estimate when the counts change a little. 20 of them are DEGs in both, 19 in
  Gill OA and Gill DO, whose previous estimate was shrunk to about -0.1 and is now -1.2 to -3.7
  (the Wald test, which calls the DEGs, does not shrink).
- **The realignment alone (StringTie control against the previous record).** Per-library
  totals are 0.992 of the previous record's (0.988 to 1.000), the mitochondrial share per
  library ranks almost identically (Spearman 0.992), and day-3 gene totals correlate at 0.997.
  The control keeps 91 to 99% of the previous DEGs per contrast (fold-change Spearman 0.993 to
  0.996) and adds 6 to 25. About 5% of the DEG calls therefore move with the alignment alone;
  differences smaller than that say nothing about the annotation or the counter.
- **The counter matters more than the annotation.** On the same alignment and the RefSeq
  annotation, featureCounts finds 16 to 96% more DEGs than StringTie + prepDE. Its counts are
  about 1.16 times higher on a typical gene (prepDE estimates read bases over 75, while the reads
  average about 63 bases) and its gene-wise dispersions about 9% lower (Foot DO: median 0.065
  against 0.072), so it has more power; it also counts on the sense strand only and leaves out
  reads aligned to several places. Within featureCounts, `ext3` gains 9 to 53 DEGs and loses 11
  to 39 per contrast (fold-change Spearman 0.986 to 0.991 with RefSeq); `full` gains 31 to 119
  and loses 24 to 68.
- **3' extension and full models under StringTie + prepDE.** `ext3` against its control: 12 to
  49 DEGs gained and 2 to 16 lost; `full`: 50 to 178 gained (7 to 37 on novel loci, which have
  no names or GO terms), 15 to 62 lost.
- **Genes whose fold change moves with `ext3`.** Under featureCounts, 6 gene-contrast pairs
  shift by more than one log2 unit between RefSeq and `ext3`. Four have the same Wald padj and
  mean count on both (apeglm's shrunk estimate flipping between near 0 and the full estimate
  for genes driven by a few animals, e.g. LOC134691537 in Gill DO, -2.78 and -0.20); two move
  with the extension: LOC134722878 in Foot OA (2.7 times the reads on `ext3`, most of them in a
  108-base extension past the RefSeq end; -1.27 and padj 0.037 on RefSeq, -0.12 and 0.14 on
  `ext3`) and LOC134717707 in Gill OW. Under StringTie + prepDE more pairs move (30 against the
  previous record), partly because StringTie redistributes coverage among lengthened
  transcripts.
- **Mitochondrial reads.** With StringTie + prepDE the mitochondrial share of each library
  matches the previous record (Spearman 0.985 to 0.992). featureCounts keeps about
  a third of them (median share 1.0% against 2.3%; Spearman 0.85), because reads of the
  mitochondrial genes also align to their nuclear copies and featureCounts leaves out reads
  aligned to several places. The genome matrix's mitochondrial rows are therefore not what
  `06` step 13 tests (`06` README).
- **Byssal genes** (`09`'s definition): the three ACDC genes are DEGs in Gill OA and Gill DO
  in the previous record and in every matrix. Foot protein 11 (LOC134706002) is a DEG in Foot DO
  on every featureCounts matrix, the record included, with an apeglm fold change of -0.05: its
  reads come from three animals (3,760 and 1,688 in two controls, 1,680 in one DO animal; under
  70 in every other library), the pattern of a gene expressed in a few animals, not a treatment
  effect (as foot proteins 4, 10 and 11 on the Iso-Seq reference, step 03). Foot protein 13 is a
  DEG in Foot DO on featureCounts `full` only.
