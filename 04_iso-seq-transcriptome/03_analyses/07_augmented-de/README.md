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
| `counts_per_gene.csv` | per gene, the total count over the day-3 libraries in the record and on each matrix |
| `new_DEGs.csv` | DEGs of `ext3` and `full` (StringTie + prepDE) that are not DEGs in the record: whether on a novel locus, whether `ext3` extended the gene, whether the record tested it, both mean counts, fold changes and padj |
| `byssal_genes.csv` | the byssal structural genes (best BLAST hit named as in `09`) that are a DEG in any matrix or in the record, by contrast and source |
| `FIG_deg_counts.png` | DEGs per contrast: the record and the three annotations, per counter |
| `FIG_lfc_vs_record.png` | fold change on each annotation (StringTie + prepDE) against the record, genes tested in both |
| `RUN_provenance.txt` | settings, input MD5s |
