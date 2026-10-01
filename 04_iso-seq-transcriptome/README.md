# 04_iso-seq-transcriptome

QC of the PacBio Iso-Seq *M. trossulus* transcriptome (`Mtros-hq_transcripts.fasta`), and
the superseded early attempt that used the isoseq transcriptome (rather than the genome) as
the DE reference.

The final pipeline uses the genome as reference (see `05_sequence-alignment` and
`06_differential-expression`). The isoseq-as-reference work is retained under `_superseded/`.

## Layout

```
04_iso-seq-transcriptome/
├── iso-seq-transcriptome.Rproj
├── 01_code/
│   ├── 01_isoseq_transcriptome_check.Rmd   transcriptome length-distribution QC (runnable)
│   ├── 01_isoseq_transcriptome_check.md    knitted output
│   ├── 01_isoseq_transcriptome_check_files/  its figure (knitr)
│   └── _superseded/isoseq-as-reference/     early kallisto-on-isoseq DE attempt (not final)
├── 02_data/                                 transcriptome FASTA is downloaded; not committed
└── 03_analyses/
    └── _superseded/14-kallisto-ng/          output placeholder of the superseded attempt
```

## Input and runnability

`01_isoseq_transcriptome_check.Rmd` downloads `Mtros-hq_transcripts.fasta` from owl
(`https://owl.fish.washington.edu/halfshell/genomic-databank/`) into `02_data/` (git-ignored)
when it is not already there, then runs locally. It is a single side analysis, so the folder
has no batch runner and the repository-level runner does not call it. The
`_superseded/isoseq-as-reference/` scripts are kept as a record and are not maintained.

## Running the expression pipeline on the Iso-Seq transcriptome

The analysis of record maps reads to the genome. To run the same chain against the Iso-Seq
transcriptome instead, as a sensitivity branch:

1. **Reference.** Collapse redundant isoforms (the transcriptome holds several per gene), then
   build a transcript-to-gene table, either by clustering or by mapping each isoform to its
   genome LOC (`03_blast/01_code/03_isoseq_vs_genome_blast.Rmd` already blasts one against the
   other).
2. **Quantification.** Pseudo-align against the collapsed transcriptome with salmon, using
   `--noLengthCorrection` (meant for 3' tag protocols such as QuantSeq, whose counts do not
   scale with transcript length), and summarise to genes with tximport. Then pass the raw
   counts (`txi$counts`, `countsFromAbundance = "no"`) to `DESeqDataSetFromMatrix()` without a
   length offset, as the tximport vignette recommends for 3' tagged RNA-seq (Soneson, Love and
   Robinson 2015, F1000Research 4:1521). The superseded `01_code/_superseded/isoseq-as-reference/`
   kallisto runs are a starting point.
3. **Differential expression.** Write the gene matrix to
   `04_iso-seq-transcriptome/03_analyses/counts/` and point `06_differential-expression` at it
   through one parameter (the count file of `01_clean_count_matrix.Rmd`), with an output root
   per reference so the genome results are not overwritten.
4. **Annotation and GO.** Annotate the transcriptome genes from
   `03_blast/03_analyses/transcriptome-uniprot/` into the same `gene_annotation.tsv` layout
   `07_enrichment` uses; the enrichment, GO slim and gene-mechanics scripts then run unchanged.
5. **Comparison.** Report, per contrast, the overlap of DEGs between references (through the
   isoform-to-LOC map) and whether the enriched terms agree.

None of this is implemented yet; it is the design if the branch is wanted.
