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

The analysis of record maps reads to the genome. Running the same chain against the Iso-Seq
transcriptome, as a sensitivity branch, is planned but **not implemented**: three inputs are
not in the repository (checked 2026-10-01). Steven Roberts' Iso-Seq notebooks
(https://sr320.github.io/iso/, https://sr320.github.io/Myt-GO/) may hold the details; they
were not reachable from the cloud environment used here.

| needed | where it is | status |
|---|---|---|
| trimmed Tag-seq reads | owl, Grace Leuchtenberger's `Github/byssus-exp-analysis/data/raw-trimmed/*_L099_R1_cmb.trim.fastq.gz` (confirmed 2026-10-01; paths in `01_code/_superseded/isoseq-as-reference/01-kallisto-genome.Rmd`) | not in the repository and not reachable from the cloud environment used here; the trimming command is probably in gannet `panopea/PSMFC-mytilus-byssus-pilot/sbatch_scripts/` |
| `Mtros-hq_transcripts.fasta` | owl (`https://owl.fish.washington.edu/halfshell/genomic-databank/`) | downloadable by `01_isoseq_transcriptome_check.Rmd` where owl is reachable |
| isoform-to-gene table | none | must be built (step 1 below); no script here maps isoforms to genome LOCs. `03_blast/01_code/03_isoseq_vs_genome_blast.Rmd` blasts the byssal foot-protein sequences against each reference separately, not one reference against the other |

The design, once those exist:

1. **Isoform-to-gene table.** Align the isoforms to the genome (GCF_036588685.1) with a
   spliced long-read aligner (minimap2 `-ax splice:hq`; Li 2018, Bioinformatics 34:3094) and
   assign each isoform to the reference gene (LOC) it overlaps in the genome annotation. This
   gives genes comparable with the genome branch. Isoforms that overlap no LOC form their own
   genes (clustered, for example, by shared alignment locus).
2. **Quantification.** Pseudo-align the trimmed reads to the isoforms with salmon
   (Patro et al. 2017, Nature Methods 14:417) using `--noLengthCorrection`, which is meant for
   3' tag protocols whose counts do not scale with transcript length, and summarise to genes
   with tximport. Pass the raw counts (`txi$counts`, `countsFromAbundance = "no"`) to
   `DESeqDataSetFromMatrix()` without a length offset, as the tximport vignette recommends
   for 3' tagged RNA-seq (Soneson, Love and Robinson 2015, F1000Research 4:1521).
3. **Differential expression.** Write the gene matrix to
   `04_iso-seq-transcriptome/03_analyses/counts/` and point `06_differential-expression` at it
   through one parameter (the count file of `01_clean_count_matrix.Rmd`), with an output root
   per reference so the genome results are not overwritten. The mitochondrial transcripts
   are separated as in the genome branch (`tools/mt_encoded.R`).
4. **Annotation and GO.** Annotate the genes from the isoform BLAST already here
   (`03_blast/03_analyses/transcriptome-uniprot/Mtros-hq-uniprot_blastx.tab`, best hit per
   gene) into the `gene_annotation.tsv` layout `07_enrichment` uses; the enrichment, GO slim
   and gene-mechanics scripts then run unchanged.
5. **Comparison.** Report, per contrast, the overlap of DEGs between references (through the
   isoform-to-LOC table) and whether the enriched terms agree.

The superseded kallisto runs in `01_code/_superseded/isoseq-as-reference/` are a starting
point for step 2 only; they used the full transcriptome without an isoform-to-gene table.
