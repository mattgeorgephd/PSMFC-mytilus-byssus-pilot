# 04_iso-seq-transcriptome

Two things. Steps 01 to 04 are a sensitivity branch: the differential expression of record
(`06`, reads aligned to the genome) repeated with the Tag-seq reads quantified against the
PacBio Iso-Seq *M. trossulus* transcriptome instead, and the two compared gene by gene; it asks
whether the DEG results depend on the reference. Steps 05 to 07 are option B: the reads
realigned to the genome and recounted on the RefSeq annotation augmented by the Iso-Seq
isoforms, built as a parallel analysis and adopted on 2026-10-02 in its featureCounts `ext3`
form. **The count matrix of record comes from here**: step 06's
`featurecounts_ext3_gene_counts.csv.gz`, which `05_sequence-alignment` step 04 passes to `06`.
Step 07 keeps the comparison with the previous record (StringTie + prepDE on the HPC
alignment) that the choice rested on.

## How to run

Open `iso-seq-transcriptome.Rproj` and knit `01_code/00_run_isoseq.Rmd`. With `online: false`
(the default, and what `00_run_pipeline.Rmd` does as its stage 04) steps 04 and 07 run, from
the committed gene counts. `online: true` also rebuilds steps 02, 03, 05 and 06, which download
the transcriptome, the genome, its GTF and the reads and need minimap2, salmon, HISAT2,
StringTie, Subread and samtools (`minimap2`, `salmon` and `bin` parameters; minimap2 2.31,
salmon 1.10.3, HISAT2 2.2.1, StringTie 2.2.1, Subread 2.1.1 and samtools 1.24 from bioconda
were used).

| step | script | writes to `03_analyses/` | needs |
|---|---|---|---|
| 01 | `01_isoseq_transcriptome_check.Rmd` | nothing in `03_analyses/` (a length QC knit by hand, not by the runner; it downloads the FASTA if missing and rewrites its own `01_isoseq_transcriptome_check.md`) | owl |
| 02 | `02_isoform_gene_map.Rmd` | `02_isoform-gene-map/`: each isoform's genome gene | owl, NCBI, minimap2; about 1.5 hours and 9 GB of memory |
| 03 | `03_salmon_quant.Rmd` | `03_salmon/`: mapping summary and gene counts | gannet, salmon; about 15 minutes for the index plus 3 to 4 minutes per library |
| 04 | `04_isoseq_de_comparison.Rmd` | `04_isoseq-de/`: the TC contrasts on the Iso-Seq counts and the comparison with the genome | committed files only |
| 05 | `05_augmented_annotation.Rmd` | `05_augmented-annotation/`: option B's annotations (RefSeq; 3' extension; full models) | step 02's alignment; a few minutes |
| 06 | `06_genome_recount.Rmd` | `06_genome-recount/`: the reads realigned (HISAT2) and counted on each annotation (StringTie + prepDE, featureCounts) | gannet, NCBI, the tools above; 15 minutes for the index plus about 4 minutes per library |
| 07 | `07_augmented_de_comparison.Rmd` | `07_augmented-de/`: the TC contrasts on each recount and on the previous record, compared; checks that the record (featureCounts `ext3`) reproduces `06` | committed files only |

## Design

1. **Isoforms to genes (step 02).** The transcriptome's 411,251 isoforms were clustered but
   never collapsed to genes. Each is aligned to the genome (RefSeq GCF_036588685.1) with
   minimap2's spliced preset for accurate full-length transcripts (`-x splice:hq -uf`; Li
   2018, *Bioinformatics* 34:3094) and takes the annotated gene (annotation release
   RS_2024_02, the one `05` counts against) whose exons it overlaps most on its strand.
   Isoforms aligned well but outside every gene form novel loci; unaligned ones stay single
   features. Isoforms on mitochondrial loci are flagged.
2. **Quantification (step 03).** The 131 trimmed read files the genome branch aligned are
   pseudo-aligned to the isoforms with salmon in selective-alignment mode (Patro et al. 2017,
   *Nature Methods* 14:417) with `--noLengthCorrection`, because 3' Tag-seq yields about one
   read per transcript whatever its length. tximport sums the isoforms of each gene
   (`countsFromAbundance = "no"`; Soneson, Love and Robinson 2015, *F1000Research* 4:1521),
   and the summed read counts go to DESeq2 without a length offset, as the tximport vignette
   advises for 3' tagged RNA-seq.
3. **Differential expression and comparison (step 04).** The six TC contrasts of `06`, with
   its samples, model, filter, order of operations and apeglm shrinkage; the mitochondrial
   loci left out as in `06`. Genes are matched to the genome results through `gene_key()`.
   Given the genome counts instead of the Iso-Seq counts, step 04 reproduces `06` exactly
   (checked 2026-10-01: same genes tested and same DEGs in all six contrasts).

**Option B (steps 05 to 07), adopted 2026-10-02.** Steps 02 to 04 showed that the RefSeq gene
models often end before the 3' ends that Tag-seq reads. Option B keeps the genome as the
reference and changes only the annotation the reads are counted on, so its genes keep their
keys, names and GO annotation and the rest of the pipeline runs on it unchanged. It was built
beside the record and compared with it (step 07); featureCounts on the `ext3` annotation then
became the count matrix of record.

4. **Annotations (step 05).** Three versions of RS_2024_02: `refseq` unchanged (the control);
   `ext3`, each transcript's last exon extended to the 3' end of same-gene isoforms whose last
   block overlaps it, stopping before the next gene on the strand; and `full`, the isoforms
   added as transcripts of their genes and the novel loci as new genes (a locus overlapping
   one gene's exons joins it). Readthrough isoforms (exons on two genes), ambiguous
   placements and mitochondrial isoforms are not used.
5. **Realignment and counting (step 06).** The record's alignments were not kept, so the 131
   libraries are realigned with HISAT2 2.2.1 and the record's settings. The record's index
   held the RefSeq splice sites and exons, which takes more memory than this machine has; here
   the same splice sites are given at alignment time. Each alignment is counted on the three
   annotations with StringTie + prepDE, as the record counted (Pertea et al. 2015, *Nature
   Biotechnology* 33:290), and with featureCounts (unique alignments, sense strand; Liao,
   Smyth and Shi 2014, *Bioinformatics* 30:923), the usual counter for 3' Tag-seq.
6. **Comparison (step 07).** The six TC contrasts with `06`'s rules on each of the six
   matrices, compared with the record (what adopting a matrix would change) and with the same
   counter's RefSeq control (the effect of the annotation alone). The RefSeq control against
   the record measures what the realignment alone changes. Given the record's own counts,
   step 07 reproduces `06` exactly (checked 2026-10-02).

## Limitations

- **No decoys.** salmon's recommended index adds the genome as decoy sequence so that reads
  from unannotated loci are not forced onto transcripts (Srivastava et al. 2020, *Genome
  Biology* 21:239). The index here has none: a decoy index of this genome needs more memory
  than the 15 GB machine this was run on. Reads from loci the Iso-Seq set lacks may be
  placed on similar isoforms.
- **Different animals.** The transcriptome comes from pooled tissue (`MT_Pool`) whose animals
  are not recorded here, and the genome from another animal; allelic differences between them
  and the experimental animals lower the share of reads each reference can place.
- **Isoform redundancy.** The isoforms were not collapsed, so many are near-duplicates;
  salmon's EM shares reads among them, and step 03 sums them per gene, so gene counts are
  unaffected, but isoform-level estimates are not meaningful and are not used.

## Results (2026-10-02)

**Isoforms to genes (step 02).** 311,399 of the 411,251 isoforms (76%) fall on 19,856 annotated
genes, 62,811 form 16,685 novel loci outside every gene, 2,456 are mitochondrial and 34,585
stay unassigned (20,317 of them do not align to the genome). Where this map and the retired
CDS-based one both put an isoform on a nuclear gene, they agree for 97.0%
(`03_analyses/02_isoform-gene-map/README.md`).

**Quantification (step 03).** salmon assigns a median 66.2% of each library's reads to the
isoforms, about what HISAT2 aligns to the genome (foot 70%, gill 63% on both). Of the assigned
reads, a median 52% fall on isoforms that end within their gene's annotation, 18% on isoforms
of annotated genes that run more than 100 bases past the last annotated exon, 14.5% on novel
loci, 7% on unassigned isoforms and 7% on mitochondrial features. The two libraries removed at
QC are 64% (T051G) and 49% (T051F) mitochondrial, against at most 10% in the others.

**Differential expression on the two references (step 04),** against the genome analysis of
record (featureCounts on the 3'-extended RefSeq annotation, since 2026-10-02):

| contrast | DEGs, Iso-Seq / genome | genome DEGs that are also Iso-Seq DEGs | fold-change Spearman |
|---|---|---|---|
| Foot OA | 194 / 161 | 93 (58%) | 0.87 |
| Foot OW | 262 / 211 | 124 (59%) | 0.85 |
| Foot DO | 867 / 678 | 483 (71%) | 0.89 |
| Gill OA | 884 / 584 | 381 (65%) | 0.84 |
| Gill OW | 366 / 227 | 157 (69%) | 0.84 |
| Gill DO | 600 / 404 | 284 (70%) | 0.85 |

**What this means for the results of record.**

- **The genome results hold.** Most genome DEGs are DEGs on the Iso-Seq reference too (58 to
  71% per contrast), all in the same direction but one in Gill OA, and fold changes over the
  genes tested on both correlate at 0.84 to 0.89. The genome DEGs that are not Iso-Seq DEGs
  are mostly just short of the threshold there (median Iso-Seq padj 0.10 to 0.13 per
  contrast; 65 to 77% below 0.2), and 658 of those 665 keep their direction.
- **The mitochondrial signal does not depend on the reference.** Each library's mitochondrial
  share of reads (`06` step 13, counted on the mitochondrial genome alone) correlates with the
  Iso-Seq reference's at Spearman 0.91 (day-3 foot 0.90, day-3 gill 0.84; 0.87 with the
  previous genome count). No step writes these correlations: they were computed by hand from
  `03_analyses/03_salmon/read_classes_by_library.csv` (`pct_mitochondrial`) and
  `06_differential-expression/03_analyses/mitochondrial/mt_share_by_sample.csv`, and recomputed
  from the committed files on 2026-10-03 (all 129 libraries in both, 0.914; day-3 foot 0.899,
  day-3 gill 0.836).
- **Byssal genes.** Named as `09` names the byssal structural genes (best BLAST hit a foot
  protein, plaque protein, precollagen, ACDC or byssal tyrosinase; 72 genes), the genome
  analysis has the three ACDC genes down about two- to four-fold in both Gill OA and Gill DO,
  DEGs on the Iso-Seq reference too with the same fold changes, and foot protein 11 in Foot DO,
  whose reads come from three animals and whose shrunk fold change is near 0. The Iso-Seq
  reference adds foot proteins 4, 10 and 11 in Foot OA and Foot DO with the same few-animal
  pattern. None of these foot proteins should be read as a treatment effect. (An earlier version
  of this README gave "Foot OA 1 of 1, Foot OW 1 of 1, Foot DO 7 of 8" from a broader,
  unrecorded definition; it is replaced by this count, which step 07 reproduces.)
- **The Iso-Seq reference still finds more DEGs,** 1.2 to 1.6 times as many, and a median 1.25
  times the reads on genes counted by both; most Iso-Seq-only DEGs on annotated genes are genes
  the genome analysis tested with fewer reads and a near-miss padj (median 0.11 to 0.17). Before
  2026-10-02 the gap was about twofold, and most of it was the counter, not the 3' ends: prepDE
  estimated read bases over 75 for reads averaging 63 bases, and its dispersions were higher
  (option B, below).

**Caveats.** No decoy sequence (reads from loci the Iso-Seq set lacks may be placed on similar
isoforms); salmon shares reads among near-duplicate isoforms by EM, so gene sums can borrow
from paralogs; novel loci have no names or GO terms; featureCounts leaves out reads with several
placements, which salmon shares.

**Option B, adopted on 2026-10-02 (steps 05 to 07).**

- **Annotations.** `ext3` extends 18,470 transcripts of 9,781 genes by a median 76 bases (4.6%
  more exonic sequence); `full` adds 303,668 isoforms to 18,908 genes and keeps 12,487 novel
  loci as genes (54% more exonic sequence; `03_analyses/05_augmented-annotation/README.md`).
- **Realignment.** HISAT2 here aligns a median 0.45 percentage points fewer reads than the
  HPC alignment; counted on the unchanged RefSeq annotation with StringTie + prepDE (the
  control), the per-library totals are 0.992 of the previous record's and the DEG calls keep 91
  to 99% of its DEGs (fold-change Spearman 0.993 to 0.996). About 5% of DEG calls move with the
  realignment alone.
- **DEGs** (previous record; StringTie + prepDE control / `ext3` / `full`; featureCounts control /
  **`ext3`, the record** / `full`): Foot OA 75; 83 / 93 / 118; 163 / **161** / 170. Foot OW 165;
  180 / 202 / 208; 208 / **211** / 208. Foot DO 363; 377 / 417 / 489; 636 / **678** / 705. Gill OA
  423; 410 / 431 / 528; 588 / **584** / 596. Gill OW 180; 175 / 198 / 228; 227 / **227** / 222.
  Gill DO 310; 299 / 306 / 340; 398 / **404** / 441.
- **The counter matters more than the annotation.** featureCounts (unique reads, sense strand)
  finds 16 to 96% more DEGs than StringTie + prepDE on the same alignment and RefSeq
  annotation: its counts are about 1.16 times higher (prepDE divides read bases by 75; the reads
  average 63) and its dispersions about 9% lower. Within featureCounts, `ext3` gains 9 to 53 DEGs
  and loses 11 to 39 per contrast (fold-change Spearman 0.986 to 0.991 with RefSeq), and only 6
  gene-contrast pairs move by more than one log2 unit, 2 of them because of the extension; `full`
  changes more (31 to 119 gained, 24 to 68 lost, some on unnamed novel loci).
- **What the adoption changed.** The record keeps 70 to 90% of the previous record's DEGs per
  contrast, all shared ones but two in the same direction, and adds 87 to 353, mostly near
  misses before (`03_analyses/07_augmented-de/record_change.csv`).
- **Mitochondrial genes.** featureCounts keeps only about a third of the mitochondrial reads in
  the genome alignment (reads with several placements are left out, and most mitochondrial
  reads also align to the nuclear copies), so `06` step 13 counts the mitochondrial genes on the
  mitochondrial genome alone (`05_sequence-alignment` step 05).
- **Unchanged by any version:** the byssal DEGs (the three ACDC genes, Gill OA and Gill DO).

## Inputs

| Input | Source |
|---|---|
| Iso-Seq transcriptome, `Mtros-hq_transcripts.fasta` | owl, https://owl.fish.washington.edu/halfshell/genomic-databank/ |
| Genome and annotation, `GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.fna.gz`, `..._genomic.gff.gz`, `..._genomic.gtf.gz` | NCBI, https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/588/685/GCF_036588685.1_PNRI_Mtr1.1.1.hap1/ (MD5s checked against NCBI's) |
| Trimmed Tag-seq reads, `*_L099_R1_cmb.trim.fastq.gz` | gannet, `panopea/PSMFC-mytilus-byssus-pilot/byssus-exp-analysis/data/raw-trimmed/` (how they were made: `05_sequence-alignment/README.md`, "The reads") |
| Genome-branch results and samples | `06_differential-expression/03_analyses/` |

All downloads go to `02_data/` and are git-ignored (`02_data/README.md`). Steven Roberts'
notebooks on this transcriptome (https://sr320.github.io/iso/, https://sr320.github.io/Myt-GO/)
cover the length QC (step 01) and the isoform BLAST annotation (`03_blast/03_analyses/transcriptome-uniprot/`);
neither maps isoforms to genome genes.

## Layout

```
04_iso-seq-transcriptome/
├── iso-seq-transcriptome.Rproj
├── 01_code/
│   ├── 00_run_isoseq.Rmd                     batch runner
│   ├── 01_isoseq_transcriptome_check.Rmd     length QC (+ .md and figure, knitted)
│   ├── 02_isoform_gene_map.Rmd               isoforms to genome genes
│   ├── 03_salmon_quant.Rmd                   salmon + tximport gene counts
│   ├── 04_isoseq_de_comparison.Rmd           TC contrasts and comparison with 06
│   ├── 05_augmented_annotation.Rmd           option B: RefSeq augmented with the isoforms
│   ├── 06_genome_recount.Rmd                 option B: realign and count on each annotation
│   ├── 07_augmented_de_comparison.Rmd        option B: TC contrasts on each recount vs the record
│   ├── _recount_library.sh                   step 06's per-library work (align, count, clean up)
│   ├── _prepde_sample.R                      prepDE's gene counts for one StringTie sample
│   ├── _paths.R                              shared paths (sourced by every step)
│   └── _superseded/                          retired scripts (README inside)
├── 02_data/                                  downloads, git-ignored (README)
└── 03_analyses/
    ├── 02_isoform-gene-map/                  step 02 (README)
    ├── 03_salmon/                            step 03 (README)
    ├── 04_isoseq-de/                         step 04 (README)
    ├── 05_augmented-annotation/              step 05 (README)
    ├── 06_genome-recount/                    step 06 (README)
    ├── 07_augmented-de/                      step 07 (README)
    └── _superseded/                          retired outputs (README)
```
