# 04_iso-seq-transcriptome

A sensitivity branch: the differential expression of record (`06`, reads aligned to the
genome) repeated with the Tag-seq reads quantified against the PacBio Iso-Seq *M. trossulus*
transcriptome instead, and the two compared gene by gene (steps 01 to 04). It asks whether the
DEG results depend on the reference. Steps 05 to 07 are option B, a parallel analysis: the
reads recounted on the genome with the RefSeq annotation augmented by the Iso-Seq isoforms,
built so that it could become the method of record if its results justify it. Nothing of
record rests on this folder.

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
| 01 | `01_isoseq_transcriptome_check.Rmd` | nothing (a length QC; its chunks are not evaluated by default) | owl |
| 02 | `02_isoform_gene_map.Rmd` | `02_isoform-gene-map/`: each isoform's genome gene | owl, NCBI, minimap2; about 1.5 hours and 9 GB of memory |
| 03 | `03_salmon_quant.Rmd` | `03_salmon/`: mapping summary and gene counts | gannet, salmon; about 15 minutes for the index plus 3 to 4 minutes per library |
| 04 | `04_isoseq_de_comparison.Rmd` | `04_isoseq-de/`: the TC contrasts on the Iso-Seq counts and the comparison with the genome | committed files only |
| 05 | `05_augmented_annotation.Rmd` | `05_augmented-annotation/`: option B's annotations (RefSeq; 3' extension; full models) | step 02's alignment; a few minutes |
| 06 | `06_genome_recount.Rmd` | `06_genome-recount/`: the reads realigned (HISAT2) and counted on each annotation (StringTie + prepDE, featureCounts) | gannet, NCBI, the tools above; 15 minutes for the index plus about 4 minutes per library |
| 07 | `07_augmented_de_comparison.Rmd` | `07_augmented-de/`: the TC contrasts on each recount, compared with the record | committed files only |

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

**Option B (steps 05 to 07), a parallel analysis.** Steps 02 to 04 showed that the RefSeq gene
models often end before the 3' ends that Tag-seq reads. Option B keeps the genome as the
reference and changes only the annotation the reads are counted on, so its genes keep their
keys, names and GO annotation and the rest of the pipeline could run on it unchanged.

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

**Differential expression on the two references (step 04).**

| contrast | DEGs, Iso-Seq / genome | genome DEGs that are also Iso-Seq DEGs | fold-change Spearman |
|---|---|---|---|
| Foot OA | 194 / 75 | 58 (77%) | 0.84 |
| Foot OW | 262 / 165 | 89 (54%) | 0.82 |
| Foot DO | 867 / 363 | 291 (80%) | 0.87 |
| Gill OA | 884 / 423 | 304 (72%) | 0.81 |
| Gill OW | 366 / 180 | 123 (68%) | 0.81 |
| Gill DO | 600 / 310 | 227 (73%) | 0.82 |

**What this means for the results of record.**

- **The genome results hold.** Most genome DEGs are DEGs on the Iso-Seq reference too (54 to
  80% per contrast), all in the same direction but three in Gill OA, and fold changes over the
  genes tested on both correlate at 0.81 to 0.87. The genome DEGs that are not Iso-Seq DEGs
  are mostly just short of the threshold there (median Iso-Seq padj 0.09 to 0.14 per
  contrast; 67 to 84% below 0.2), and 359 of those 368 keep their direction.
- **The mitochondrial signal does not depend on the reference.** Each library's mitochondrial
  share of reads correlates between the references at Spearman 0.87 (day-3 foot 0.85, day-3
  gill 0.81).
- **Byssal genes.** Named as `09` names the byssal structural genes (best BLAST hit a foot
  protein, plaque protein, precollagen, ACDC or byssal tyrosinase; 72 genes), the genome
  analysis has six byssal DEGs: the three ACDC genes, down about two- to four-fold in both
  Gill OA and Gill DO. All six are DEGs on the Iso-Seq reference too, with the same fold
  changes. The Iso-Seq reference adds foot proteins 4, 10 and 11 in Foot OA and Foot DO, with
  shrunk fold changes near 0 and a small padj, the sign of a gene expressed in a few animals
  only; they should not be read as treatment effects. (An earlier version of this README gave
  "Foot OA 1 of 1, Foot OW 1 of 1, Foot DO 7 of 8" from a broader, unrecorded definition; it is
  replaced by this count, which step 07 reproduces for option B.)
- **The genome analysis is conservative.** The Iso-Seq reference finds about twice as many DEGs.
  On genes counted by both it places about 1.3 times the reads on a typical gene (2.4 times or
  more on a quarter of genes), and most Iso-Seq-only DEGs on annotated genes are genes the
  genome analysis tested with fewer reads and a near-miss padj (median 0.15 to 0.23). The
  extra reads come largely from 3' ends beyond the RefSeq gene models (18% of assigned reads)
  and from unannotated loci (14.5%), which the genome analysis, counting annotated exons only,
  cannot see.

**Caveats.** No decoy sequence (reads from loci the Iso-Seq set lacks may be placed on similar
isoforms); salmon shares reads among near-duplicate isoforms by EM, so gene sums can borrow
from paralogs; novel loci have no names or GO terms; the genome counts are coverage-based
prepDE estimates rather than read counts. Whether to carry the Iso-Seq gene models into the
genome analysis (extend the annotation with them and recount) is an open decision in
`../tasks.md`; this branch stays a sensitivity analysis until then.

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
