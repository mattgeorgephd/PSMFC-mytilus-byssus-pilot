# 04_iso-seq-transcriptome

A sensitivity branch: the differential expression of record (`06`, reads aligned to the
genome) repeated with the Tag-seq reads quantified against the PacBio Iso-Seq *M. trossulus*
transcriptome instead, and the two compared gene by gene. It asks whether the DEG results
depend on the reference. Nothing of record rests on it.

## How to run

Open `iso-seq-transcriptome.Rproj` and knit `01_code/00_run_isoseq.Rmd`. With `online: false`
(the default, and what `00_run_pipeline.Rmd` does as its stage 04) only step 04 runs, from the
committed isoform map and gene counts. `online: true` also rebuilds steps 02 and 03, which
download the transcriptome, the genome and the reads and need minimap2 and salmon
(`minimap2`, `salmon` parameters; minimap2 2.31 and salmon 1.10.3 from bioconda were used).

| step | script | writes to `03_analyses/` | needs |
|---|---|---|---|
| 01 | `01_isoseq_transcriptome_check.Rmd` | nothing (a length QC; its chunks are not evaluated by default) | owl |
| 02 | `02_isoform_gene_map.Rmd` | `02_isoform-gene-map/`: each isoform's genome gene | owl, NCBI, minimap2; about 1.5 hours and 9 GB of memory |
| 03 | `03_salmon_quant.Rmd` | `03_salmon/`: mapping summary and gene counts | gannet, salmon; about 15 minutes for the index plus 3 to 4 minutes per library |
| 04 | `04_isoseq_de_comparison.Rmd` | `04_isoseq-de/`: the TC contrasts on the Iso-Seq counts and the comparison with the genome | committed files only |

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

## Results

To be filled after the run.

## Inputs

| Input | Source |
|---|---|
| Iso-Seq transcriptome, `Mtros-hq_transcripts.fasta` | owl, https://owl.fish.washington.edu/halfshell/genomic-databank/ |
| Genome and annotation, `GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.fna.gz`, `..._genomic.gff.gz` | NCBI, https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/588/685/GCF_036588685.1_PNRI_Mtr1.1.1.hap1/ (MD5s checked against NCBI's) |
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
│   ├── _paths.R                              shared paths (sourced by every step)
│   └── _superseded/                          retired scripts (README inside)
├── 02_data/                                  downloads, git-ignored (README)
└── 03_analyses/
    ├── 02_isoform-gene-map/                  step 02 (README)
    ├── 03_salmon/                            step 03 (README)
    ├── 04_isoseq-de/                         step 04 (README)
    └── _superseded/                          retired outputs (README)
```
