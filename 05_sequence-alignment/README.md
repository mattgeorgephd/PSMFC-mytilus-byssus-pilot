# 05_sequence-alignment

Read QC, alignment and quantification of the Tag-seq reads against the *Mytilus trossulus*
genome (GenBank GCA_036588685.1 / RefSeq GCF_036588685.1), ending in the count matrices that
`06_differential-expression` starts from.

The alignment itself, HISAT2 + StringTie (`01_code/01_hisat_stringtie.Rmd`), ran on a
collaborator HPC workstation (paths under `/home/shared/...`) against trimmed reads and a genome
that are not stored in this repository; it is kept verbatim as the method record, with its
summary outputs committed in `03_analyses/`. The count-matrix step (`02_prepDE.Rmd`) runs
anywhere.

## How to run

Open `sequence-alignment.Rproj` and knit `01_code/00_run_sequence_alignment.Rmd` (or let the
repository-level `00_run_pipeline.Rmd` do it). By default it runs steps 02 and 03.

| step | script | what it does |
|---|---|---|
| 01 | `01_hisat_stringtie.Rmd` | HISAT2 alignment and StringTie quantification (`-e -B` against the reference annotation); writes `03_analyses/hisat/` and its `sample_list.txt`. HPC only |
| 02 | `02_prepDE.Rmd` | gene and transcript count matrices in `03_analyses/prepDE/` (read length 75) |
| 03 | `03_read_trimming.Rmd` | how the trimmed reads were made: per-library read retention from the committed FastQC tables, and (with `online: true`) the trimming recipe reproduced on a sample of reads; writes `03_analyses/read_trimming/` |

`02_prepDE.Rmd` counts from the per-sample StringTie GTFs when step 01's outputs are present,
with `01_code/_prepde.R`, an R port of StringTie's `prepDE.py3` checked to give byte-identical
matrices. Without them, which is the case in this repository, it rebuilds the gene matrix from
the committed transcript matrix (the original `prepDE.py` output) through the reference
transcript table (`03_analyses/hisat/t_data.ctab`) and `02_data/strg_gene_ids.csv`, and checks
the sums. `01_code/_derive_strg_gene_ids.R` is the one-off record of how that ID map was
recovered.

## The reads

GSAF (University of Texas at Austin, job JA22078) sequenced the 131 libraries as 3' Tag-seq,
single-end 100 bp on a NovaSeq S1, over two lanes. The reads aligned in step 01 are the files
`*_L099_R1_cmb.trim.fastq.gz`: one per library, the two lanes combined (`cmb`, lane `L099`),
trimmed. Who trimmed them was not recorded, but step 03 reproduces the recipe: 99.8% of the
trimmed reads in a sample of library T001F (32,296 of 32,349) are identical, base for base and
in quality, to the output of

    tagseq_clipper.pl raw.fastq \
      | fastx_clipper -a AAAAAAAA -l 20 -Q33 \
      | fastx_clipper -a AGATCGGAAG -l 20 -Q33

that is, the Matz-lab Tag-seq clipper (https://github.com/z0on/tag-based_RNAseq; Meyer,
Aglyamova and Matz 2011, *Molecular Ecology* 20:3599), which cuts the 5' leader (degenerate
bases and the template-switching G run) and drops PCR duplicates (reads sharing the leader and
the first 20 bases), then FASTX-Toolkit 0.0.14 `fastx_clipper` for the poly-A tail and the
Illumina adapter, keeping reads of 20 bases or more, with no quality filter. The other recipes
tried match far fewer reads (`03_analyses/read_trimming/recipe_check.csv`). As a result the
analysed reads are deduplicated: in the 72 libraries whose raw reads are on gannet with both
lanes, a median of 43% of raw reads remain (35 to 58%), and the 131 trimmed libraries hold 0.78
to 3.86 million reads (median 2.72 million) of 20 to 95 bases (mean about 63).

Two consequences for anyone rerunning from raw reads:

- **Raw reads exist on gannet for 73 of the 131 libraries only** (`20220405-tagseq/`: T001 to
  T030, T131G and T132 to T137, and T030G has lane 1 only). The raw reads of the other 58
  (T031 to T058 and T110 to T131F) have not been located; their trimmed files are the earliest
  copy known.
- **The August 2022 trimming is a different one.** Matt George's 2022 script
  (`01_code/_superseded/1_1_process-tagseq-data-mytilus.Rmd`) trimmed the first 73 libraries with
  cutadapt (a fixed 15-base 5' cut, no deduplication) and aligned them to other *Mytilus*
  genomes; `03_analyses/fastqc/multiqc_report_trimmed_merged_2022-08-09.html` is the QC of those
  reads, not of the reads analysed here.

## Layout

```
05_sequence-alignment/
├── sequence-alignment.Rproj
├── 01_code/
│   ├── 00_run_sequence_alignment.Rmd   batch runner
│   ├── 01_hisat_stringtie.Rmd          HISAT2 + StringTie (HPC record)
│   ├── 02_prepDE.Rmd                   count matrices
│   ├── 03_read_trimming.Rmd            how the trimmed reads were made (retention, recipe check)
│   ├── _prepde.R                       R port of prepDE.py3 (sourced by 02)
│   ├── _derive_strg_gene_ids.R         one-off: recovers the StringTie gene IDs
│   └── _superseded/                    records (README inside)
│       ├── 1_1_process-tagseq-data-mytilus.Rmd   2022 cutadapt trimming + HISAT2 to other genomes
│       ├── 07-HiSat_GL.Rmd             earlier HISAT2 attempt (different assembly + augustus)
│       └── 07-kallisto.Rmd(.md)        kallisto pseudo-alignment, superseded by HISAT2
├── 02_data/
│   ├── sample-submission/              Tag-seq sequencing submission paperwork
│   └── strg_gene_ids.csv               StringTie gene IDs of 284 reference transcripts
└── 03_analyses/
    ├── hisat/                          reference StringTie tables + MultiQC alignment reports
    ├── prepDE/                         transcript and gene count matrices (the DE input)
    ├── fastqc/{trimmed,untrimmed}/     FastQC per-sample read QC
    ├── read_trimming/                  read retention and the trimming recipe check
    ├── _superseded/kallisto/           kallisto quant per sample
    └── knit_html/                      runner reports and logs (git-ignored)
```

## External inputs (not in repo)

| Input | Location |
|-------|----------|
| Raw (untrimmed) Tag-seq reads, 73 of the 131 libraries | gannet: https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/20220405-tagseq/ |
| Trimmed reads (`*_L099_R1_cmb.trim.fastq.gz`), all 131 libraries, read by `01_hisat_stringtie.Rmd` and by `04_iso-seq-transcriptome` | gannet: https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/byssus-exp-analysis/data/raw-trimmed/ (recipe above) |
| Genome assembly + annotation | NCBI `GCF_036588685.1` (downloaded by the script into `02_data/ncbi_dataset/`, git-ignored) |

`01_hisat_stringtie.Rmd` and `_superseded/07-HiSat_GL.Rmd` invoke HISAT2 / StringTie at fixed
`/home/shared/...` paths and read inputs not committed here, so they do not run as they are
outside that HPC environment. Packages for steps 02 and 03: base R, here, rmarkdown; step 03
with `online: true` also needs curl, perl, FASTX-Toolkit 0.0.14 and cutadapt (5.2 used).
