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
repository-level `00_run_pipeline.Rmd` do it). By default it runs step 02 only.

| step | script | what it does |
|---|---|---|
| 01 | `01_hisat_stringtie.Rmd` | HISAT2 alignment and StringTie quantification (`-e -B` against the reference annotation); writes `03_analyses/hisat/` and its `sample_list.txt`. HPC only |
| 02 | `02_prepDE.Rmd` | gene and transcript count matrices in `03_analyses/prepDE/` (read length 75) |

`02_prepDE.Rmd` counts from the per-sample StringTie GTFs when step 01's outputs are present,
with `01_code/_prepde.R`, an R port of StringTie's `prepDE.py3` checked to give byte-identical
matrices. Without them, which is the case in this repository, it rebuilds the gene matrix from
the committed transcript matrix (the original `prepDE.py` output) through the reference
transcript table (`03_analyses/hisat/t_data.ctab`) and `02_data/strg_gene_ids.csv`, and checks
the sums. `01_code/_derive_strg_gene_ids.R` is the one-off record of how that ID map was
recovered.

## Layout

```
05_sequence-alignment/
├── sequence-alignment.Rproj
├── 01_code/
│   ├── 00_run_sequence_alignment.Rmd   batch runner
│   ├── 01_hisat_stringtie.Rmd          HISAT2 + StringTie (HPC record)
│   ├── 02_prepDE.Rmd                   count matrices
│   ├── _prepde.R                       R port of prepDE.py3 (sourced by 02)
│   ├── _derive_strg_gene_ids.R         one-off: recovers the StringTie gene IDs
│   └── _superseded/
│       ├── 07-HiSat_GL.Rmd             earlier HISAT2 attempt (different assembly + augustus)
│       └── 07-kallisto.Rmd(.md)        kallisto pseudo-alignment, superseded by HISAT2
├── 02_data/
│   ├── sample-submission/              Tag-seq sequencing submission paperwork
│   └── strg_gene_ids.csv               StringTie gene IDs of 284 reference transcripts
└── 03_analyses/
    ├── hisat/                          reference StringTie tables + MultiQC alignment reports
    ├── prepDE/                         transcript and gene count matrices (the DE input)
    ├── fastqc/{trimmed,untrimmed}/     FastQC per-sample read QC
    ├── _superseded/kallisto/           kallisto quant per sample
    └── knit_html/                      runner reports and logs (git-ignored)
```

## External inputs (not in repo)

| Input | Location |
|-------|----------|
| Raw (untrimmed) Tag-seq reads | gannet: https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/20220405-tagseq/ |
| Trimmed reads (`*_L099_R1_cmb.trim.fastq.gz`), read by `01_hisat_stringtie.Rmd` | produced on the HPC; the trimming step is not in this repository and `20220405-tagseq/` holds raw reads only |
| Genome assembly + annotation | NCBI `GCF_036588685.1` (downloaded by the script into `02_data/ncbi_dataset/`, git-ignored) |

`01_hisat_stringtie.Rmd` and `_superseded/07-HiSat_GL.Rmd` invoke HISAT2 / StringTie at fixed
`/home/shared/...` paths and read inputs not committed here, so they do not run as they are
outside that HPC environment. Packages for step 02: base R, here, rmarkdown.
