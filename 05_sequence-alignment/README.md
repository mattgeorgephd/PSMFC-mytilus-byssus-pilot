# 05_sequence-alignment

Read QC and alignment of Tag-seq reads to the *Mytilus trossulus* genome (GenBank
GCA_036588685.1 / RefSeq GCF_036588685.1), producing the count inputs used downstream by
`06_differential-expression`.

The authoritative pipeline is HISAT2 + StringTie (`01_code/13-Hisat.Rmd`). These scripts
ran on a collaborator HPC workstation (paths under `/home/shared/...`) against raw reads
and a genome that are NOT stored in this repo; they are kept verbatim as the method record.
Only their committed outputs are in `03_analyses/`.

## Layout

```
05_sequence-alignment/
├── sequence-alignment.Rproj
├── 01_code/
│   ├── 13-Hisat.Rmd                 authoritative HISAT2 + StringTie (HPC)
│   └── _superseded/
│       ├── 07-HiSat_GL.Rmd          earlier HISAT2 attempt (different assembly + augustus)
│       └── 07-kallisto.Rmd(.md)     kallisto pseudo-alignment, superseded by HISAT2
├── 02_data/
│   └── sample-submission/           Tag-seq sequencing submission paperwork
└── 03_analyses/
    ├── hisat/                       StringTie ctabs + MultiQC alignment reports
    ├── fastqc/{trimmed,untrimmed}/  FastQC per-sample read QC
    └── _superseded/kallisto/        kallisto quant per sample
```

## External inputs (not in repo)

| Input | Location |
|-------|----------|
| Raw (untrimmed) Tag-seq reads | gannet: https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/20220405-tagseq/ |
| Trimmed reads (`*_L099_R1_cmb.trim.fastq.gz`), read by `13-Hisat.Rmd` | produced on the HPC; the trimming step is not in this repository and `20220405-tagseq/` holds raw reads only |
| Genome assembly + annotation | NCBI `GCF_036588685.1` (downloaded in the script) |

`13-Hisat.Rmd` and `07-HiSat_GL.Rmd` invoke HISAT2/StringTie at fixed `/home/shared/...`
paths and read inputs not committed here, so they do not run as-is outside that HPC
environment. They are retained as the documented method.
