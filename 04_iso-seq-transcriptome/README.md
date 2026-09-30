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
│   ├── 05-IsoSeq-transcriptome-check.Rmd   transcriptome length-distribution QC (runnable)
│   ├── 05-IsoSeq-transcriptome-check.md     knitted output + figure
│   └── _superseded/isoseq-as-reference/     early kallisto-on-isoseq DE attempt (not final)
├── 02_data/                                 transcriptome FASTA is downloaded; not committed
└── 03_analyses/
    └── _superseded/14-kallisto-ng/          output placeholder of the superseded attempt
```

## Input and runnability

`05-IsoSeq-transcriptome-check.Rmd` downloads `Mtros-hq_transcripts.fasta` from owl
(`https://owl.fish.washington.edu/halfshell/genomic-databank/`) into `02_data/` when it is
not already there, then runs locally. The
`_superseded/isoseq-as-reference/` scripts are kept as a record and are not maintained.
