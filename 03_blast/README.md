# 03_blast

Sequence-similarity searches used to annotate genes and to compare the Iso-Seq
transcriptome to the genome. Two lines of work:

- Genome CDS vs Swiss-Prot + Mytilus-foot proteins (`01_genome_blast.Rmd`), producing the
  gene-to-GO mapping (`LOC_GO_list.txt`, `g.spid.txt`). `LOC_GO_list.txt` is read by
  `04_sequence-alignment` steps 04 and 05 and `05_differential-expression` step 01 (the
  mitochondrial loci, through `tools/mt_encoded.R`), `05_differential-expression` step 06,
  `06_iso-seq-transcriptome` step 04, `07_enrichment` step 01 and
  `09_gene-mechanics-correlation` steps 01 and 04. The search ran on the lab server in June to
  August 2024 against Swiss-Prot release 2024_04 (identified from its hits) and the committed
  foot proteins; the committed tables in `03_analyses/genome-foot/` are its outputs. The 2024
  script did not run as written (it joined an object, `masterID`, that it never made, among
  other faults) and is kept in `01_code/_superseded/`. The fixed script runs: from the 2024
  blastx table (gannet) its R steps rebuild `g.spid.txt` byte for byte and, in its `as_2024`
  mode, `LOC_GO_list.txt` byte for byte (given the 2024 UniProt records, which UniProt no longer
  serves and which were rebuilt from the committed file for the test). With `run: false`, the
  default, it only checks the committed tables against a blastx table given as `blastx_tab`.
  `02_genome_blast_uniprot_check.Rmd`, formerly `08_gene-annotation/01_code/Annotation.Rmd`,
  inspects a SwissProt-only blastx of the same CDS and the UniProt annotation download; it
  writes only HPC intermediates, nothing in this repository.
- Iso-Seq transcriptome vs UniProt (`_uniprot_retrieval.py` + outputs) and Iso-Seq vs genome
  (`03_isoseq_vs_genome_blast.Rmd`).

Grace's BLAST scripts ran on an HPC workstation (`/home/shared/...`) against databases and
FASTAs not stored in this repo. Steps 02 and 03 are kept verbatim as the method record; step 01
was fixed on 2026-10-03 (its 2024 version is in `01_code/_superseded/`). Only committed
outputs are in `03_analyses/`.

## Layout

```
03_blast/
├── blast.Rproj
├── 01_code/
│   ├── 01_genome_blast.Rmd                genome CDS vs Swiss-Prot+foot (HPC, run: true); writes LOC_GO/g.spid
│   ├── 02_genome_blast_uniprot_check.Rmd  checks of a SwissProt-only genome blastx (HPC); writes nothing
│   ├── 03_isoseq_vs_genome_blast.Rmd      Iso-Seq transcriptome vs genome (HPC)
│   ├── _uniprot_retrieval.py              UniProt retrieval for transcriptome blast hits
│   └── _superseded/                       the 2024 version of step 01 (README inside)
├── 02_data/                         pointers only; databases/FASTAs are external
└── 03_analyses/
    ├── transcriptome-uniprot/       isoseq-transcriptome vs UniProt blastx + GO/SPID tables
    └── genome-foot/                 LOC_GO_list.txt, g.spid.txt (genome blast GO mapping)
```

## External inputs (not in repo)

Swiss-Prot (`uniprot_sprot`), the UniProt Mytilus-foot query, the genome CDS
(`GCF_036588685.1`), and the Iso-Seq transcriptome (`Mtros-hq_transcripts.fasta`). The scripts
download the UniProt files; the genome CDS and the Iso-Seq FASTA are read from copies already
on the HPC.

Step 01 needs BLAST+ 2.15.0 and the genome CDS; its search of the 53,269 CDSs would take about
a day and a half on 4 threads (estimated from a 100-CDS sample on 2026-10-03, which reproduced
the 2024 hits exactly; the 2024 run used 20 threads). Steps 02 and 03 run only in that HPC
environment. This folder has no
batch runner; the repository-level `00_run_pipeline.Rmd` reads its committed outputs.
