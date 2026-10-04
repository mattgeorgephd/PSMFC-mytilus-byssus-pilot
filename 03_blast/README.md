# 03_blast

Sequence-similarity searches used to annotate genes and to compare the Iso-Seq
transcriptome to the genome. Two lines of work:

- Genome CDS vs Swiss-Prot + Mytilus-foot proteins (`01_genome_blast.Rmd`), producing the
  gene-to-GO mapping (`LOC_GO_list.txt`, `g.spid.txt`). `LOC_GO_list.txt` of the 2026 search is read by
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
  The analysis reads the search of 2026: the same script, run in the cloud environment on
  2026-10-03 and 04 against Swiss-Prot release 2026_03 plus 241 foot and byssal proteins (the
  196 of the current "(mytilus foot)" query and 45 byssal proteins it misses, found in the
  literature and listed in `02_data/byssal_additions_2026_03.tsv`), in 200 resumable parts
  (`_blastx_parts.sh`), with the UniProt records of 2026_03; its tables are in
  `03_analyses/genome-foot-sprot2026_03/` (41,588 hits; README inside). Against the 2024 hits,
  16,712 genes have a best hit instead of 16,675, and 226 have a new or different one, 53 of them
  a byssal addition (among them precollagen-P and preCOL-NG for two genes that had no hit, DEGs
  in foot vs gill, and PTMP1). Of the genes that are DEGs in a stressor contrast, 12 change best
  protein and none to a byssal addition.
  `04_refresh_uniprot_records.Rmd` keeps the comparison that separates the new search from the
  new records: the 2024 blastx table (committed) joined to the records of UniProt release
  2026_03 for its 10,740 hit proteins, in `03_analyses/genome-foot-uniprot2026_03/` (every
  gene keeps its 2024 best hit; 670 of the hit proteins' names and 8,991 of their GO ID sets
  changed since 2024). It runs offline; nothing downstream reads it.
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
│   ├── 04_refresh_uniprot_records.Rmd     the genome hits with UniProt 2026_03 records (offline)
│   ├── _blastx_parts.sh                   step 01's blastx in resumable parts
│   ├── _build_foot_fasta.R                one-off: the 2026 search's foot proteins (196 + 45 byssal)
│   ├── _uniprot_retrieval.py              UniProt retrieval for transcriptome blast hits
│   └── _superseded/                       the 2024 version of step 01 (README inside)
├── 02_data/                         pointers only; databases/FASTAs are external
└── 03_analyses/
    ├── transcriptome-uniprot/       isoseq-transcriptome vs UniProt blastx + GO/SPID tables
    ├── genome-foot/                 LOC_GO_list.txt, g.spid.txt, blastx table (genome blast of 2024, 2024 records)
    ├── genome-foot-sprot2026_03/    the same tables of the genome blast of 2026 (read by the analysis)
    └── genome-foot-uniprot2026_03/  the 2024 hits with UniProt 2026_03 records (step 04, a comparison)
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
