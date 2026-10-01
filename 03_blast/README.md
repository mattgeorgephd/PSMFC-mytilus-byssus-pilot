# 03_blast

Sequence-similarity searches used to annotate genes and to compare the Iso-Seq
transcriptome to the genome. Two lines of work:

- Genome CDS vs SwissProt + Mytilus-foot proteins (`01_genome_blast.Rmd`), producing the
  gene-to-GO mapping (`LOC_GO_list.txt`, `g.spid.txt`) consumed by `06_differential-expression`
  (`06_join_annotation.Rmd`), `07_enrichment` and `09_gene-mechanics-correlation`.
  `02_genome_blast_uniprot_check.Rmd`, formerly `08_gene-annotation/01_code/Annotation.Rmd`,
  inspects a SwissProt-only blastx of the same CDS and the UniProt annotation download; it
  writes nothing.
- Iso-Seq transcriptome vs UniProt (`_uniprot_retrieval.py` + outputs) and Iso-Seq vs genome
  (`03_isoseq_vs_genome_blast.Rmd`).

Grace's blast scripts ran on an HPC workstation (`/home/shared/...`) against databases and
FASTAs not stored in this repo; they are kept verbatim as the method record. Only committed
outputs are in `03_analyses/`.

## Layout

```
03_blast/
├── blast.Rproj
├── 01_code/
│   ├── 01_genome_blast.Rmd                genome CDS vs SwissProt+foot (HPC); writes LOC_GO/g.spid
│   ├── 02_genome_blast_uniprot_check.Rmd  checks of a SwissProt-only genome blastx (HPC); writes nothing
│   ├── 03_isoseq_vs_genome_blast.Rmd      Iso-Seq transcriptome vs genome (HPC)
│   └── _uniprot_retrieval.py              UniProt retrieval for transcriptome blast hits
├── 02_data/                         pointers only; databases/FASTAs are external
└── 03_analyses/
    ├── transcriptome-uniprot/       isoseq-transcriptome vs UniProt blastx + GO/SPID tables
    └── genome-foot/                 LOC_GO_list.txt, g.spid.txt (genome blast GO mapping)
```

## External inputs (not in repo)

SwissProt (`uniprot_sprot`), the UniProt Mytilus-foot query, the genome CDS
(`GCF_036588685.1`), and the Iso-Seq transcriptome (`Mtros-hq_transcripts.fasta`) are all
downloaded from UniProt / NCBI / owl by the scripts.

None of these scripts runs outside that HPC environment, so this folder has no batch runner;
the repository-level `00_run_pipeline.Rmd` reads its committed outputs.
