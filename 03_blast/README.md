# blast

Sequence-similarity searches used to annotate genes and to compare the Iso-Seq
transcriptome to the genome. Two lines of work:

- Genome CDS vs SwissProt + Mytilus-foot proteins (`Mtros-genome-blast.Rmd`), producing the
  gene-to-GO mapping (`LOC_GO_list.txt`, `g.spid.txt`) consumed by gene-annotation.
- Iso-Seq transcriptome vs UniProt (`uniprot-retrieval.py` + outputs) and isoseq vs genome
  (`Isoseq_vs_genome_blast.Rmd`).

Grace's blast scripts ran on an HPC workstation (`/home/shared/...`) against databases and
FASTAs not stored in this repo; they are kept verbatim as the method record. Only committed
outputs are in `03_analyses/`.

## Layout

```
blast/
├── blast.Rproj
├── 01_code/
│   ├── Mtros-genome-blast.Rmd       genome CDS vs SwissProt+foot (HPC); writes LOC_GO/g.spid
│   ├── Isoseq_vs_genome_blast.Rmd   isoseq transcriptome vs genome (HPC)
│   └── uniprot-retrieval.py         UniProt retrieval for transcriptome blast hits
├── 02_data/                         pointers only; databases/FASTAs are external
└── 03_analyses/
    ├── transcriptome-uniprot/       isoseq-transcriptome vs UniProt blastx + GO/SPID tables
    └── genome-foot/                 LOC_GO_list.txt, g.spid.txt (genome blast GO mapping)
```

## External inputs (not in repo)

SwissProt (`uniprot_sprot`), the UniProt Mytilus-foot query, the genome CDS
(`GCF_036588685.1`), and the Iso-Seq transcriptome (`Mtros-hq_transcripts.fasta`) are all
downloaded from UniProt / NCBI / owl by the scripts.
