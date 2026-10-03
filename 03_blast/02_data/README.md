# 02_data

Inputs of the BLAST runs recorded in `../01_code/` (they ran on the lab server; nothing here
is read by the pipeline).

| Input | Source | Stored here? |
|-------|--------|--------------|
| `uniprotkb_mytilus_foot.fasta` | UniProt REST query `(mytilus foot)`; 202 proteins, the mussel part of the BLAST database `swissprot_n_mytilusfoot`; read by `01_genome_blast.Rmd` | yes: the copy that run used, from `byssus-exp-analysis/data/` (dated 2024-08-05). The query returns different entries as UniProt changes (196 on 2026-10-03), so this copy is the record |
| `foot_proteins.fasta` | 135 nucleotide CDSs of mussel foot and byssal proteins (GenBank), the query of `03_isoseq_vs_genome_blast.Rmd` | yes: the copy that run used, from `byssus-exp-analysis/data/` (dated 2024-06-13) |
| Swiss-Prot database | UniProt FTP. The genome BLAST used release 2024_04 (its hits match that release and no other between 2023_05 and 2024_05; the 2024 script had saved its download as `uniprot_sprot_r2024_01.fasta`); `01_genome_blast.Rmd` downloads it from UniProt's archive of previous releases. The Iso-Seq annotation used release 2023_02 | no (public, by release) |
| Genome CDS (`cds_from_genomic.fasta`) | NCBI `GCF_036588685.1`; a copy is on gannet, `panopea/PSMFC-mytilus-byssus-pilot/byssus-exp-analysis/data/ncbi_dataset/` | no (120 MB) |
| Iso-Seq transcriptome (`Mtros-hq_transcripts.fasta`) | owl: `halfshell/genomic-databank/` | no (852 MB) |
