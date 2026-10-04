# 03_analyses

| Subfolder | Produced by | Contents |
|-----------|-------------|----------|
| `transcriptome-uniprot/` | `01_code/_uniprot_retrieval.py` (+ blastx) | isoseq-transcriptome vs UniProt blastx (`Mtros-hq-uniprot_blastx.tab`), GO and SPID tables, UniProt id-mapping |
| `genome-foot/` | `01_code/01_genome_blast.Rmd` (HPC, 2024) | `LOC_GO_list.txt`, `g.spid.txt`: gene-to-GO / SwissProt-ID mapping from the genome blast, with the UniProt records of 2024; `genome_n_foot_blastx.tab`: the blastx table of that search (41,432 hits; the gannet copy, `byssus-exp-analysis/output/`, dated 2024-08-06, MD5 521f4d6729cfc433f63875909413b640) |
| `genome-foot-sprot2026_03/` | `01_code/01_genome_blast.Rmd` (`run: true`, cloud environment, 2026-10-04) | the genome blast of 2026: Swiss-Prot 2026_03 plus 241 foot and byssal proteins, with the UniProt records of 2026_03: `LOC_GO_list.txt`, `g.spid.txt`, `genome_n_foot_blastx.tab` (41,588 hits), `RUN_provenance.txt` (README inside) |
| `genome-foot-uniprot2026_03/` | `01_code/04_refresh_uniprot_records.Rmd` | the 2024 hits with the UniProt records of release 2026_03, a comparison: `LOC_GO_list.txt`, `uniprot_records_2026_03.tsv`, `RUN_provenance.txt` (README inside) |

The genome-foot tables are the bridge from blast hits to GO terms. The analysis reads
`genome-foot-sprot2026_03/LOC_GO_list.txt`: `tools/mt_encoded.R` (the mitochondrial loci:
`04_sequence-alignment` steps 04 and 05, `05_differential-expression` step 01),
`05_differential-expression/01_code/06_join_annotation.Rmd`, `06_iso-seq-transcriptome` step 04,
`07_enrichment/01_code/01_go_inputs.Rmd` (the GO annotation of every enrichment method) and
`09_gene-mechanics-correlation` scripts 01 and 04. `genome-foot/` holds the 2024 run's tables
as the record of that search; `../01_code/01_genome_blast.Rmd` rebuilds them (see
`../README.md`), and `g.spid.txt` and the blastx table are inputs of step 04, whose table
(`genome-foot-uniprot2026_03/`) is kept as a comparison and read by nothing downstream.
