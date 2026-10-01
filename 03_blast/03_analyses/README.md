# 03_analyses

| Subfolder | Produced by | Contents |
|-----------|-------------|----------|
| `transcriptome-uniprot/` | `01_code/_uniprot_retrieval.py` (+ blastx) | isoseq-transcriptome vs UniProt blastx (`Mtros-hq-uniprot_blastx.tab`), GO and SPID tables, UniProt id-mapping |
| `genome-foot/` | `01_code/01_genome_blast.Rmd` | `LOC_GO_list.txt`, `g.spid.txt`: gene-to-GO / SwissProt-ID mapping from the genome blast |

`genome-foot/` is the bridge from blast hits to GO terms. `LOC_GO_list.txt` is read by
`06_differential-expression/01_code/06_join_annotation.Rmd`, `07_enrichment/01_code/01_go_inputs.Rmd`
(the GO annotation of every enrichment method) and `09_gene-mechanics-correlation` scripts 01
and 04.
