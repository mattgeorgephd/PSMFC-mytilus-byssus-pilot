# 03_analyses

| Subfolder | Produced by | Contents |
|-----------|-------------|----------|
| `transcriptome-uniprot/` | `uniprot-retrieval.py` (+ blastx) | isoseq-transcriptome vs UniProt blastx (`Mtros-hq-uniprot_blastx.tab`), GO and SPID tables, UniProt id-mapping |
| `genome-foot/` | `Mtros-genome-blast.Rmd` | `LOC_GO_list.txt`, `g.spid.txt`: gene-to-GO / SwissProt-ID mapping from the genome blast |

`genome-foot/` is the bridge from blast hits to GO terms. `LOC_GO_list.txt` is read by
`06_differential-expression/01_code/04-File_joining.Rmd`, `07_enrichment` scripts 07 and 09
(`g.spid.txt`, `LOC_GO_list.txt`) and `09_gene-mechanics-correlation` scripts 20 and 23.
