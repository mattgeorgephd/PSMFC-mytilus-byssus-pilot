# 03_analyses

| Subfolder | Produced by | Contents |
|-----------|-------------|----------|
| `transcriptome-uniprot/` | `uniprot-retrieval.py` (+ blastx) | isoseq-transcriptome vs UniProt blastx (`Mtros-hq-uniprot_blastx.tab`), GO and SPID tables, UniProt id-mapping |
| `genome-foot/` | `Mtros-genome-blast.Rmd` | `LOC_GO_list.txt`, `g.spid.txt`: gene-to-GO / SwissProt-ID mapping from the genome blast |

`genome-foot/` is the bridge from blast hits to GO terms; gene-annotation reads it.
