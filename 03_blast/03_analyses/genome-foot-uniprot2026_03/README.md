# genome-foot-uniprot2026_03

Written by `../../01_code/04_refresh_uniprot_records.Rmd`: the genome BLAST hits of 2024 (every
CDS of RefSeq GCF_036588685.1 against Swiss-Prot release 2024_04 plus the UniProt "mytilus
foot" proteins, `../genome-foot/genome_n_foot_blastx.tab`) with the UniProt records of release
2026_03. A comparison, read by nothing downstream: it holds the 2024 hits fixed and so
separates what the new records changed from what the search of 2026
(`../genome-foot-sprot2026_03/`, the table the analysis reads since 2026-10-04) changed. The
analysis read this table from 2026-10-03 to 04.

| File | Contents |
|---|---|
| `LOC_GO_list.txt` | one row per blastx hit (41,432), with the columns of `../genome-foot/LOC_GO_list.txt` in the same order: the hit (`transcript`, `sp`, `uniprot_accession.x`, the alignment columns), the record of release 2026_03 (`Reviewed` to `Interacts.with`), and from `../genome-foot/g.spid.txt` the CDS's gene (`gene`; 24,636 rows have one), its first hit's e-value (`V11`) and the 2024 subject ID (`V2.x`). `geneID` and `species` are the parts of the record's entry name, which UniProt has changed for 74 of the hit proteins since 2024 |
| `uniprot_records_2026_03.tsv` | the records of the 10,740 hit proteins as UniProt's REST service returned them on 2026-10-03 (release 2026_03), with the fields step 01 requested; the input of the offline run |
| `RUN_provenance.txt` | the release, the code commit and the input MD5s |

Against the 2024 records: every gene keeps its best hit (highest bitscore), the hit proteins'
names changed for 670 of them and their GO ID sets for 8,991. 149 of the 41,432 rows are hits on 33
proteins whose sequence UniProt has since updated (`Length` is the new length; the alignment
refers to the 2024 sequence). The 2024 table lost the first hit of the search and repeated 33
rows (`../../01_code/01_genome_blast.Rmd`, "Reproducing the 2024 tables"); this one has each hit
once.
