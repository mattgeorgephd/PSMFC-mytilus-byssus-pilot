# genome-foot-sprot2026_03

Written by `../../01_code/01_genome_blast.Rmd` with `run: true` (knit 2026-10-04 in the cloud
environment, R 4.6.1): every CDS of RefSeq GCF_036588685.1 (53,269) searched with blastx
(BLAST+ 2.15.0, `-evalue 1E-20 -max_target_seqs 1`, in 200 parts with `_blastx_parts.sh`) against
UniProtKB/Swiss-Prot release 2026_03 plus the 241 Mytilus foot and byssal proteins of
`../../02_data/uniprotkb_mytilus_foot_2026_03_byssal.fasta` (575,989 sequences), each hit joined
to its UniProt record of release 2026_03 and each CDS to its gene. The search took about 31
hours on 4 threads (2026-10-03 to 04).

| File | Contents |
|---|---|
| `LOC_GO_list.txt` | one row per blastx hit (HSP; 41,588 rows, 0 missing, 0 extra, 0 repeated by the script's check), with the columns of `../genome-foot/LOC_GO_list.txt` in the same order; 29,686 CDSs and 10,776 proteins (41,138 rows on Swiss-Prot entries, 450 on the foot proteins' TrEMBL entries); 24,746 rows carry their CDS's gene (`gene`), one CDS per gene and protein, as in 2024 |
| `g.spid.txt` | one row per gene and UniProt accession (17,427): the gene, the first hit's e-value, the subject ID, the CDS and the accession |
| `genome_n_foot_blastx.tab` | the blastx table the two are built from (41,588 hits) |
| `RUN_provenance.txt` | the release, the database, BLAST+ and its settings, the code commit (3644b38, clean), R and package versions, the input MD5s |

Against the 2024 hits with the same records (`../genome-foot-uniprot2026_03/`, step 04): 16,712
genes have a best hit (highest bitscore) instead of 16,675; 188 have a different best protein, 38
have a hit for the first time and 1 has none. For 53 of the 226 changed or new best hits the
protein is one of the 45 byssal additions (`../../02_data/byssal_additions_2026_03.tsv`), among
them precollagen-P for LOC134718846 and preCOL-NG ("Nongradient byssal") for LOC134718847, two
genes that had no hit (DEGs in foot vs gill, `FG_TC`, padj 0.042 and 0.0005, and in no stressor
contrast), and proximal thread matrix protein 1 for LOC134717947 (before: collagen alpha-5(VI);
too few reads to be tested). The other changes follow from Swiss-Prot's 4,675 entries new or
changed since 2024_04: a search of every CDS against those entries alone, with the database size
of 2026_03, had predicted 147 changed best proteins and 15 first hits, and 142 and 15 of them
happened. Of the genes that are DEGs in a stressor contrast, 12 change best protein (none to a
byssal addition, none from no hit); of the foot-vs-gill DEGs, 40, 10 of them to a byssal
addition.
