# 01_code/_superseded

Retired BLAST scripts, kept as records. Do not run them.

| Script | What it did | Replaced by |
|---|---|---|
| `01_genome_blast_2024.Rmd` | the genome CDS blastx against Swiss-Prot plus the UniProt "mytilus foot" proteins, as run on the lab server in June to August 2024 (moved here unchanged on 2026-10-03). It wrote `../../03_analyses/genome-foot/LOC_GO_list.txt` and `g.spid.txt`, but as kept it does not run: it joins an object, `masterID`, that it never makes, reads the foot-protein table into `mytfoot_uniprot_info` and then uses `mytilus_foot_uniprot_info`, saves the UniProt downloads in its working folder and reads them from `data/`, and loads no packages | `../01_genome_blast.Rmd`, which runs and, from the 2024 blastx table, rebuilds `g.spid.txt` byte for byte, and `LOC_GO_list.txt` byte for byte with `as_2024: true` when given the 2024 UniProt record tables (rebuilt from the committed file for the test, since UniProt no longer serves them) |
