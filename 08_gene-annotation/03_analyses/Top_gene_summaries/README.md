# Top_gene_summaries

Written by `../../01_code/02_uniprot_summaries.Rmd` (NCBI Entrez) and
`../../01_code/03_ortholog_lists.Rmd` (OrthoDB). Both need network access, so the runner skips
them by default (`online: false`); the tables here are from the original analysis.

| File | Contents |
|---|---|
| `<code>_topgene_summs.csv` | NCBI gene summary for each UniProt accession in the contrast's top-50 list (`06_differential-expression/03_analyses/top_DEGs/Top_50_genes/`): `UniProt_ID`, `Summary` |
| `<code>_topgene_summs_ortho.csv` | the same accessions with the bivalve species that have an ortholog in OrthoDB (`bivalve_validation`) |
| `ortho_species.tab.gz` | the OrthoDB species table used to name those species |

`<code>` is one of the six TC contrasts: FOA, FOW, FDO, GOA, GOW, GDO.

## Before relying on these tables

- **Two are from an earlier top-50 list.** `FOA_topgene_summs.csv` has 16 accessions that are
  not in the current FOA list (all 32 current ones are there), and `GOA_topgene_summs.csv` lacks
  4 of the 48 current GOA accessions and has 4 others. The `_ortho` tables carry the same
  accessions. The top-50 lists themselves are unchanged since the repository's first commit, so
  the mismatch is older than that. The other four contrasts match their lists exactly.
  Rerunning steps 02 and 03 with network access refreshes all of them.
- **Failed requests used to vanish.** Until 2026-10-01 an accession whose NCBI request failed
  was dropped from the table without a trace (the error handler assigned to a copy of the
  list). Step 02 now records it as `Error: <message>`. Whether any request failed in the
  original run cannot be told from these tables.
- `FOW_topgene_summs_ortho.csv` lists one accession twice.
