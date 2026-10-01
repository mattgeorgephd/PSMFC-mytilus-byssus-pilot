# Top_gene_summaries

Written by `../../01_code/02_uniprot_summaries.Rmd` (NCBI E-utilities) and
`../../01_code/03_ortholog_lists.Rmd` (OrthoDB). Both need network access, so the runner skips
them by default (`online: false`). The tables here were refreshed on 2026-10-01 from the
current top-50 lists (`06_differential-expression/03_analyses/top_DEGs/Top_50_genes/`).

| File | Contents |
|---|---|
| `<code>_topgene_summs.csv` | one row per UniProt accession in the contrast's top-50 list: `UniProt_ID`, `uniprot_gene`, the NCBI Gene record of that protein (`ncbi_gene_id`, `ncbi_symbol`, `ncbi_organism`), `symbol_matches` (the NCBI symbol is one of the UniProt gene names) and `Summary` (NCBI's summary, or why there is none) |
| `<code>_topgene_summs_ortho.csv` | the same rows plus `ortho_group` (the OrthoDB 12.2 group of the accession at the Metazoa level) and `bivalve_validation` (the bivalve species in that group, or why there are none) |
| `ortho_species.tab.gz` | OrthoDB 12.2's species table (`odb12v2_species.tab.gz`), used to name those species |
| `RUN_provenance_summaries.txt`, `RUN_provenance_orthologs.txt` | when and from which commit each step ran, settings, input MD5s |

`<code>` is one of the six TC contrasts: FOA, FOW, FDO, GOA, GOW, GDO.

## The 2026-10-01 refresh

- **Summaries.** 223 distinct accessions across the six lists; 195 are linked to an NCBI Gene
  record, 25 proteins have no gene link and 3 have no NCBI protein record. No request failed.
  `symbol_matches` is FALSE for 21 accessions, all naming differences for the same gene
  (Xenopus homeolog suffixes such as `edem3.L`, renamed symbols such as NARS to NARS1 or
  ZK1073.1 to ndrr-2, LOC placeholders), checked one by one.
- **Why the summaries changed.** The earlier tables searched NCBI Gene with the accession as
  free text and took the first hit, which can be another gene whose record mentions the
  accession (O95238, SPDEF, gave AR). Step 02 now follows the accession to its NCBI Protein
  record and that protein's linked gene. The FOA and GOA tables had also been fetched from an
  earlier top-50 list; all six now match their lists.
- **Orthologs.** Step 03 could not run as written (it joined OrthoDB's headerless species
  table on a column it does not have), so the earlier `_ortho` tables came from an older
  version of the script and OrthoDB release. They are being regenerated with OrthoDB 12.2.
