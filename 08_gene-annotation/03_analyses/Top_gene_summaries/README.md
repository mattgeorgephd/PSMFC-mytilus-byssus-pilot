# Top_gene_summaries

Written by `../../01_code/02_uniprot_summaries.Rmd` (NCBI E-utilities) and
`../../01_code/03_ortholog_lists.Rmd` (OrthoDB). Both need network access, so the runner skips
them by default (`online: false`). The tables here were refreshed on 2026-10-02 from the
current top-50 lists (`05_differential-expression/03_analyses/top_DEGs/Top_50_genes/`).

| File | Contents |
|---|---|
| `<code>_topgene_summs.csv` | one row per UniProt accession in the contrast's top-50 list: `UniProt_ID`, `uniprot_gene`, the NCBI Gene record of that protein (`ncbi_gene_id`, `ncbi_symbol`, `ncbi_organism`), `symbol_matches` (the NCBI symbol is one of the UniProt gene names) and `Summary` (NCBI's summary, or why there is none) |
| `<code>_topgene_summs_ortho.csv` | the same rows plus `ortho_group` (the OrthoDB 12.2 group of the accession at the Metazoa level) and `bivalve_validation` (the bivalve species in that group, or why there are none) |
| `ortho_species.tab.gz` | OrthoDB 12.2's species table (`odb12v2_species.tab.gz`), used to name those species |
| `RUN_provenance_summaries.txt`, `RUN_provenance_orthologs.txt` | when and from which commit each step ran, settings, input MD5s |

`<code>` is one of the six TC contrasts: FOA, FOW, FDO, GOA, GOW, GDO.

## The 2026-10-04 refresh

Both steps were last run on 2026-10-04 from the top-50 lists of the count matrix of record
(featureCounts on the Iso-Seq-extended annotation) annotated with the genome BLAST of 2026
(Swiss-Prot 2026_03 plus the foot and byssal proteins, UniProt 2026_03 records); the new best
hits put 6 accessions into the lists and took 5 out.

- **Summaries.** 239 distinct accessions across the six lists; 201 are linked to an NCBI Gene
  record, 34 proteins have no gene link and 4 have no NCBI protein record. No request failed.
  `symbol_matches` is FALSE for 16 accessions, all naming differences for the same gene
  (Xenopus homeolog suffixes such as `edem3.L`, renamed symbols such as ZK1073.1 to ndrr-2,
  MBL to MBL2 or UO to Uro, species prefixes such as `Cbr-tbh-1`, LOC placeholders), checked
  one by one.
- **Why the summaries changed.** The earlier tables searched NCBI Gene with the accession as
  free text and took the first hit, which can be another gene whose record mentions the
  accession (O95238, SPDEF, gave AR). Step 02 now follows the accession to its NCBI Protein
  record and that protein's linked gene. The FOA and GOA tables had also been fetched from an
  earlier top-50 list; all six now match their lists. The lists themselves changed on
  2026-10-01, when the 167 mitochondrial pseudogene copies left the 05 fits (7 accessions
  entered, 9 left), on 2026-10-02 with the new count matrix, and on 2026-10-04 with the new
  BLAST annotation.
- **Orthologs.** Step 03 could not run as written (it joined OrthoDB's headerless species
  table on a column it does not have), so the earlier `_ortho` tables came from an older
  version of the script and OrthoDB release. With OrthoDB 12.2, 179 of the 239 accessions have
  a Metazoa-level group with at least one bivalve member, 11 have a group without one, and 49
  are in no OrthoDB 12.2 group (proteins of species OrthoDB does not hold, such as
  *Lymnaea stagnalis* CYP10, P48416); the earlier tables listed bivalves for nearly every
  accession, which this release does not support. Each species is listed once (some have two
  assemblies in OrthoDB).
- **Provenance.** Both provenance files name the commit each step ran from, with no
  uncommitted changes, and the MD5s of the top-50 lists they read.
