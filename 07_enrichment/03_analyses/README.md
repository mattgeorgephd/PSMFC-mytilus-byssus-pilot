# 03_analyses

| Item | Produced by | Contents |
|------|-------------|----------|
| `Func_annot_DAVID/` | DAVID web tool, from the lists `07-GOenrichment_listcreationforDAVID` writes | DAVID functional annotation chart/cluster results; `08-func_enrichment_DAVID_visual` reads the `*_DAVID_chart.xlsx` files |
| `Revigo_results/` | REVIGO web tool, from the lists `09-GOenrichment_listcreation_REVIGO` writes | REVIGO summary tables (`*_revigo_table.tsv`) and copies of the submitted `*_ID_REVIGO.txt` lists |
| `uniprot_BG_DAVID.txt` | an earlier run of `07-GOenrichment_listcreationforDAVID` | UniProt background list for DAVID, not de-duplicated (6,365 lines, 5,547 unique IDs). The current script writes a de-duplicated 5,581-ID list to `DAVID_lists/uniprot_BG_DAVID.txt`; every ID here is in it |

The scripts write the submission lists to
`../../06_differential-expression/03_analyses/DEG_lists/DAVID_lists/` and `REVIGO_lists/`. The
`*_ID_REVIGO.txt` files in `Revigo_results/` are identical to the committed `REVIGO_lists/`.
