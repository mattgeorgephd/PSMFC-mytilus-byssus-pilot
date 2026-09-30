# 07_enrichment

GO-term enrichment of the DEG lists via DAVID and REVIGO: builds the gene/GO lists for DAVID,
visualizes DAVID functional annotation, and prepares REVIGO inputs.

Absorbed from Grace's repo. Paths resolve through `01_code/_paths.R` (`here::here()` anchored
on `enrichment.Rproj`).

## Layout

```
07_enrichment/
├── enrichment.Rproj
├── 01_code/
│   ├── _paths.R                                  shared paths (cross-folder reads and writes)
│   ├── 07-GOenrichment_listcreationforDAVID.Rmd  background + per-contrast UniProt lists for DAVID
│   ├── 08-func_enrichment_DAVID_visual.Rmd       plots of the DAVID functional-annotation charts
│   └── 09-GOenrichment_listcreation_REVIGO.Rmd   background + per-contrast GO-ID lists for REVIGO
├── 02_data/                 no stored inputs; reads cross-folder (below)
└── 03_analyses/
    ├── Func_annot_DAVID/    DAVID functional annotation results (read by 08)
    ├── Revigo_results/      REVIGO outputs and copies of the submitted lists
    └── uniprot_BG_DAVID.txt an earlier UniProt background list for DAVID (see 03_analyses/README.md)
```

## Inputs and outputs (cross-folder)

| path | read / written by |
|---|---|
| `../06_differential-expression/02_data/gene_count_matrix_clean.csv` | read by 07 |
| `../06_differential-expression/02_data/gene_count_matrix_clean` (space-delimited copy) | read by 09 |
| `../06_differential-expression/03_analyses/DEG_lists/GOterms_genome/*_sigs_ID.csv` | read by 07 and 09 (comma-separated, written by `04-File_joining`) |
| `../03_blast/03_analyses/genome-foot/g.spid.txt` | read by 07 |
| `../03_blast/03_analyses/genome-foot/LOC_GO_list.txt` | read by 09 |
| `../06_differential-expression/03_analyses/DEG_lists/DAVID_lists/` | written by 07 |
| `../06_differential-expression/03_analyses/DEG_lists/REVIGO_lists/` | written by 09 |

DAVID and REVIGO are web tools: the lists are submitted by hand and the committed results in
`03_analyses/` are what they returned.

## Runnability

All three scripts knit from a fresh session. 07 and 09 rewrite the DAVID / REVIGO input lists
from the current DEG tables. The committed lists (and so the DAVID / REVIGO results made from
them) were built from earlier versions of the `*_sigs_ID.csv` tables: re-running 09 changes
`FDO_ID_REVIGO.txt` and `GDO_ID_REVIGO.txt` substantially (hundreds more GO terms each) and
`GOA_ID_REVIGO.txt` slightly, and 07 reorders `GOA_ID_uniprot_DAVID.txt` (same IDs). The
other lists reproduce exactly. Re-submit to DAVID / REVIGO if you regenerate them.
