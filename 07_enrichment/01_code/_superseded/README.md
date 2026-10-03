# _superseded

The DAVID / REVIGO workflow, kept as a record; replaced by `../01_go_inputs.Rmd` to
`../06_method_comparison.Rmd`.

| Script | What it did |
|---|---|
| `07-GOenrichment_listcreationforDAVID.Rmd` | wrote the UniProt background and per-contrast accession lists submitted to DAVID |
| `08-func_enrichment_DAVID_visual.Rmd` | plotted the DAVID functional-annotation charts |
| `09-GOenrichment_listcreation_REVIGO.Rmd` | wrote the GO-ID lists submitted to REVIGO |

Why they were replaced: the submitted lists were built from an older version of the DEG
tables (re-running the scripts changed several lists substantially), used one pooled
background for every contrast rather than each contrast's tested genes, mixed up- and
down-regulated genes, and depended on hand submissions to web tools.
