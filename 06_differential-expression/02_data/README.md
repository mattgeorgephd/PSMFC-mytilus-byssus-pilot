# 02_data

Raw sample sheets; nothing here is written by a script.

| Item | Description | Read by |
|------|-------------|---------|
| `PSMFC-mytilus-byssus-pilot-RNA-tagseq_raw.csv` | Tag-seq sample sheet: library ID, treatment group (`trt`, e.g. `T_OA_d3`), RNA box and well, RNA concentration, volume and yield. Its last three, unnamed columns are free text; for the FX libraries they say "gill" and "foot_control", both wrong (see `../README.md`), and no script uses them | `01_clean_count_matrix.Rmd` |
| `psmfc_mussel_rna_summary.csv` | RNA isolation log: sample (`T01-F_PG`, `T01-F`, `T01-G`), concentration, volume, yield, tissue label, isolation date | `01_clean_count_matrix.Rmd` (crosswalk check) |

The gene count matrix comes from `../../05_sequence-alignment/03_analyses/featurecounts/`
(`05` step 04, the count matrix of record); the clean
matrix and the sample table built from these sheets are written to
`../03_analyses/count_matrix/`.
