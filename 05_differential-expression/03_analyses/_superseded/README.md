# _superseded

Outputs of the previous scripts (`../../01_code/_superseded/`) and of retired analyses, kept as a
record; no current script reads them.

| Item | What replaced it |
|---|---|
| `DEG_lists/` | per-contrast count matrices and sample tables (`*_countmatrix.csv`, `*_treatmentinfo.csv`), earlier annotated GOA/FDO lists and `DEG_provenance_check.csv`. The sample sets are now defined by rule in `../DEG_lists/contrast_samples.csv` and the counts read from `../count_matrix/` |
| `gene_count_matrix_clean_preQC_space-delimited` | a space-delimited copy of the clean matrix read by the old REVIGO script; replaced by `../count_matrix/gene_count_matrix_clean.csv` |
| `LC_contrasts/` | the day-0 lab-control (LC) contrasts and their figures, retired because the day-0 feet were dissected differently (README inside); no replacement, the day-3 treatment-control contrasts stay the contrasts of record |
