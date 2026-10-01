# 02_extract-tensometer-data

Written by `../../01_code/02_extract_tensometer_data.Rmd`.

| Item | Contents |
|---|---|
| `thread-summary-raw-output.xlsx` | every extracted trace: sheet `data` (one row per trace: folder, thread condition, phase, day, mussel, thread, mussel-level labels from the key, peak force, integral, QC counts; no extension, see the folder README), `coverage` (traces per folder) and `pairing` (animals pulled before and after exposure). Read by script 03 and by `09_gene-mechanics-correlation` script 03 |
| `QC_plots/<source_folder>/` | one loess QC plot per trace |
