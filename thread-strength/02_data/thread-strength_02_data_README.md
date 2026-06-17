# 02_data

Inputs consumed by the scripts in `../01_code/`. Do not write outputs here;
outputs belong in `../03_analyses/`.

| Item | Description | Read by |
|------|-------------|---------|
| `tensometer_output/control/` | Raw force/displacement `.txt` traces, baseline threads | `1_extract_tensometer_data.Rmd` |
| `tensometer_output/treatment/` | Raw force/displacement `.txt` traces, post-stress threads | `1_extract_tensometer_data.Rmd` |
| `pictures/control/`, `pictures/treatment/` | Microscope images used to measure plaque area (`pad_area`) | (manual measurement) |
| `summary.xlsx` | Curated table (force + AUC + plaque area + treatment); hand-assembled, see folder README | `2_analyze_thread_strength.Rmd` |
| `instrument-reference/` | Force-gauge LabVIEW VI, OMEGA manual, wiring map, original folder README | reference only |

Master data for `summary.xlsx`: project Google Sheet (linked in the folder README).
