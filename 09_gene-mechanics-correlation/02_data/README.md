# 02_data

| Item | Description |
|------|-------------|
| `expected_animals.csv` | Every day-3 animal and whether it is in the foot / gill fits, or the reason it is out; script 01 checks the fitted animal set against it |
| `_superseded/HIF_GCM.csv`, `HSP_GCM.csv`, `perox_GCM.csv`, `foot_byss_GCM.csv` | Gene-family count matrices written by the legacy `01_code/_superseded/11-byssal_thread_by_sample.Rmd`; no current script reads them |

The full gene count matrix, sample table and DEG lists are read cross-folder from
`../../06_differential-expression/03_analyses/`; the genome annotation from
`../../03_blast/03_analyses/genome-foot/`; thread measurements from
`../../02_thread-strength/03_analyses/` (script 03's thread summary, script 05's per-animal
response classification, and script 02's extraction output).
