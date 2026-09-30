# 02_data

| Item | Description |
|------|-------------|
| `expected_animals.csv` | Every day-3 animal and whether it is in the foot / gill fits, or the reason it is out; script 20 checks the fitted animal set against it |
| `HIF_GCM.csv` | Hypoxia-inducible factor gene-family count matrix |
| `HSP_GCM.csv` | Heat-shock protein gene-family count matrix |
| `perox_GCM.csv` | Peroxidase gene-family count matrix |
| `foot_byss_GCM.csv` | Foot/byssus protein gene-family count matrix |

The full gene count matrix and DEG lists are read cross-folder from
`../../06_differential-expression/`; the genome annotation from
`../../03_blast/03_analyses/genome-foot/`; thread measurements from
`../../02_thread-strength/03_analyses/` (script 2's thread summary, script 4's per-animal
response classification, and script 1's extraction output).
