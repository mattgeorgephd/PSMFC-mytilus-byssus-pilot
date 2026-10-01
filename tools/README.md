# tools

Shared R helpers sourced by the analysis folders. Base R, plus rmarkdown and ggplot2 where noted.

| File | Used by | What it does |
|---|---|---|
| `run_steps.R` | every `00_run_*.Rmd` and the repository-level `00_run_pipeline.Rmd` | renders each numbered script in its own R process with an HTML report and a log, keeps `run_log.csv`, and stops on the first failure (`render_step()`, `list_steps()`, `run_steps()`, `finish_run()`) |
| `plot_style.R` | every script that draws a figure (02, 06, 07, 08, 09) | the one set of colours and the shared ggplot theme: treatments (control grey, OA green, OW orange, DO purple), thread condition (adds baseline blue and lab reference light grey), direction (red up, blue down), tissue and foot region. Validated for colour-vision deficiency |
| `mt_encoded.R` | `07_enrichment`, `08_gene-annotation` | flags LOCs whose best BLAST hit is one of the 13 mtDNA-encoded proteins (`is_mt_encoded()`); about 140 such LOCs carry one mitochondrial signal many times over |
| `pipeline_checks.R` | `02_thread-strength` (04, 05), `09_gene-mechanics-correlation` (01, 02) | `warn_unless()` checks counted into `RUN_provenance*.txt`, `provenance_lines()` (settings, git state, R and package versions, input MD5s), `psmfc_repo_root()` |

Scripts reach this folder through their folder's `here::here()` root, `file.path(here::here(), "..", "tools", ...)`,
or through a `repo_root` they resolve first.
