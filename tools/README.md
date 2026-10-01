# tools

Shared R helpers sourced by the analysis folders. Base R, plus rmarkdown and ggplot2 where noted.

| File | Used by | What it does |
|---|---|---|
| `run_steps.R` | every `00_run_*.Rmd` and the repository-level `00_run_pipeline.Rmd` | renders each numbered script in its own R process with an HTML report and a log, keeps `run_log.csv`, and stops on the first failure (`render_step()`, `list_steps()`, `run_steps()`, `finish_run()`) |
| `plot_style.R` | every script that draws a figure (02, 06, 07, 08, 09) | the one set of colours and the shared ggplot theme: treatments (control grey, OA green, OW orange, DO purple), thread condition (adds baseline blue and lab reference light grey), direction (red up, blue down), tissue and foot region. Validated for colour-vision deficiency. Also `short_name()`, which shortens a UniProt protein name for a figure label: it drops the alternative names in parentheses and the `[Cleaved into: ...]` / `[Includes: ...]` sections, keeps brackets that belong to the recommended name (`Amine oxidase [flavin-containing] A`), and cuts at a word boundary |
| `gene_ids.R` | `06_differential-expression` (01, 06), `07_enrichment` (01), `09_gene-mechanics-correlation` (01, 03, 04, 05), `mt_encoded.R` | `gene_key()`: the join key of a count-matrix gene name, the LOC identifier when the name holds one (`gene-LOC1\|LOC1`, `gene-LOC1`, `STRG.10\|LOC1`), otherwise the gene ID without `gene-` (`ND2`, `Trnaa-agc-10`, `STRG.12`). Unique over the 47,806 matrix rows; it replaces four slightly different local versions, one of which collapsed the tRNA genes onto their shared name and one of which left `gene-` on the 476 names without `\|` |
| `mt_encoded.R` | `06_differential-expression` (01) | the mitochondrial loci of the count matrix: the 12 protein genes and 5 RNA genes of the mitochondrial genome (NC_007687.1) and the LOCs on unplaced scaffolds whose best BLAST hit is a mitochondrially encoded protein (`mitochondrial_loci()`, written to `06/03_analyses/count_matrix/mitochondrial_loci.csv`); `is_mt_encoded()` and `mt_protein_symbol()` |
| `pipeline_checks.R` | `02_thread-strength` (04, 05), `07_enrichment` (01), `09_gene-mechanics-correlation` (01, 02, 05) | `warn_unless()` checks counted into `RUN_provenance*.txt`, `provenance_lines()` (settings, git state, R and package versions, input MD5s), `psmfc_repo_root()` |

Scripts reach this folder through their folder's `here::here()` root, `file.path(here::here(), "..", "tools", ...)`,
or through a `repo_root` they resolve first.
