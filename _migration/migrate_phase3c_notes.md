# Phase 3c migration notes (final structural phase)

Completes the absorb of Grace's expression-analysis repo and dissolves
`Byssus-expression-analysis/`, `analyses/`, `code/`, and `data/`. After this, the repo's
structural reorganization is done: every analysis lives in a `<name>/` folder with
`01_code` / `02_data` / `03_analyses` and a `.Rproj`.

## What it does

- `differential-expression/` <- 15 DESeq2/DEG scripts; count matrices + treatment table +
  sample metadata into `02_data/`; `DEG_lists/` (174) into `03_analyses/`.
- `enrichment/` <- 3 DAVID/REVIGO scripts; `Func_annot_DAVID/`, `Revigo_results/`,
  `uniprot_BG_DAVID.txt` into `03_analyses/`.
- `gene-annotation/` <- 4 annotation scripts (your `06-annotation` to `_superseded/`);
  `Foot_proteins.txt` into `02_data/`; `Top_gene_summaries/` into `03_analyses/`.
- `gene-mechanics-correlation/` <- your `20-` + Grace's `11-`; the four `*_GCM.csv` into
  `02_data/`; the two gene plots into `03_analyses/`.
- Grace's three READMEs -> `_migration/byssus-provenance/`. Her `.Rhistory`,
  `project-template.Rproj`, `.ipynb_checkpoints`, and the dissolving `data/`/`analyses/`
  gitignores/READMEs -> `_to_delete/`.

## Apply

1. `git switch -c reorg/phase3c`
2. `bash migrate_phase3c.sh` from the repo root.
3. Copy the READMEs into place (mapping below).
4. `git status && git diff`, commit, push.

No script rewrites this phase, so nothing to copy into `01_code/` beyond the READMEs.

## README mapping

| File | Destination |
|------|-------------|
| `differential-expression_README.md` | `differential-expression/README.md` |
| `differential-expression_02_data_README.md` | `differential-expression/02_data/README.md` |
| `differential-expression_03_analyses_README.md` | `differential-expression/03_analyses/README.md` |
| `enrichment_README.md` | `enrichment/README.md` |
| `enrichment_02_data_README.md` | `enrichment/02_data/README.md` (also makes the empty folder exist) |
| `enrichment_03_analyses_README.md` | `enrichment/03_analyses/README.md` |
| `gene-annotation_README.md` | `gene-annotation/README.md` |
| `gene-annotation_02_data_README.md` | `gene-annotation/02_data/README.md` |
| `gene-annotation_03_analyses_README.md` | `gene-annotation/03_analyses/README.md` |
| `gene-mechanics-correlation_README.md` | `gene-mechanics-correlation/README.md` |
| `gene-mechanics-correlation_02_data_README.md` | `gene-mechanics-correlation/02_data/README.md` |
| `gene-mechanics-correlation_03_analyses_README.md` | `gene-mechanics-correlation/03_analyses/README.md` |

## Judgment calls flagged

- `13-Hisat` treated as authoritative, `07-HiSat_GL` superseded (3b).
- `Top_gene_summaries/` placed in `gene-annotation/03_analyses/` though `16-top_DEGs`
  (differential-expression) is a co-writer; it is a shared DE/annotation product.
- Scripts moved as a faithful archive carrying original paths.

## Phase 4: make the absorbed R analyses runnable

The 25 transcriptomics scripts read/write the committed CSVs but via Grace's original paths
(absolute `/home/shared/.../byssus-exp-analysis/output|data/...` or repo-root-relative
`output/...`, `data/...`). To run them in this layout each needs repointing to
`here()` / `02_data` / `03_analyses`, with cross-folder reads (DEG lists, GO mapping, thread
data) via a `repo_root` pointer. Recommended order, done one folder per push like the earlier
phases:

1. `differential-expression` (foundation: count matrix -> DESeq2 -> DEG lists).
2. `enrichment` and `gene-annotation` (read DE outputs + blast GO).
3. `gene-mechanics-correlation` (cross-cutting: GCMs + thread-strength + DE).

Separately, a top-level `README.md` mapping the canonical structure would be a good capstone;
I can draft it whenever you want.
