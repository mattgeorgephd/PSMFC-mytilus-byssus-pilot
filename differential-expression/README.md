# differential-expression

DESeq2 differential expression of Tag-seq counts across treatments (OA, OW, DO) in foot and
gill tissue, using the HISAT2 + StringTie genome-based count matrix. This is the core
expression analysis behind the manuscript DEG results.

Absorbed from Grace Leuchtenberger's expression-analysis repo (now canonical here). Paths
resolve through `01_code/_paths.R`; see Runnability for which scripts knit headless.

## Layout

```
differential-expression/
├── differential-expression.Rproj
├── 01_code/                 count-matrix build, per-tissue/contrast DESeq2, shrinkage/filtration,
│                            provenance check, secretion state, sensitivity fits, file joining,
│                            DEG venn/volcano, counts, top DEGs, and the batch driver (20-*)
├── 02_data/                 count matrices, treatment design table, sample metadata
└── 03_analyses/
    └── DEG_lists/           significant-DEG tables per tissue x contrast (Foot/, Gill/, GOterms_genome/)
```

## Script order

`01_5-gene_count_matrix` (assemble counts) then `02_5_DESeq_*` (per tissue x contrast) then
`03-*_Shrinkage_filtration` (apeglm shrinkage + filtering) then
`03_5-DEG_table_provenance_check` (re-derives each contrast, TC and LC, from the committed
inputs and flags or rewrites stale result tables) then `04-File_joining` (merge with GO;
writes `DEG_join_summary.csv`, the source of the per-contrast DEG counts) then `12-DEG_venn`,
`12-Volcano-plots`, `15-number_DEGS`, `16-top_DEGs`, `19-DEG_list_cleanup`.

`01_code/20-run_differential_expression.Rmd` runs every script that knits from a fresh
session as a batch, each in its own R process, with a log per step in
`03_analyses/knit_html/` (git-ignored): `01_7`, `03_5` (with `rewrite_stale`), `02_6`, `02_7`,
`04`, `16`, `12-Volcano-plots`, `12-DEG_venn`, `15`, `19`. Run it after any change to the
inputs or the DESeq scripts, then the gene-mechanics driver.

Side scripts, none of which changes the primary lists:

- `01_7-secretion_state` labels every foot sample `on` / `off` for the byssal
  plaque-protein module (mfp-2, mfp-4, foot proteins 10/12/15, tyrosinase-like 1), writes
  `02_data/secretion_state.csv`, tests whether "on" differs by arm (Fisher) and summarises
  the module at day 0 vs day 3 (`DEG_lists/secretion_state_*`). Run it after `01_5`.
- `02_6-DESeq_fourlevel_sensitivity` fits one four-level model per tissue (`~ treatment`
  over all day-3 samples, shared dispersion) and compares each stressor's DEG list with the
  pairwise fit of record (`DEG_lists/sensitivity_fourlevel/`).
- `02_7-DESeq_foot_TC_secretion_sensitivity` refits the foot stressor contrasts with the
  secretion state as a covariate (`~ secretion_state + treatment`). **Gated:
  `USE_SECRETION_STATE <- FALSE` at the top of the script; nothing runs until it is set
  TRUE.** The same flag in `gene-mechanics-correlation/01_code/20-*.Rmd` controls whether the
  state enters the gene-mechanics regressions; both default to off.

## The two controls

The design has two controls and they answer different questions. Read every contrast
against this table.

| control | animals | what a contrast against it measures | used for |
|---|---|---|---|
| **TC, treatment control** (`*_TC_*`; `treatment == "control" & day == 3`) | T126-T137, held three days under ambient conditions in the same system as the stressor arms | the stressor effect, net of time in the system, handling and secretion state | the DEG lists of record (`<X>_TC_siggene*.csv`), the gene-mechanics DEG union, enrichment |
| **LC, lab control** (`*_LC_*`; `treatment == "control" & day == 0`) | T001-T012, sampled before entering the system | the stressor **plus** three days in the system, handling, and the byssal secretion state (day-0 animals are almost all actively secreting thread; day-3 animals mostly are not, in every arm) | the time axis only: `FTC_LC` / `GTC_LC` (day-3 control vs day 0) isolate the non-stressor component; the stressor-vs-day-0 tables are context, not stressor effects |

The LC output prefixes keep their historical names; renaming them (e.g. `*_vsDay0_*`) would
break nothing in the active pipeline (no downstream script reads them) but would break the
provenance file and the paths in `03-LC_Shrinkage_filtration`, so the roles are documented
here and in the script headers instead. On 17 September 2026 the LC tables were re-derived
for the first time (`03_5`, family `LC`): six of eight reproduce exactly; `FDO_LC` and
`GDO_LC` on disk were pre-QC fits that still included T051F / T051G (the same stale-fit
signature as `GDO_TC` on 16 September: five extra genes after the count filter) and were
rewritten from the committed 22- and 23-sample inputs (DEG counts 1,047 -> 1,563 and
1,994 -> 2,458).

`03_5` and `04` are documented in `01_code/DEG-lists-provenance_DOC.md`, including the
16 September 2026 findings: the gill-hypoxia full table was a stale pre-QC fit (rewritten;
the 307-gene DEG list was already correct) and the gill-OA join had been dropping the five
mitochondrial DEGs (fixed; 711 genes now carry through to every downstream file).

The gene count matrix is produced upstream by sequence-alignment (HISAT2 + StringTie) and
placed here as the DE input, per the established handoff.

## Runnability

Paths are resolved through `01_code/_paths.R` (`here::here()` anchored on the `.Rproj`).
`03_5`, `04`, `16`, `19`, `12-Volcano-plots`, `15-number_DEGS`, `01_7`, `02_6` and `02_7`
knit cleanly from a fresh session (the batch driver runs them); `12-DEG_venn` needs `ggvenn`
and installs it from CRAN when it is missing. `01_5` reads the raw StringTie matrix, which is not in the
repository, so the clean matrix it produced is the committed input. The `02_5_DESeq_*`
scripts are the original interactive DESeq2 runs: every `DESeqDataSetFromMatrix()` in them
is preceded by an explicit `factor(..., levels = c("control", ...))` (or `day` levels `0`,
`3`; `tissue` levels `F`, `G`) and followed by a `stopifnot(resultsNames(dds)[2] == ...)`
guard, so the reference level does not depend on locale sort order, and `03_5` re-derives
and verifies every table they produced. `DEG_lists/` is both
written by the DESeq scripts and read back by the joining and summary scripts, so it is an
intermediate hub; `DEG_provenance_check.csv` records whether its result tables match their
inputs.
