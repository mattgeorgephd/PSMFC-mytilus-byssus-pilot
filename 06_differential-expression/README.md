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
│                            provenance check, four-level sensitivity fit, file joining,
│                            DEG venn/volcano, counts, top DEGs, and the batch driver (20-*)
├── 02_data/                 count matrices, treatment design table, sample metadata
└── 03_analyses/
    └── DEG_lists/           significant-DEG tables per tissue x contrast (Foot/, Gill/, GOterms_genome/)
```

## Script order

`01_5-gene_count_matrix` (assemble counts) then `02_5_DESeq_*` (per tissue x contrast) then
`03-*_Shrinkage_filtration` (apeglm shrinkage + filtering) then
`03_5-DEG_table_provenance_check` (re-derives each TC contrast of record from the committed
inputs and flags or rewrites stale result tables; LC on request) then `04-File_joining` (merge with GO;
writes `DEG_join_summary.csv`, the source of the per-contrast DEG counts) then `12-DEG_venn`,
`12-Volcano-plots`, `15-number_DEGS`, `16-top_DEGs`, `19-DEG_list_cleanup`.

`01_code/20-run_differential_expression.Rmd` runs every script that knits from a fresh
session as a batch, each in its own R process, with a log per step in
`03_analyses/knit_html/` (git-ignored): `03_5` (with `rewrite_stale`, TC contrasts), `02_6`,
`04`, `16`, `12-Volcano-plots`, `12-DEG_venn`, `15`, `19`. Run it after any change to the
inputs or the DESeq scripts, then the gene-mechanics driver.

Side script, which does not change the primary lists:
`02_6-DESeq_fourlevel_sensitivity` fits one four-level model per tissue (`~ treatment` over
all day-3 samples, shared dispersion) and compares each stressor's DEG list with the
pairwise fit of record (`DEG_lists/sensitivity_fourlevel/`).

## Control of record

The **day-3 treatment control** is the control of record for gene expression. Every
stressor effect is a contrast against it.

| control | animals | what a contrast against it measures | used for |
|---|---|---|---|
| **TC, treatment control, of record** (`*_TC_*`; `treatment == "control" & day == 3`) | T126-T137, held three days under ambient conditions in the same system as the stressor arms | the stressor effect, net of time in the system and handling | the DEG lists of record (`<X>_TC_siggene*.csv`), the gene-mechanics DEG union, enrichment |
| **LC, lab control, not of record** (`*_LC_*`; `treatment == "control" & day == 0`) | T001-T012, sampled before entering the system | the stressor **plus** three days in the system and handling | nothing downstream; `FTC_LC` / `GTC_LC` (day-3 control vs day 0) describe the non-stressor component only |

The LC scripts (`02_5_DESeq_*_LC_genome`, `03-LC_Shrinkage_filtration`) and their tables in
`DEG_lists/` are kept, but no downstream script reads them and the batch driver does not
verify them (`03_5` checks them only with `families` including `"LC"`).

`03_5` and `04` are documented in `01_code/DEG-lists-provenance_DOC.md`.

The gene count matrix is produced upstream by sequence-alignment (HISAT2 + StringTie) and
placed here as the DE input, per the established handoff.

## Runnability

Paths are resolved through `01_code/_paths.R` (`here::here()` anchored on the `.Rproj`).
`03_5`, `04`, `16`, `19`, `12-Volcano-plots`, `15-number_DEGS` and `02_6`
knit cleanly from a fresh session (the batch driver runs them); `12-DEG_venn` needs `ggvenn`
and installs it from CRAN when it is missing. `01_5` reads the raw StringTie matrix, which is not in the
repository, so the clean matrix it produced is the committed input. The `02_5_DESeq_*`
scripts are the original interactive DESeq2 runs: every `DESeqDataSetFromMatrix()` in them
is preceded by an explicit `factor(..., levels = c("control", ...))` (or `day` levels `0`,
`3`; `tissue` levels `F`, `G`) and followed by a `stopifnot(resultsNames(dds)[2] == ...)`
guard, so the reference level does not depend on locale sort order, and `03_5` re-derives
and verifies every TC table they produced. `DEG_lists/` is both
written by the DESeq scripts and read back by the joining and summary scripts, so it is an
intermediate hub; `DEG_provenance_check.csv` records whether its result tables match their
inputs.
