# differential-expression

DESeq2 differential expression of Tag-seq counts across treatments (OA, OW, DO) in foot and
gill tissue, using the HISAT2 + StringTie genome-based count matrix. This is the core
expression analysis behind the manuscript DEG results.

Absorbed from Grace Leuchtenberger's expression-analysis repo (now canonical here). The
scripts carry their original paths and are not yet rewritten for this layout (see Runnability).

## Layout

```
differential-expression/
├── differential-expression.Rproj
├── 01_code/                 15 scripts: count-matrix build, per-tissue/contrast DESeq2,
│                            shrinkage/filtration, file joining, DEG venn/volcano, counts, top DEGs
├── 02_data/                 count matrices, treatment design table, sample metadata
└── 03_analyses/
    └── DEG_lists/           significant-DEG tables per tissue x contrast (Foot/, Gill/, GOterms_genome/)
```

## Script order

`01_5-gene_count_matrix` (assemble counts) then `02_5_DESeq_*` (per tissue x contrast) then
`03-*_Shrinkage_filtration` (apeglm shrinkage + filtering) then
`03_5-DEG_table_provenance_check` (re-derives each contrast from the committed inputs and
flags or rewrites stale result tables) then `04-File_joining` (merge with GO; writes
`DEG_join_summary.csv`, the source of the per-contrast DEG counts) then `12-DEG_venn`,
`12-Volcano-plots`, `15-number_DEGS`, `16-top_DEGs`, `19-DEG_list_cleanup`.

`03_5` and `04` are documented in `01_code/DEG-lists-provenance_DOC.md`, including the
16 September 2026 findings: the gill-hypoxia full table was a stale pre-QC fit (rewritten;
the 307-gene DEG list was already correct) and the gill-OA join had been dropping the five
mitochondrial DEGs (fixed; 711 genes now carry through to every downstream file).

The gene count matrix is produced upstream by sequence-alignment (HISAT2 + StringTie) and
placed here as the DE input, per the established handoff.

## Runnability

Paths are resolved through `01_code/_paths.R` (`here::here()` anchored on the `.Rproj`).
`03_5`, `04` and `16` knit cleanly from a fresh session. The `02_5_DESeq_*` scripts still
carry interactive-session assumptions (see `DEG-lists-provenance_DOC.md`: in particular they
rely on locale sort order to make `control` the reference level). `DEG_lists/` is both
written by the DESeq scripts and read back by the joining and summary scripts, so it is an
intermediate hub; `DEG_provenance_check.csv` records whether its result tables match their
inputs.
