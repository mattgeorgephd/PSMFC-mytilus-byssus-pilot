# DEG lists: provenance check and annotation join

Documentation for `03_5-DEG_table_provenance_check.Rmd` and the rewritten
`04-File_joining.Rmd` (16 September 2026), plus the one-line fix to `16-top_DEGs.Rmd`.
Both scripts run inside `differential-expression.Rproj`; paths come from `01_code/_paths.R`.

## Why these exist

Two problems were found in the committed DEG tables:

1. **`Gill/GDO_TC_apeglm.csv` was from a different DESeq2 fit than `Gill/GDO_TC_siggene_apeglm.csv`.**
   The full table had 11,007 rows and 199 genes at `padj < 0.05`; the DEG list had 307 genes,
   and none of its baseMean values appeared in the full table. Re-deriving the contrast from
   the committed 23-sample count matrix reproduces the **307-gene list exactly** (same genes,
   log2FC within 5e-8, padj within 1e-11). The 11,007-row table is reproduced only when
   T051G is added back, and `01_5-gene_count_matrix.Rmd` (line 83) removed T051F/T051G for
   QC. So the DEG list was current and the full table was the stale, pre-QC file. The full
   table has been regenerated from the 23-sample fit; the DEG list and the manuscript count
   (307) are unchanged.
2. **The gill-OA annotation join dropped the five mitochondrial DEGs.** `04-File_joining.Rmd`
   built its join key as `str_extract(gene, "LOC\\d+")`, which is `NA` for the mitochondrial
   genes (`gene-ND2|ND2`, ...), and the GOA chunk alone filtered `!is.na(LOC_ID)` before the
   join. ND1, ND2, ND6, CYTB and ATP6 (all `padj < 0.05` in gill OA) were therefore absent
   from `GOA_sigs_merged.csv`, from the annotated file, from the zenodo file and from the
   gill-OA top-50 table, even though the master annotation table carries all 13
   mitochondrial protein-coding genes keyed by symbol.

## `03_5-DEG_table_provenance_check.Rmd`

**Purpose.** Re-derive every treatment-vs-treatment-control contrast from the committed
inputs and say whether the two result files on disk came from that fit.

**Inputs.** `02_data/gene_count_matrix_clean.csv` (master matrix, 129 samples) and
`03_analyses/DEG_lists/<Foot|Gill>/<X>_TC_treatmentinfo.csv` for the six contrasts FOA, FOW,
FDO, GOA, GOW, GDO. The per-contrast `<X>_TC_countmatrix.csv` files are compared with the
master-matrix subset but not used as input, because `FOA` and `GOA` were written after a
`data.frame()` conversion that replaced the gene IDs with row numbers.

**Steps per contrast** (identical to `02_5_DESeq_*` and `03-TC_shrinkage_filtration`):
`DESeqDataSetFromMatrix(~ treatment)` with `treatment` a factor whose first level is
`control`; `DESeq()` on all 47,806 genes; keep genes with >= 10 counts in at least a third
of samples; `lfcShrink(coef = 2, type = "apeglm")`; DEGs are `padj < 0.05`.

**Outputs.** `03_analyses/DEG_lists/DEG_provenance_check.csv`, one row per contrast:

| column | meaning |
|---|---|
| `n_samples`, `genes_after_filter` | sample count and genes passing the count filter |
| `n_DEG_rederived` | genes at `padj < 0.05` in the re-derived fit (the number to report) |
| `countmatrix_matches_master` | the committed per-contrast count matrix holds the same counts as the master subset |
| `countmatrix_gene_ids_lost` | that file lost its gene IDs (`FOA`, `GOA`) |
| `full_*`, `siggene_*` | rows on disk vs re-derived, whether the gene sets match, and the largest baseMean / log2FC / padj differences on shared genes |
| `full_table_consistent`, `siggene_consistent` | same gene set and baseMean agreement within 1e-6 |

With `params$rewrite_stale = TRUE` any inconsistent file is overwritten with the re-derived
table in script 03's format (`write.table(row.names = FALSE)`). `gene_counts_<X>_TC.csv` is
not touched (its `normal` and `ashr` rows need estimators this script does not run).

**Result on 16 September 2026.** All six DEG lists reproduce exactly; five full tables
reproduce exactly; `GDO_TC_apeglm.csv` was stale and was rewritten (now 11,002 rows, 307 at
`padj < 0.05`, and the DEG list is a subset of it).

**Run.** Knit in RStudio (check only), or
`rmarkdown::render("01_code/03_5-DEG_table_provenance_check.Rmd", params = list(rewrite_stale = TRUE))`
to also rewrite stale files. About 30 s per contrast. Needs `DESeq2` and `apeglm`.

**A locale trap in the original scripts.** `02_5_DESeq_*` pass `treatment` as character, so
DESeq2 builds the factor with R's sort order. In an English locale `control` sorts before
`DO`/`OA`/`OW` and is the reference; in a C locale (many servers, `LANG=C`, some CI runners)
the uppercase codes sort first, coefficient 2 becomes `treatment_control_vs_<trt>`, and every
log2 fold change flips sign with no error. The check script sets the levels explicitly and
asserts the coefficient name; the same two lines should go into `02_5_DESeq_*` before they
are next run.

## `04-File_joining.Rmd` (rewritten)

**Purpose.** Join each DEG list to the BLAST/UniProt/GO master table and write the
`GOterms_genome/<X>_sigs_{merged,ID,unID}.csv` trio that the venn, volcano, count, top-DEG,
zenodo, DAVID, REVIGO and GO-slim scripts read.

**What changed.**

- One function, `join_contrast()`, applied to a six-row contrast table, instead of six
  near-identical chunks with hand-typed counts in comments.
- Join key `deg_key()`: the `LOC` identifier when present, otherwise the symbol between
  `gene-` and `|`. Mitochondrial DEGs now join in every contrast.
- `left_join(..., na_matches = "never")`: the master table has ~16,800 rows with an empty
  key (isoform hits already represented under another transcript of the same gene); a default
  join would attach all of them to any DEG whose key were missing.
- Hard stops if a DEG name cannot be keyed or a DEG list has duplicate names.
- New summary `03_analyses/DEG_lists/DEG_join_summary.csv`:

| column | meaning |
|---|---|
| `n_DEG_genes` | rows in the DEG list = genes; **this is the DEG count to report** |
| `n_merged_rows` | rows in `<X>_sigs_merged.csv`, one per gene x BLAST hit |
| `n_ID_rows`, `n_unID_rows` | merged rows with / without an annotation |
| `n_genes_annotated`, `n_genes_unannotated` | distinct genes with / without an annotation |
| `n_mito_DEGs`, `mito_DEGs` | mitochondrial genes in the DEG list |

**Counts after the fix (16 September 2026):**

| contrast | DEG genes | merged rows | annotated genes | unannotated genes | mito |
|---|---|---|---|---|---|
| FDO | 351 | 460 | 252 | 99 | 0 |
| GDO | 307 | 409 | 199 | 108 | 0 |
| FOA | 80 | 92 | 60 | 20 | 0 |
| GOA | 711 | 835 | 463 | 248 | 5 (ND2, ND1, CYTB, ATP6, ND6) |
| FOW | 153 | 201 | 104 | 49 | 0 |
| GOW | 175 | 251 | 122 | 53 | 0 |

The manuscript draft quoted the merged-row counts (409, 251, 830 gill; 460, 201, 92 foot) as
DEG numbers. The gene counts are 307, 175, 711 (gill hypoxia, warming, acidification) and
351, 153, 80 (foot).

**Downstream effects of the GOA change.** Only gill-OA files changed:
`GOA_sigs_merged.csv` and `GOA_sigs_ID.csv` (+5 rows), `clean_zenodo_files/GOA_sigs_zenodo`
(+5 rows), and `gene-annotation/.../Top_50_genes/GOA_topgenes.csv`, where ATP6 (log2FC 2.26),
ND1 (1.92), ND6 (1.79) and ND2 (1.20) enter the 25 most up-regulated genes and displace
S-crystallin SL11 (0.96), RECQ4 (0.93), PDIA4 (0.88) and Mucin-like protein (0.86).
`DAVID_lists/GOA_ID_uniprot_DAVID.txt` already contained the five mitochondrial accessions
(it was made before the filter was added), so the DAVID results stand. Not regenerated here,
because the packages are not available in the sandbox: `goslims_genome/GOA_sigs_ID*.tab`
(re-run the GOA section of `gene-annotation/01_code/06-get_GOSlims.Rmd`).

## `16-top_DEGs.Rmd`

Added `library(readr)`: the script calls `read_csv()` but loaded only `dplyr`, `ggplot2` and
`gridExtra`, so it only ran in a session where `tidyverse` was already attached.

## Known issues left in place

- `enrichment/01_code/07-*DAVID.Rmd` and `09-*REVIGO.Rmd` read four of the six
  `_sigs_ID.csv` files with `read.csv(sep = " ")`; the files are comma-separated, so those
  chunks no longer parse. The committed DAVID lists reproduce from the current files for all
  six contrasts (GOA in a different order). The committed REVIGO lists reproduce for FOA,
  FOW and GOW (FOA and FOW differ only by one blank line) but **not** for the two hypoxia
  contrasts: FDO has 3,462 GO IDs on disk against 4,459 from the current file, GDO 2,460
  against 3,963, so `Revigo_results/` for hypoxia reflects older DEG files.
- `09-*REVIGO.Rmd` builds its background from `02_data/gene_count_matrix_clean` (no
  extension), the pre-QC 131-sample matrix that still contains T051F/T051G.
- `with_GO_terms/` holds an older copy of the `_sigs_*` trio that nothing reads.
- The zenodo files are the merged (gene x hit) tables, so a gene with several BLAST hits is
  deposited several times; a `distinct(gene)` version would be the honest deposit.
