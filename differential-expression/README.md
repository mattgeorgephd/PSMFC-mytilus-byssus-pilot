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
`03-*_Shrinkage_filtration` (apeglm shrinkage + filtering) then `04-File_joining` (merge with GO)
then `12-DEG_venn`, `12-Volcano-plots`, `15-number_DEGS`, `16-top_DEGs`, `19-DEG_list_cleanup`.

The gene count matrix is produced upstream by sequence-alignment (HISAT2 + StringTie) and
placed here as the DE input, per the established handoff.

## Runnability (phase 4 pending)

These scripts read/write the committed CSVs here, but using Grace's original paths: some
absolute (`/home/shared/.../byssus-exp-analysis/output/...`), some repo-root-relative
(`output/...`, `data/...`). They will not run in this layout until repointed to `here()` /
`02_data` / `03_analyses`. That rewrite is phase 4. `DEG_lists/` is both written by the DESeq
scripts and read back by the joining and summary scripts, so it is an intermediate hub.
