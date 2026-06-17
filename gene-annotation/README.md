# gene-annotation

Functional annotation of genes: GOSlim assignment, UniProt summaries, ortholog lists, and the
top-gene summary tables. Sits downstream of the blast GO mapping.

Authoritative annotation is Grace/Sam's (`Annotation.Rmd` and the numbered scripts); your earlier
`06-annotation.Rmd` is kept under `01_code/_superseded/`.

## Layout

```
gene-annotation/
├── gene-annotation.Rproj
├── 01_code/
│   ├── Annotation.Rmd, 06-get_GOSlims.Rmd, 17-uniprot_summaries.Rmd, 18-ortholog-lists.Rmd
│   └── _superseded/06-annotation.Rmd(.md)   earlier annotation attempt (yours)
├── 02_data/
│   └── Foot_proteins.txt                    byssal foot protein reference list
└── 03_analyses/
    └── Top_gene_summaries/                  top-gene summary tables (incl. Top_50_genes/)
```

## Inputs and outputs (cross-folder)

Annotation reads the gene-to-GO mapping in `../blast/03_analyses/genome-foot/`
(`LOC_GO_list.txt`, `g.spid.txt`) and DEG lists from `../differential-expression/`.
`Top_gene_summaries/` is written by both `16-top_DEGs` (in differential-expression) and the
UniProt-summary script here, so it is a shared DE/annotation product kept in this folder.

## Runnability (phase 4 pending)

Scripts carry original paths and are not yet repointed to this layout.
