# 08_gene-annotation

Functional annotation of genes: GOSlim assignment, UniProt summaries, ortholog lists, and the
top-gene summary tables. Sits downstream of the BLAST GO mapping, which reaches it through the
annotated DEG tables in `06_differential-expression`.

Authoritative annotation is Grace/Sam's (`Annotation.Rmd` and the numbered scripts); your earlier
`06-annotation.Rmd` is kept under `01_code/_superseded/`.

## Layout

```
08_gene-annotation/
├── gene-annotation.Rproj
├── 01_code/
│   ├── _paths.R                             shared paths (cross-folder reads and writes)
│   ├── 06-get_GOSlims.Rmd                   GO-slim (BP) assignment per DEG, per contrast
│   ├── 17-uniprot_summaries.Rmd             NCBI gene summaries for the top-50 genes
│   ├── 18-ortholog-lists.Rmd                bivalve orthologs for the top-50 genes (OrthoDB)
│   ├── Annotation.Rmd                       HPC blast archive (does not source _paths.R)
│   └── _superseded/06-annotation.Rmd(.md)   earlier annotation attempt (yours)
├── 02_data/
│   └── Foot_proteins.txt                    byssal foot protein reference list
└── 03_analyses/
    └── Top_gene_summaries/                  top-gene summary tables (incl. Top_50_genes/)
```

## Inputs and outputs (cross-folder)

| path | read / written by |
|---|---|
| `../06_differential-expression/03_analyses/DEG_lists/GOterms_genome/*_sigs_ID.csv` | read by 06 (comma-separated, written by `04-File_joining`) |
| `../06_differential-expression/03_analyses/DEG_lists/goslims_genome/` | written by 06 |
| `03_analyses/Top_gene_summaries/Top_50_genes/` | written by `06_differential-expression/01_code/16-top_DEGs.Rmd`, read by 17 |
| `03_analyses/Top_gene_summaries/*_topgene_summs.csv` | written by 17, read by 18 |
| `03_analyses/Top_gene_summaries/*_topgene_summs_ortho.csv`, `ortho_species.tab.gz` | written by 18 |

`Top_gene_summaries/` is written by both `16-top_DEGs` (in `06_differential-expression`) and
the scripts here, so it is a shared DE/annotation product kept in this folder.

## Runnability

Every path resolves through `01_code/_paths.R`. The scripts need network access:

- `06-get_GOSlims` needs the Bioconductor packages `GSEABase` and `GO.db` (installed on first
  run if missing) and downloads `goslim_generic.obo` from the Gene Ontology Consortium.
- `17-uniprot_summaries` queries NCBI Entrez through `rentrez`.
- `18-ortholog-lists` downloads the OrthoDB species table and queries the OrthoDB API.

`Annotation.Rmd` is the HPC blast record (`/home/shared/...` paths) and does not run outside
that environment.
