# 08_gene-annotation

Functional annotation of the treatment-control (TC) DEGs: their GO slim (biological process)
profile, NCBI gene summaries and bivalve orthologs of the top DEGs. It reads the annotated DEG
tables and top-50 lists of `06_differential-expression` and writes only to its own
`03_analyses/`.

Paths resolve through `01_code/_paths.R` (`here::here()` anchored on `gene-annotation.Rproj`).
The genome BLAST that produced the annotation is in `03_blast` (its record script
`02_genome_blast_uniprot_check.Rmd` used to sit here as `Annotation.Rmd`); an earlier
annotation attempt is kept under `01_code/_superseded/`.

## How to run

Open `gene-annotation.Rproj` and knit `01_code/00_run_gene_annotation.Rmd` (or let the
repository-level `00_run_pipeline.Rmd` do it, after 06).

| step | script | writes to `03_analyses/` | network |
|---|---|---|---|
| 01 | `01_go_slims.Rmd` | `goslims/`: GO slim (BP) tables per TC contrast, a summary table and heatmap | none |
| 02 | `02_uniprot_summaries.Rmd` | `Top_gene_summaries/<code>_topgene_summs.csv`: NCBI gene summaries of the top-50 DEGs | NCBI Entrez (`rentrez`) |
| 03 | `03_ortholog_lists.Rmd` | `Top_gene_summaries/<code>_topgene_summs_ortho.csv`, `ortho_species.tab.gz`: bivalve orthologs | OrthoDB |

By default the runner runs step 01 only (`online: false`); steps 02 and 03 need network access,
and their committed tables are kept. They were last run on 2026-10-01, from the current top-50
lists (`03_analyses/Top_gene_summaries/README.md`). Step 02 finds each UniProt accession's NCBI
Gene record through NCBI Protein and the protein-to-gene link (it used to take the first
free-text hit in NCBI Gene, which can be another gene), and spaces and retries its requests;
step 03 uses OrthoDB release 12.2. Both record a failed request as `Error: <message>` and write
a provenance file. An NCBI API key, if you use one, goes in the `ENTREZ_KEY` environment
variable (for example in `~/.Renviron`), never in a script; without one, NCBI's per-address
limit is shared with other users of the same address.

Packages: GSEABase, GO.db, tidyverse (step 01); rentrez (02); httr, jsonlite (03); here,
rmarkdown.

## GO slims (step 01)

Each TC DEG with a UniProt hit is mapped onto the generic GO slim: it belongs to a slim term
when any of its GO IDs is that term or a descendant. The slim is pinned in
`02_data/goslim_generic.obo` (GO release 2023-07-27, the release of the `GO.db` used here); the
GO graph comes from `GO.db`. The earlier script (`06-get_GOSlims.Rmd`) lost most of the
annotation: GO IDs kept a leading space after splitting on ";", which `GSEABase::GOCollection()`
silently drops, so each gene contributed only its first-listed GO ID, and genes were then looked
up through only the first GO ID of each slim term. Its tables are kept in
`03_analyses/_superseded/goslims_genome/`; they listed 36-65% of the gene-to-slim links implied
even by their own GO-ID column.

Each gene takes the GO IDs of its best BLAST hit (highest bitscore), as in `07_enrichment`;
until 2026-10-01 it took its first-listed hit, which gave 33 TC DEGs a different set of GO IDs
from the one 07 tested. The mitochondrial genes and their nuclear copies are no longer in the
DEG lists (they are analysed on their own in `06_differential-expression` step 13); before
they were removed they made up 89 of the 103 genes in Gill OA's "generation of precursor
metabolites and energy" cell, one mitochondrial signal counted many times.

## Layout

```
08_gene-annotation/
├── gene-annotation.Rproj
├── 01_code/
│   ├── 00_run_gene_annotation.Rmd    batch runner
│   ├── 01_go_slims.Rmd               GO slim (BP) per TC contrast, heatmap
│   ├── 02_uniprot_summaries.Rmd      NCBI gene summaries for the top-50 DEGs
│   ├── 03_ortholog_lists.Rmd         bivalve orthologs for the top-50 DEGs (OrthoDB)
│   ├── _paths.R                      shared paths
│   └── _superseded/06-annotation.Rmd(.md)   earlier annotation attempt
├── 02_data/
│   ├── Foot_proteins.txt             byssal foot-protein coding sequences (reference)
│   └── goslim_generic.obo            the pinned generic GO slim
└── 03_analyses/
    ├── goslims/                      step 01
    ├── Top_gene_summaries/           steps 02-03
    ├── _superseded/goslims_genome/   the earlier GO slim tables
    └── knit_html/                    runner reports and logs (git-ignored)
```
