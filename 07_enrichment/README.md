# 07_enrichment

GO enrichment of the DEG lists of every contrast in `06_differential-expression`: topGO as
the method of record, goseq and clusterProfiler for comparison, rrvgo to collapse redundant
terms the way REVIGO does, and a comparison of the methods. It replaces the DAVID and REVIGO
web workflow, whose lists were built from an older version of the DEG tables and whose results
are kept in `03_analyses/_superseded/`.

Paths resolve through `01_code/_paths.R` (`here::here()` anchored on `enrichment.Rproj`). The
folder reads from `06_differential-expression`, `03_blast` and `05_sequence-alignment` and
writes only to its own `03_analyses/`.

## How to run

Open `enrichment.Rproj` and knit `01_code/00_run_enrichment.Rmd` (or let the repository-level
`00_run_pipeline.Rmd` do it, after 06). About nine minutes on four cores.

| step | script | writes to `03_analyses/` |
|---|---|---|
| 01 | `01_go_inputs.Rmd` | `01_go-inputs/`: one annotation table for every method; gene-set sizes |
| 02 | `02_topgo.Rmd` | `02_topgo/`: topGO `weight01` Fisher, **the enrichment of record** |
| 03 | `03_goseq.Rmd` | `03_goseq/`: goseq Wallenius with transcript-length bias; the bias diagnostic |
| 04 | `04_clusterprofiler.Rmd` | `04_clusterprofiler/`: `enricher` per run, merged into a compareCluster view |
| 05 | `05_rrvgo.Rmd` | `05_rrvgo/`: topGO terms reduced to clusters of similar terms |
| 06 | `06_method_comparison.Rmd` | `06_method-comparison/`: how far the methods agree; consensus terms |

`01_code/_go_helpers.R` holds what the scripts share (gene sets, annotation readers, ancestor
propagation, figure labels, the dot plot). Packages: topGO, goseq, clusterProfiler,
enrichplot, rrvgo, GOSemSim, GO.db, GSEABase, org.Hs.eg.db (only for a column rrvgo insists
on; step 05 replaces it), tidyverse, patchwork, here, rmarkdown.

## Design

- **Gene sets.** For each of the 16 contrasts (`06 .../DEG_lists/contrasts.csv`) the universe is
  the genes DESeq2 gave an adjusted p (independent filtering leaves the others at NA, so they
  could never be DEGs). Up- and down-regulated DEGs (padj < 0.05) are tested separately against
  that universe.
- **Annotation.** The genome-wide BLAST (`03_blast/03_analyses/genome-foot/LOC_GO_list.txt`)
  can give several hits per LOC; the highest bitscore is kept, as in `09`. GO IDs are trimmed
  and checked against `GO.db`. topGO takes the direct annotation and walks the GO graph itself;
  goseq and clusterProfiler get each gene's terms plus all their ancestors.
- **Terms tested.** At least 10 annotated genes; goseq and clusterProfiler also cap at 500.
- **Enriched.** topGO `weight01` p < 0.01 (conventionally not FDR-adjusted, because
  decorrelated p-values are not independent); goseq and clusterProfiler BH p < 0.05.
- **Length bias.** goseq weights by median reference-transcript length. In 3' Tag-seq, one
  tag per transcript, there is little to correct: `03_goseq/goseq_pwf_TC_BP.png` shows no
  consistent trend, and goseq and clusterProfiler p-values agree almost perfectly (Spearman
  0.996).

## Results in brief (TC contrasts, biological process)

| run | topGO `weight01` | goseq | clusterProfiler |
|---|---|---|---|
| Foot OA up | 14 (tRNA aminoacylation, amino-acid transport, glucose starvation) | 23 | 31 |
| Foot DO down | 19 (axonemal dynein assembly, cilium movement) | 34 | 59 |
| Gill OA up | 12 (mitochondrial electron transport) | 14 | 14 |
| the other nine runs | 0-20 each | 0 | 0-1 |

Where an FDR-controlled method also finds terms, topGO's terms agree in part (Jaccard 0.09 to
0.30), as expected from `weight01` preferring specific terms over their parents. In the runs
where neither goseq nor clusterProfiler finds anything, read topGO's list as exploratory.
`06_method-comparison/consensus_terms_TC_BP.csv` lists the 30 terms that topGO and at least
one other method call.

The Gill OA up terms are carried by LOCs annotated as mitochondrially encoded proteins (87 of
the 88 DEGs in the top term, "ATP synthesis coupled electron transport"; see
`tools/mt_encoded.R`): one mitochondrial signal counted many times. The dot plots mark such
terms with a triangle, and every enriched-term table has `n_mt_encoded`.

## Layout

```
07_enrichment/
├── enrichment.Rproj
├── 01_code/
│   ├── 00_run_enrichment.Rmd    batch runner
│   ├── 01_...Rmd ... 06_...Rmd  the steps above
│   ├── _paths.R, _go_helpers.R  shared paths and helpers
│   └── _superseded/             the DAVID / REVIGO list scripts
├── 02_data/                     no stored inputs (README)
└── 03_analyses/
    ├── 01_go-inputs/ ... 06_method-comparison/   one folder per step
    ├── _superseded/             DAVID and REVIGO lists and results
    └── knit_html/               runner reports and logs (git-ignored)
```

Each folder has its own README.
