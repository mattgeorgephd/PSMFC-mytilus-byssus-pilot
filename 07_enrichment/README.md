# 07_enrichment

GO enrichment of the DEG lists of every contrast in `06_differential-expression` (the six
treatment-control contrasts and foot vs gill in the day-3 controls): topGO as
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
| 01 | `01_go_inputs.Rmd` | `01_go-inputs/`: one annotation table for every method; gene-set sizes; the GO release and package versions (`RUN_provenance.txt`) |
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

- **Gene sets.** For each of the 7 contrasts (`06 .../DEG_lists/contrasts.csv`) the universe is
  the genes DESeq2 gave an adjusted p (independent filtering leaves the others at NA, so they
  could never be DEGs). Up- and down-regulated DEGs (padj < 0.05) are tested separately against
  that universe.
- **Annotation.** The genome-wide BLAST (`03_blast/03_analyses/genome-foot/LOC_GO_list.txt`)
  can give several hits per LOC; the highest bitscore is kept, as in `08` and `09`. Genes are
  joined to it through `gene_key()` (`tools/gene_ids.R`). GO IDs are trimmed and checked
  against `GO.db`. The mitochondrial loci are in no universe (06 leaves them out). topGO takes the direct annotation and walks the GO graph itself;
  goseq and clusterProfiler get each gene's terms plus all their ancestors.
- **Terms tested.** At least 10 annotated genes; goseq and clusterProfiler also cap at 500.
- **Enriched.** topGO `weight01` p < 0.01 (conventionally not FDR-adjusted, because
  decorrelated p-values are not independent); goseq and clusterProfiler BH p < 0.05.
- **Ontologies.** Biological process, molecular function and cellular component, each
  tested and drawn separately.
- **Length bias.** goseq weights by median reference-transcript length. In 3' Tag-seq, one
  tag per transcript, there is little to correct: `03_goseq/goseq_pwf_TC_BP.png` shows no
  consistent trend, and goseq and clusterProfiler p-values agree almost perfectly (Spearman
  0.94 to 1.00 per run, median 0.99).

## Results in brief (TC contrasts, biological process)

On the count matrix of record (featureCounts on the Iso-Seq-extended annotation, since
2026-10-02):

| run | topGO `weight01` | goseq | clusterProfiler |
|---|---|---|---|
| Foot OA up | 5 (tRNA aminoacylation, ER unfolded protein response, regulation of translational initiation) | 26 | 28 (tRNA aminoacylation, amino-acid activation) |
| Foot DO down | 38 (cilium movement, axoneme and dynein arm assembly) | 46 | 59 (cilium movement, cilium assembly) |
| Gill OA up | 17 (glutathione metabolism, TCA cycle, protein folding, response to unfolded protein, proton-motive-force-driven ATP synthesis, cellular detoxification) | 2 | 3 (glutathione metabolism, response to ER stress, protein folding) |
| Gill DO up | 15 (ERAD pathway, glycine transport) | 2 | 3 (response to ER stress, protein N-linked glycosylation, ERAD) |
| Gill DO down | 11 (neuroblast division, adaptive immune response) | 2 | 2 (adaptive immune response, neuroblast division) |
| Foot OW down | 8 | 2 | 10 (protein glycosylation) |
| Gill OA down | 16 | 1 (extracellular structure organization) | 0 |
| the other five runs | 4-13 each | 0 | 0 |

topGO reports terms in all 12 runs. Where an FDR-controlled method also finds terms, topGO's
terms agree in part (median Jaccard 0.12 with goseq, 0.18 with clusterProfiler), as expected
from `weight01` preferring specific terms over their parents. In the runs where neither goseq
nor clusterProfiler finds anything, read topGO's list as exploratory.
`06_method-comparison/consensus_terms_TC_<ontology>.csv` lists the terms topGO and at least one
FDR-controlled method call: 34 in BP, 27 in MF and 46 in CC.

**What the count matrix of record changed** (2026-10-02, against StringTie + prepDE): the three
FDR-supported runs of before remain (Foot OA up, Foot DO down, Gill OA up), Foot DO down with
more terms (46 and 59 against 26 and 53), and Foot OW down and both Gill DO runs gain a few.
Gill OA up keeps glutathione metabolism and adds the ER stress response under FDR control, but
its nuclear ATP synthase and detoxification terms are now topGO only (FDR-supported before with
8 and 9 terms, now 2 and 3); 11 of its 19 previous topGO terms remain.

Before the mitochondrial loci were separated, Gill OA up was dominated by mitochondrial
electron transport (87 of the 88 DEGs in "ATP synthesis coupled electron transport" were
copies of mitochondrial genes); the mitochondrial genes themselves are tested in
`06_differential-expression` step 13.

The manuscript GO figure is not chosen yet; every option is drawn for every family and
ontology: the topGO, goseq and clusterProfiler dot plots, the rrvgo parent-term view and the
method comparison.

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
