# 07_enrichment

GO enrichment of the DEG lists of every contrast in `05_differential-expression` (the six
treatment-control contrasts and foot vs gill in the day-3 controls): topGO as
the method of record, goseq and clusterProfiler for comparison, rrvgo to collapse redundant
terms the way REVIGO does, and a comparison of the methods. It replaces the DAVID and REVIGO
web workflow, whose lists were built from an older version of the DEG tables and whose results
are kept in `03_analyses/_superseded/`.

Paths resolve through `01_code/_paths.R` (`here::here()` anchored on `enrichment.Rproj`). The
folder reads from `05_differential-expression`, `03_blast` and `04_sequence-alignment` and
writes only to its own `03_analyses/`.

## How to run

Open `enrichment.Rproj` and knit `01_code/00_run_enrichment.Rmd` (or let the repository-level
`00_run_pipeline.Rmd` do it, after 05). About nine minutes on four cores.

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
on; step 05 replaces it), tidyverse, patchwork, here, rmarkdown. `GO.db` must be 3.23.1 (GO
release 2026-01-23, Bioconductor 3.23): `01_code/_paths.R` stops otherwise
(`check_go_release()`, `tools/pipeline_checks.R`; see `AGENTS.md`, How to run).

## Design

- **Gene sets.** For each of the 7 contrasts (`05 .../DEG_lists/contrasts.csv`) the universe is
  the genes DESeq2 gave an adjusted p (independent filtering leaves the others at NA, so they
  could never be DEGs). Up- and down-regulated DEGs (padj < 0.05) are tested separately against
  that universe.
- **Annotation.** The genome-wide BLAST of 2026 (Swiss-Prot 2026_03 plus the Mytilus foot and
  byssal proteins) with the UniProt records of release 2026_03, plus `03_blast` step 05's six
  genes (`03_blast/03_analyses/genome-foot-sprot2026_03-noseg/LOC_GO_list.txt`), can give several hits
  per LOC; the highest bitscore is kept, as in `08` and `09`. Genes are
  joined to it through `gene_key()` (`tools/gene_ids.R`). GO IDs are trimmed and checked
  against `GO.db`. The mitochondrial loci are in no universe (05 leaves them out). topGO takes the direct annotation and walks the GO graph itself;
  goseq and clusterProfiler get each gene's terms plus all their ancestors.
- **Terms tested.** At least 10 annotated genes; goseq and clusterProfiler also cap at 500.
- **Enriched.** topGO `weight01` p < 0.01 (conventionally not FDR-adjusted, because
  decorrelated p-values are not independent); goseq and clusterProfiler BH p < 0.05.
- **Ontologies.** Biological process, molecular function and cellular component, each
  tested and drawn separately.
- **Length bias.** goseq weights by median reference-transcript length. In 3' Tag-seq, one
  tag per transcript, there is little to correct: `03_goseq/goseq_pwf_TC_BP.png` shows no
  consistent trend, and goseq and clusterProfiler p-values agree almost perfectly (Spearman
  0.95 to 1.00 per run, median 0.99).

## Results in brief (TC contrasts, biological process)

On the count matrix of record (featureCounts on the Iso-Seq-extended annotation, since
2026-10-02), with the genome BLAST of 2026, the UniProt records of release 2026_03 and GO
release 2026-01-23 (since 2026-10-04):

| run | topGO `weight01` | goseq | clusterProfiler |
|---|---|---|---|
| Foot OA up | 4 (tRNA aminoacylation, regulation of translational initiation, protein folding) | 11 | 11 (tRNA aminoacylation, amino-acid activation) |
| Foot DO down | 30 (sperm motility, outer and inner dynein arm assembly, cilium movement, axoneme assembly) | 43 | 56 (cilium movement, microtubule-based movement) |
| Gill OA up | 20 (glutathione metabolism, carboxylic acid metabolism, protein folding, NADPH regeneration) | 6 | 6 (glutathione metabolism, sulfur compound metabolism, carboxylic acid metabolism) |
| Foot OW down | 9 (positive regulation of the ERK1 and ERK2 cascade, intracellular calcium homeostasis) | 5 | 10 (cellular homeostasis, positive regulation of the ERK1 and ERK2 cascade) |
| Gill DO up | 15 (protein folding, glycine transport, quality control of misfolded proteins) | 0 | 2 (protein transport) |
| the other seven runs | 3-12 each | 0 | 0 |

topGO reports terms in all 12 runs. Where an FDR-controlled method also finds terms, topGO's
terms agree in part (median Jaccard 0.12 with goseq, 0.08 with clusterProfiler), as expected
from `weight01` preferring specific terms over their parents. goseq and clusterProfiler
p-values agree almost perfectly (Spearman 0.95 to 1.00 per run, median 0.99). In the runs where
neither goseq nor clusterProfiler finds anything, read topGO's list as exploratory.
`06_method-comparison/consensus_terms_TC_<ontology>.csv` lists the terms topGO and at least one
FDR-controlled method call: 24 in BP, 28 in MF and 49 in CC.

**What the annotation and GO release of 2026 changed** (2026-10-04, against the 2024 BLAST
hits with their 2024 UniProt records and GO release 2023-07-27; the DEGs are the same). Of the
296 TC terms of record (all ontologies) before, 177 remain among 285; in BP the consensus terms
go from 34 to 24, 20 of them kept. Most of the change comes from UniProt's re-annotation of the
hit proteins (8,991 of their 10,740 GO ID sets changed since 2024): the new records alone keep
186 of the 296 terms, the new GO release alone 271 (17 of the 25 it loses are terms GO made
obsolete), and the new search, with the records and GO release held at 2026, 275 of 285. The
FDR-supported runs of Foot OA up, Foot DO down, Gill OA up and Foot OW down remain, with the same
themes; Foot OW down gains "positive regulation of the ERK1 and ERK2 cascade". The ER stress
signal weakens: under FDR control, Foot OA up loses the ER unfolded protein response (topGO p
0.00016 before, 0.060 now), Gill OA up "response to endoplasmic reticulum stress" (0.0016,
0.083) and protein folding (still in topGO), Gill DO up its ER stress, ERAD and N-linked
glycosylation terms (ER stress 0.0078, 0.32; N-linked glycosylation and protein folding still
in topGO), and Gill DO down its two terms (adaptive immune response, neuroblast division).

On the count matrix of record against the previous one (StringTie + prepDE, 2026-10-02): the
FDR-supported runs of before remained (Foot OA up, Foot DO down, Gill OA up), and Foot OW down
and both Gill DO runs gained a few terms; Gill OA up's nuclear ATP synthase and detoxification
terms became topGO only.

Before the mitochondrial loci were separated, Gill OA up was dominated by mitochondrial
electron transport (87 of the 88 DEGs in "ATP synthesis coupled electron transport" were
copies of mitochondrial genes); the mitochondrial genes themselves are tested in
`05_differential-expression` step 13.

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
