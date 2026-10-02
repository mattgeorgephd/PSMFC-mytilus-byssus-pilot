# 06_differential-expression

DESeq2 differential expression of the Tag-seq counts across treatments (OA, OW, DO) in foot
and gill, from the HISAT2 + StringTie genome-based count matrix of `05_sequence-alignment`.
This is the expression analysis behind the manuscript DEG results, and it feeds GO enrichment
(`07`), annotation (`08`) and the gene-mechanics associations (`09`).

Absorbed from Grace Leuchtenberger's expression-analysis repository (now canonical here).
Paths resolve through `01_code/_paths.R` (`here::here()` anchored on
`differential-expression.Rproj`). The folder reads only from `05_sequence-alignment` and
`03_blast` and writes only to its own `03_analyses/`.

## How to run

Open `differential-expression.Rproj` and knit `01_code/00_run_differential_expression.Rmd`
(or let the repository-level `00_run_pipeline.Rmd` do it). It runs the numbered scripts in
order, each in a fresh R process, with an HTML report and a log per step in
`03_analyses/knit_html/` (git-ignored) and `run_log.csv`; its `steps` parameter runs a subset.
About fifteen minutes, six of them in script 03 and four in script 13.

| step | script | reads | writes to `03_analyses/` |
|---|---|---|---|
| 01 | `01_clean_count_matrix.Rmd` | `05_sequence-alignment/03_analyses/prepDE/gene_count_matrix.csv` and `hisat/t_data.ctab`, `02_data/` sample sheets, the BLAST table | `count_matrix/`, including `mitochondrial_loci.csv` |
| 02 | `02_define_contrasts.Rmd` | the sample table | `DEG_lists/contrasts.csv`, `contrast_samples.csv` |
| 03 | `03_deseq_contrasts.Rmd` | counts (without the mitochondrial loci), contrasts | `dds/` (fitted objects, git-ignored), `DEG_lists/filter_summary.csv`, PCA plots in `figures/` |
| 04 | `04_shrinkage_filtration.Rmd` | `dds/` | `DEG_lists/<Foot,Gill,Foot_vs_Gill>/`: apeglm tables, DEG lists, MA plots; `DEG_lists/DEG_counts.csv` |
| 05 | `05_fourlevel_sensitivity.Rmd` | counts, TC DEG lists | `DEG_lists/sensitivity_fourlevel/` |
| 06 | `06_join_annotation.Rmd` | TC DEG lists, `03_blast/03_analyses/genome-foot/LOC_GO_list.txt` | `DEG_lists/GOterms_genome/`, `DEG_lists/DEG_join_summary.csv` |
| 07 | `07_top_degs.Rmd` | annotated TC DEGs | `top_DEGs/Top_50_genes/`: the top-50 tables and bar plots labelled by gene symbol |
| 08-10 | `08_deg_venn.Rmd`, `09_volcano_plots.Rmd`, `10_number_degs.Rmd` | annotated TC DEGs | the manuscript TC figures in `figures/` (`TC_venn_*`, `TC_volcano_*`, `TC_DEG_numbers.png`) |
| 11 | `11_deg_figures_all_contrasts.Rmd` | every contrast's tables | `figures/DEG_counts_all_contrasts.png`, `volcano_TC.png`, `volcano_FG.png` |
| 12 | `12_deg_list_cleanup.Rmd` | annotated TC DEGs | `DEG_lists/GOterms_genome/clean_zenodo_files/`: one row per DEG with its best-hit protein |
| 13 | `13_mitochondrial_expression.Rmd` | counts, `mitochondrial_loci.csv`, contrasts | `mitochondrial/` and `figures/MT_mitochondrial_expression.png` (manuscript figure), `figures/MT_nuclear_copy_share.png` |

Packages: DESeq2, apeglm, ashr, tidyverse, gridExtra, ggvenn (installed by step 08 if
missing), here, rmarkdown.

## Contrasts

Every contrast is defined in one place, `02_define_contrasts.Rmd`, from the sample table:

| family | contrasts | samples | design | of record? |
|---|---|---|---|---|
| **TC** | `<T><X>_TC`, X = OA, OW, DO | stressor at day 3 + treatment control at day 3 | `~ treatment`, reference `control` | **yes**: the stressor effect |
| FG | `FG_TC` | foot + gill of the day-3 treatment controls | `~ tissue`, reference foot | no |

T is `F` (foot) or `G` (gill). Each fit keeps genes with at least 10 counts in at least a third
of its samples, shrinks log2 fold changes with apeglm, and calls a DEG at padj < 0.05. The
mitochondrial loci are not in these fits (see below).

| contrast | DEGs (up / down) | | contrast | DEGs (up / down) |
|---|---|---|---|---|
| Foot OA (TC) | 75 (58 / 17) | | Gill OA (TC) | 423 (173 / 250) |
| Foot OW (TC) | 165 (122 / 43) | | Gill OW (TC) | 180 (99 / 81) |
| Foot DO (TC) | 363 (190 / 173) | | Gill DO (TC) | 310 (173 / 137) |
| Gill vs foot, day-3 controls | 5860 (3834 higher in gill / 2026 higher in foot) | | | |

The full table is `03_analyses/DEG_lists/DEG_counts.csv`. Leaving the mitochondrial loci out
changed the counts from the earlier run. Gill OA had 711 DEGs, 116 of them mitochondrial
loci; re-adjusting the earlier p-values without those rows gives 549, because removing 116
very small p-values moves every other gene down the Benjamini-Hochberg ranking, and the
refit gives 543. The other contrasts had no mitochondrial DEGs and moved by 1 to 11 genes
(Foot OA 80 to 87, Foot OW 153 to 164, Foot DO 351 to 361, Gill OW 175 to 173, Gill DO 307 to
306), because DESeq2 re-estimates the dispersion trend and its independent-filtering
threshold on the smaller gene set. On 2026-10-01 the 167 mitochondrial pseudogene copies,
which the BLAST route had missed (see below), were left out as well, and the counts changed
again: Gill OA 543 to 423 (81 of the 543 were pseudogene copies), Foot OA 87 to 75, Foot OW 164
to 165, Foot DO 361 to 363, Gill OW 173 to 180, Gill DO 306 to 310. The copies carry a large,
OA-raised share of the reads, so leaving them out also moves the size factors and the
dispersion trend of every fit, which is why Foot OA, with no pseudogene DEG, changes too.

**The day-3 treatment control is the only control.** The treatment controls (T126-T137) spent
three days in the same system as the stressor arms, so a contrast against them isolates the
stressor. The day-0 lab controls (T001-T012) are not used as a control anywhere: their feet
were dissected differently (two pieces at day 0, one at day 3), so a day-0 vs day-3
difference would mix the sampled region into the effect, and the experimenter confirms they
are not true lab controls. The former LC family (each stressor, and the day-3 control,
against the day-0 control) and the day-0 foot-vs-gill contrast are retired; their tables and
figures are in `03_analyses/_superseded/LC_contrasts/`.

## Samples and foot regions

Tissue was foot or gill. Two parts of the foot were sequenced, and the sample sheets name them
inconsistently: every animal has a library of the **phenol gland to the tip of the foot**
(IDs ending `F`; the RNA isolation log calls them `T01-F_PG`, "phenol gland"), and the twelve
day-0 animals also have a library of the **rest of the foot, without the phenol gland** (IDs
ending `FX`; the isolation log's `T01-F`, "foot"; the Tag-seq sheet's tissue column wrongly
says "gill"). `01_clean_count_matrix.Rmd` records both as foot with a `region` column and writes
`library_crosswalk.csv`, which matches every library to its isolation record. All contrasts
use the phenol-gland-to-tip libraries; the FX libraries enter none, because the two regions
differ strongly (in the same 12 animals 3,174 of 7,393 genes differ, among them byssal
tyrosinases and collagens over a thousand-fold), and there is no day-3 FX library to compare.

## Things to know before interpreting

- **Mitochondrial genes are analysed on their own.** The reference holds the mitochondrial
  genome (NC_007687.1, 12 protein genes) and 293 loci on unplaced scaffolds that are copies of
  those genes: 126 protein-coding LOCs and 167 pseudogenes (the pseudogenes, which have no CDS
  for the BLAST route to see, are found by NCBI's names for them;
  `05_sequence-alignment/02_data/annotation_mt_like_loci.csv`). Reads of one mitochondrial transcript are split between the gene and
  its copies, so in the genome contrasts one signal was counted many times (116 of Gill OA's
  711 DEGs, all up about 1.4-fold, and most of its top GO terms). Script 01 lists these loci
  (`count_matrix/mitochondrial_loci.csv`, 310 rows), scripts 03 and 05 leave them out, and
  script 13 tests each protein once, as the sum of its gene and copies, against the day-3
  control. Result: Gill OA raises 10 of the 12 proteins (1.2- to 1.8-fold) and their sum
  (1.4-fold), Foot OA 3 (ND1, ND2, ND5; 1.25- to 1.4-fold) and their sum (1.2-fold); OW and DO
  change none (`mitochondrial/mt_de_TC.csv`).
- **Byssal secretory genes and the day-0 dissection.** Plaque genes such as foot protein-4
  variant-1 (LOC134711106) and byssal peroxidase-like 4 (LOC134692428) are expressed in the
  day-0 foot libraries (median about 165 counts) and absent from most day-3 ones, controls
  included (zero counts in 35 of 46). Because the day-0 feet were dissected differently, this
  may reflect the region sampled as much as byssal secretion; it is one reason no contrast
  uses the day-0 libraries.
- **Outlier-replaced genes.** With 7 or more samples per group, `DESeq()` replaces a count
  with an extreme Cook's distance and refits the gene; the Wald p (and so the DEG call) comes
  from the refit, while `lfcShrink(type = "apeglm")` uses the original counts. 13 TC DEGs had
  a count replaced, and for one (Foot OW) the reported apeglm fold change is more than 1.25
  times the refit estimate; in `FG_TC` 23 and 7. This is the standard DESeq2 workflow and is
  kept; script 13 draws the refit estimate for the mitochondrial proteins, where it matters
  (Gill OA ND1 and ATP6, one outlying library, T025G).
- **Removed libraries.** T051F and T051G were removed at QC; T047 has no foot library.

## Layout

```
06_differential-expression/
├── differential-expression.Rproj
├── 01_code/
│   ├── 00_run_differential_expression.Rmd   batch runner
│   ├── 01_...Rmd ... 13_...Rmd              the steps above
│   ├── _paths.R                             paths every script sources
│   └── _superseded/                         the previous scripts (README inside)
├── 02_data/                                 the raw sample sheets
└── 03_analyses/
    ├── count_matrix/        clean counts, sample table, library crosswalk, mitochondrial loci (step 01)
    ├── DEG_lists/           contrast definitions, DESeq2 tables and DEG lists (steps 02-06, 12)
    ├── mitochondrial/       the mitochondrial proteins on their own (step 13)
    ├── dds/                 fitted DESeq2 objects (step 03; git-ignored)
    ├── figures/             PCA, DEG counts, volcano, Venn and mitochondrial figures (steps 03, 08-11, 13)
    ├── top_DEGs/            top-50 DEGs per TC contrast (step 07)
    ├── _superseded/         per-contrast inputs written by the old scripts; the retired LC contrasts
    └── knit_html/           runner reports and logs (git-ignored)
```

Each folder has its own README.
