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
About ten minutes, six of them in script 03.

| step | script | reads | writes to `03_analyses/` |
|---|---|---|---|
| 01 | `01_clean_count_matrix.Rmd` | `05_sequence-alignment/03_analyses/prepDE/gene_count_matrix.csv`, `02_data/` sample sheets | `count_matrix/` |
| 02 | `02_define_contrasts.Rmd` | the sample table | `DEG_lists/contrasts.csv`, `contrast_samples.csv` |
| 03 | `03_deseq_contrasts.Rmd` | counts, contrasts | `dds/` (fitted objects, git-ignored), `DEG_lists/filter_summary.csv`, PCA plots in `figures/` |
| 04 | `04_shrinkage_filtration.Rmd` | `dds/` | `DEG_lists/<Foot,Gill,Foot_vs_Gill>/`: apeglm tables, DEG lists, MA plots; `DEG_lists/DEG_counts.csv` |
| 05 | `05_fourlevel_sensitivity.Rmd` | counts, TC DEG lists | `DEG_lists/sensitivity_fourlevel/` |
| 06 | `06_join_annotation.Rmd` | TC DEG lists, `03_blast/03_analyses/genome-foot/LOC_GO_list.txt` | `DEG_lists/GOterms_genome/`, `DEG_lists/DEG_join_summary.csv` |
| 07 | `07_top_degs.Rmd` | annotated TC DEGs | `top_DEGs/Top_50_genes/` |
| 08-10 | `08_deg_venn.Rmd`, `09_volcano_plots.Rmd`, `10_number_degs.Rmd` | annotated TC DEGs | the manuscript TC figures in `figures/` (`TC_venn_*`, `TC_volcano_*`, `TC_DEG_numbers.png`) |
| 11 | `11_deg_figures_all_contrasts.Rmd` | every contrast's tables | `figures/DEG_counts_all_contrasts.png`, `volcano_<TC,LC,FG>.png`, `TC_vs_LC_overlap.png`; `DEG_lists/DEG_overlap_TC_LC.csv` |
| 12 | `12_deg_list_cleanup.Rmd` | annotated TC DEGs | `DEG_lists/GOterms_genome/clean_zenodo_files/` |

Packages: DESeq2, apeglm, ashr, tidyverse, gridExtra, ggvenn (installed by step 08 if
missing), here, rmarkdown.

## Contrasts

Every contrast is defined in one place, `02_define_contrasts.Rmd`, from the sample table:

| family | contrasts | samples | design | of record? |
|---|---|---|---|---|
| **TC** | `<T><X>_TC`, X = OA, OW, DO | stressor at day 3 + treatment control at day 3 | `~ treatment`, reference `control` | **yes**: the stressor effect |
| LC | `<T><X>_LC` | stressor at day 3 + lab control at day 0 | `~ treatment` | no |
| LC | `FTC_LC`, `GTC_LC` | treatment control at day 3 + lab control at day 0 | `~ day`, reference `0` | no |
| FG | `FG_TC`, `FG_LC` | foot + gill of the day-3 (TC) or day-0 (LC) controls | `~ tissue`, reference foot | no |

T is `F` (foot) or `G` (gill). Each fit keeps genes with at least 10 counts in at least a third
of its samples, shrinks log2 fold changes with apeglm, and calls a DEG at padj < 0.05.

| contrast | DEGs (up / down) | | contrast | DEGs (up / down) |
|---|---|---|---|---|
| Foot OA (TC) | 80 (58 / 22) | | Gill OA (TC) | 711 (375 / 336) |
| Foot OW (TC) | 153 (109 / 44) | | Gill OW (TC) | 175 (96 / 79) |
| Foot DO (TC) | 351 (188 / 163) | | Gill DO (TC) | 307 (170 / 137) |
| Foot, day-3 vs day-0 control | 1996 (1266 / 730) | | Gill, day-3 vs day-0 control | 3597 (1824 / 1773) |
| Gill vs foot, day-3 controls | 6169 (4009 higher in gill / 2160 higher in foot) | | Gill vs foot, day-0 controls | 6889 (4774 / 2115) |

The full table is `03_analyses/DEG_lists/DEG_counts.csv`.

**The day-3 treatment control is the control of record.** The treatment controls (T126-T137)
spent three days in the same system as the stressor arms, so a contrast against them isolates
the stressor. The lab controls (T001-T012) were sampled before entering the system, so an LC
contrast carries the stressor plus three days in the system and handling; `FTC_LC` and
`GTC_LC` measure that second part alone (1996 and 3597 DEGs, more than any TC contrast). LC
and FG results are computed and drawn so they can be inspected, but no downstream result of
record rests on them.

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

- **Mitochondrially encoded proteins.** About 140 LOCs in the genome annotation have one of the
  13 mtDNA-encoded proteins as best BLAST hit (most sit on unplaced scaffolds; see
  `tools/mt_encoded.R`). In Gill OA (TC), 116 of the 711 DEGs are such LOCs, all up by about
  1.4-fold; no other TC contrast has any. They most likely carry one mitochondrial-transcript
  signal, so they inflate that contrast's DEG count and dominate its GO results.
- **Byssal secretory genes and day.** Plaque genes such as foot protein-4 variant-1
  (LOC134711106) and byssal peroxidase-like 4 (LOC134692428) are expressed in the day-0 foot
  libraries (median about 165 counts) and absent from most day-3 ones, controls included
  (zero counts in 35 of 46). The day-0 vs day-3 difference therefore includes byssal
  secretion itself, which LC contrasts mix into the stressor effect.
- **Removed libraries.** T051F and T051G were removed at QC; T047 has no foot library.

## Layout

```
06_differential-expression/
├── differential-expression.Rproj
├── 01_code/
│   ├── 00_run_differential_expression.Rmd   batch runner
│   ├── 01_...Rmd ... 12_...Rmd              the steps above
│   ├── _paths.R                             paths every script sources
│   └── _superseded/                         the previous scripts (README inside)
├── 02_data/                                 the raw sample sheets
└── 03_analyses/
    ├── count_matrix/        clean counts, sample table, library crosswalk (step 01)
    ├── DEG_lists/           contrast definitions, DESeq2 tables and DEG lists (steps 02-06, 11-12)
    ├── dds/                 fitted DESeq2 objects (step 03; git-ignored)
    ├── figures/             PCA, DEG counts, volcano, Venn and overlap figures (steps 03, 08-11)
    ├── top_DEGs/            top-50 DEGs per TC contrast (step 07)
    ├── _superseded/         per-contrast inputs written by the old scripts
    └── knit_html/           runner reports and logs (git-ignored)
```

Each folder has its own README.
