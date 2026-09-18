# gene-mechanics-correlation

Cross-cutting analysis linking gene expression to byssal thread mechanics, the basis for the
manuscript gene/mechanics interaction results (Section 3.4). Correlates per-sample expression of
candidate gene families (HIF, HSP, peroxidase, foot/byssus proteins) with post-stress thread
measurements.

## Layout

```
gene-mechanics-correlation/
├── gene-mechanics-correlation.Rproj
├── 01_code/
│   ├── 20-gene_mechanics_correlation.Rmd    paired table, VST, per-gene association
│   ├── 21-gene_mechanics_expanded.Rmd       weighted regression (reported), mixed models, modules, permutation, best hits
│   ├── 22-rna_thread_manifest_and_expression_tables.Rmd
│   ├── 23-byssus_foot_gene_list_expression.Rmd
│   ├── 24-run_gene_mechanics_by_tissue.Rmd  renders 20-23 for foot and gill
│   ├── gene-mechanics-pipeline_DOC.md       how the four chain together; results; config
│   └── 11-byssal_thread_by_sample.Rmd       Grace's legacy per-sample joining (not in the chain)
├── 02_data/
│   └── HIF_GCM.csv, HSP_GCM.csv, perox_GCM.csv, foot_byss_GCM.csv   gene-family count matrices
├── 03_analyses/
│   ├── gene_mechanics/     scripts 20-21, tissue-suffixed
│   ├── expr_tables/        script 22
│   ├── byssus_genes/       script 23
│   └── knit_html/          rendered reports from driver 24 (git-ignored)
```

Scripts 20 to 23 take a knit parameter `tissue` ("F" default, or "G"). Knit
`24-run_gene_mechanics_by_tissue.Rmd` to run the whole chain for both tissues, after
`thread-strength` scripts 1 to 4. Each script reads the previous one's CSV handoffs; none
shares an R session with another.

## Inputs (cross-folder)

All paths resolve from a `repo_root` found by walking up from `here::here()`.

| input | from |
|---|---|
| `thread-strength/03_analyses/thread-summary.xlsx` | curated threads, scripts 1-2 |
| `thread-strength/03_analyses/decompose-adhesion/mussel_response_classification.csv` | per-animal response, script 4 |
| `thread-strength/03_analyses/extract-tensometer-data/thread-summary-raw-output.xlsx` | every extracted trace, script 1 |
| `differential-expression/02_data/gene_count_matrix_clean.csv` | counts |
| `differential-expression/03_analyses/DEG_lists/` | Tag-seq arm per sample, DEG lists, annotation |

## Status (18 September 2026)

45 foot / 46 gill paired animals (control 11, OA 12, OW 12, DO 10-11), 44 / 45 with
baselines. Peak force and plaque area and their change from baseline are the declared
primary metrics (`PRIMARY_METRICS`, `tier` column in every output); adhesion, extension and
the change in adhesion are exploratory. Candidates come from the genome-wide BLAST
annotation (226 foot / 279 gill keyword matches, of which 201 / 251 are above the detection
floor and tested). The reported test for every metric is the plaque-count-weighted
per-animal regression, and every gene set (candidates, six modules, the DEG union) carries a
10,000-shuffle permutation p per metric, family-wide and for the primary family
(`03_analyses/gene_mechanics/best_hits_<T>.csv`, `permutation_summary_<T>.csv`).

**Nothing survives the search in either tissue.** Foot: best candidate PDE8B vs the change
in force q = 0.24, family-wide permutation p = 0.74 (primary family 0.55); best module
tRNA-synthetases vs the change in adhesion q = 0.083, family-wide 0.36; DEG union 0.70.
Gill: candidates 0.49 (HSP70 12A vs plaque area, q = 0.091, per-metric permutation 0.08),
modules 0.72, DEG union 0.24 (Arp2/3 subunit vs the change in force, q = 0.038 within the
metric, permutation 0.077). The `byssal_structural` module is unrelated to every metric
(foot best p = 0.18): the plaque-protein genes mark thread secretion, not strength. The
byssal secretion state (`differential-expression/01_code/01_7-secretion_state.Rmd`) is
joined into the manifest and can be switched on as a covariate (`USE_SECRETION_STATE`,
default FALSE). Details in `01_code/gene-mechanics-pipeline_DOC.md` §4.
