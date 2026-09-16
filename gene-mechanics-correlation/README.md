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
│   ├── 21-gene_mechanics_expanded.Rmd       mixed models, modules, permutation, diagnostics
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
└── gene-mechanics-results-report.docx   results write-up, foot and gill
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

## Status (16 September 2026)

44 paired animals (control 10, OA 12, OW 12, DO 10), 42 with baselines. Peak force and
plaque area are the primary metrics; adhesion (their ratio) hides the stressor effect on this
dataset. Change-from-baseline metrics (`dlog_*`) are tested alongside the day-3 levels.

Foot: nothing on the level metrics survives correction. On the change in adhesion, Collagen
alpha-1(V) reaches q = 0.034 with a permutation-corrected p of 0.054 and the same sign and
rank without the control arm; the collagen module is q = 0.052. Gill (n = 45): nothing
survives, and the collagen signal has no counterpart there. Treat as
hypothesis-strengthening, not established. Details in `gene-mechanics-results-report.docx`
and `01_code/gene-mechanics-pipeline_DOC.md`.
