# gene-mechanics-correlation

Cross-cutting analysis linking gene expression to byssal thread mechanics, the basis for the
manuscript gene/mechanics interaction results (Section 3.4). Correlates per-sample expression of
candidate gene families (HIF, HSP, peroxidase, foot/byssus proteins) with post-stress thread
measurements.

Combines your `20-gene_mechanics_correlation.Rmd` with Grace's `11-byssal_thread_by_sample.Rmd`.

## Layout

```
gene-mechanics-correlation/
├── gene-mechanics-correlation.Rproj
├── 01_code/
│   ├── 20-gene_mechanics_correlation.Rmd    your correlation analysis
│   └── 11-byssal_thread_by_sample.Rmd       Grace's per-sample thread/expression joining
├── 02_data/
│   └── HIF_GCM.csv, HSP_GCM.csv, perox_GCM.csv, foot_byss_GCM.csv   gene-family count matrices
└── 03_analyses/
    └── foot_byss_gene_plot.pdf, gill_byss_gene_plot.pdf
```

## Inputs (cross-folder)

This folder reads the gene-family count matrices in `02_data/`, the full count matrix and DEG
lists from `../differential-expression/`, and the thread measurements from `../thread-strength/`.
Being cross-cutting, its phase-4 rewrite will use a `repo_root` pointer (like summary-plots).

## Status

Per prior analysis: 31 paired animals (foot transcriptome + post-stress thread measurements)
across three treatment groups; no gene-mechanics correlation survives FDR correction, with HIF-1a
the strongest raw signal but sensitive to thread-count weighting. Treat the transcriptome-to-
attachment link as a working hypothesis until this analysis is finalized.

## Runnability (phase 4 pending)

Scripts carry original paths and are not yet repointed to this layout.
