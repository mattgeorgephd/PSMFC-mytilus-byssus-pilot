# 09_gene-mechanics-correlation

Cross-cutting analysis linking gene expression to byssal thread mechanics, the basis for the
manuscript gene/mechanics interaction results (Section 3.4). Correlates per-sample expression of
candidate gene families (HIF, HSP, peroxidase, foot/byssus proteins) with post-stress thread
measurements.

## Layout

```
09_gene-mechanics-correlation/
├── gene-mechanics-correlation.Rproj
├── 01_code/
│   ├── 00-run_gene_mechanics_by_tissue.Rmd  driver: renders 20-23 for foot and gill
│   ├── 20-gene_mechanics_correlation.Rmd    paired table, VST, candidate set, per-gene ANCOVA (the reported test)
│   ├── 21-gene_mechanics_expanded.Rmd       modules, influence, best hits (same ANCOVA)
│   ├── 22-rna_thread_manifest_and_expression_tables.Rmd
│   ├── 23-byssus_foot_gene_list_expression.Rmd
│   ├── gene-mechanics-pipeline_DOC.md       how the four chain together; config; checks
│   └── _superseeded/11-byssal_thread_by_sample.Rmd   Grace's legacy per-sample joining (not in the chain)
├── 02_data/
│   ├── expected_animals.csv   every day-3 animal: in the foot / gill fits, or the reason it is out
│   └── HIF_GCM.csv, HSP_GCM.csv, perox_GCM.csv, foot_byss_GCM.csv   gene-family count matrices
└── 03_analyses/
    ├── gene_mechanics/     scripts 20-21, tissue-suffixed
    ├── expr_tables/        script 22
    ├── byssus_genes/       script 23
    └── knit_html/          rendered reports and run_log.csv from driver 00
```

Scripts 20 to 23 take a knit parameter `tissue` ("F" default, or "G"). Knit
`00-run_gene_mechanics_by_tissue.Rmd` to run the whole chain for both tissues, after
`02_thread-strength` scripts 1 to 4 and the `06_differential-expression` batch driver. Each script reads the previous one's CSV handoffs; none
shares an R session with another.

## Inputs (cross-folder)

All paths resolve from a `repo_root` found by walking up from `here::here()` to the first
folder that contains both `02_thread-strength/` and `06_differential-expression/`.

| input | from | read by |
|---|---|---|
| `02_thread-strength/03_analyses/02_assemble-thread-summary/thread-summary.xlsx` | thread summary, script 2 | 20, 22, 23 |
| `02_thread-strength/03_analyses/04_decompose-adhesion/mussel_response_classification.csv` | per-animal response, script 4 | 20 |
| `02_thread-strength/03_analyses/0{3,4}_*/DATA_ancova_animals.csv` | per-animal ANCOVA rows, scripts 3 and 4 | 20 (agreement check) |
| `02_thread-strength/03_analyses/01_extract-tensometer-data/thread-summary-raw-output.xlsx` | every extracted trace, script 1 | 22 |
| `06_differential-expression/02_data/gene_count_matrix_clean.csv` | counts | 20, 22, 23 |
| `06_differential-expression/03_analyses/DEG_lists/` | Tag-seq arm per sample, DEG lists, DEG annotation | 20, 22, 23 |
| `03_blast/03_analyses/genome-foot/LOC_GO_list.txt` | genome-wide BLAST / UniProt / GO annotation | 20, 23 |
| `tools/pipeline_checks.R` | `warn_unless()` checks and `RUN_provenance_<T>.txt` | 20, 21 |

## Analysis design

The reported test is one baseline-adjusted regression (ANCOVA) per gene and metric on the
per-animal values, `level_day3 ~ expression + treatment + level_baseline`, with force, area
and adhesion on the log scale (geometric means of the plaques) and extension raw. Animals
without baseline threads are not in the fits. Peak force and plaque area are the declared
primary metrics (`METRICS` in script 20, `tier` column in every output); adhesion and
extension are exploratory. BH is applied within each metric (`q_lm`) and within each tier
family (`q_family`). Candidates are the keyword matches in the genome-wide BLAST annotation
that are above the detection floor. Script 21 adds six pathway modules, the
leave-one-animal-out influence check and `best_hits_<T>.csv`, all through the same ANCOVA;
no other test is fitted. Details in `01_code/gene-mechanics-pipeline_DOC.md`.
