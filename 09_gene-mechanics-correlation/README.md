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
│   ├── 00_run_gene_mechanics_by_tissue.Rmd  batch runner: scripts 01-04 for foot and gill
│   ├── 01_gene_mechanics_correlation.Rmd    paired table, VST, candidate set, per-gene ANCOVA (the reported test)
│   ├── 02_gene_mechanics_expanded.Rmd       modules, influence, best hits (same ANCOVA)
│   ├── 03_rna_thread_manifest_and_expression_tables.Rmd
│   ├── 04_byssus_foot_gene_list_expression.Rmd
│   ├── gene-mechanics-pipeline_DOC.md       how the four chain together; config; checks
│   └── _superseded/11-byssal_thread_by_sample.Rmd   Grace's legacy per-sample joining (not in the chain)
├── 02_data/
│   ├── expected_animals.csv   every day-3 animal: in the foot / gill fits, or the reason it is out
│   └── _superseded/           gene-family count matrices written by the legacy script; not read
└── 03_analyses/
    ├── gene_mechanics/     scripts 01 and 02, tissue-suffixed
    ├── expr_tables/        script 03
    ├── byssus_genes/       script 04
    ├── knit_html/          reports, logs and run_log.csv from runner 00 (git-ignored)
    └── _superseded/        the legacy script's two PDFs
```

Scripts 01 to 04 take a knit parameter `tissue` (`"F"` default, or `"G"`). Knit
`00_run_gene_mechanics_by_tissue.Rmd` to run the whole chain for both tissues, after the
`02_thread-strength` and `06_differential-expression` runners (the repository-level
`00_run_pipeline.Rmd` runs all three in order). Each script reads the previous one's CSV
handoffs; none shares an R session with another.

Foot means the phenol gland to the tip of the foot, the region sampled in every animal. The
twelve day-0 animals also have a library of the rest of the foot (IDs ending `FX`); no script
here uses them (see `06_differential-expression/README.md`).

## Inputs (cross-folder)

All paths resolve from a `repo_root` found by walking up from `here::here()` to the first
folder that contains both `02_thread-strength/` and `06_differential-expression/`.

| input | from | read by |
|---|---|---|
| `02_thread-strength/03_analyses/03_assemble-thread-summary/thread-summary.xlsx` | thread summary, 02 script 03 | 01, 03, 04 |
| `02_thread-strength/03_analyses/05_decompose-adhesion/mussel_response_classification.csv` | per-animal response, 02 script 05 | 01 |
| `02_thread-strength/03_analyses/0{4,5}_*/DATA_ancova_animals.csv` | per-animal ANCOVA rows, 02 scripts 04 and 05 | 01 (agreement check) |
| `02_thread-strength/03_analyses/02_extract-tensometer-data/thread-summary-raw-output.xlsx` | every extracted trace, 02 script 02 | 03 |
| `06_differential-expression/03_analyses/count_matrix/gene_count_matrix_clean.csv` | counts, 06 script 01 | 01, 03, 04 |
| `06_differential-expression/03_analyses/count_matrix/treatmentinfo_clean.csv` | treatment, day and region per library, 06 script 01 | 01, 03 |
| `06_differential-expression/03_analyses/DEG_lists/<Foot or Gill>/<T><X>_TC_siggene.csv` | TC DEG lists, 06 script 04 | 01, 03, 04 |
| `06_differential-expression/03_analyses/DEG_lists/GOterms_genome/*_sigs_ID.csv` | annotated TC DEGs, 06 script 06 | 01 |
| `03_blast/03_analyses/genome-foot/LOC_GO_list.txt` | genome-wide BLAST / UniProt / GO annotation | 01, 04 |
| `tools/pipeline_checks.R` | `warn_unless()` checks and `RUN_provenance_<T>.txt` | 01, 02 |

## Analysis design

The reported test is one baseline-adjusted regression (ANCOVA) per gene and metric on the
per-animal values, `level_day3 ~ expression + treatment + level_baseline`, with force, area
and adhesion on the log scale (geometric means of the plaques) and extension raw. Animals
without baseline threads are not in the fits. Peak force and plaque area are the declared
primary metrics (`METRICS` in script 01, `tier` column in every output); adhesion and
extension are exploratory. BH is applied within each metric (`q_lm`) and within each tier
family (`q_family`). Candidates are the keyword matches in the genome-wide BLAST annotation
that are above the detection floor. Script 02 adds six pathway modules, the
leave-one-animal-out influence check and `best_hits_<T>.csv`, all through the same ANCOVA;
no other test is fitted. Details in `01_code/gene-mechanics-pipeline_DOC.md`.
