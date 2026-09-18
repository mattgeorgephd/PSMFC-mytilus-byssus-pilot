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
│   ├── 20-gene_mechanics_correlation.Rmd    paired table, VST, candidate set, per-gene ANCOVA (the reported test)
│   ├── 21-gene_mechanics_expanded.Rmd       thread-level mixed ANCOVA (sensitivity), modules, influence, best hits
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
baselines and in the fits. The reported test is one baseline-adjusted regression (ANCOVA)
per gene and metric on the per-animal values, `level_day3 ~ expression + treatment +
level_baseline`, with force, area and adhesion on the log scale (geometric means of the
plaques) and extension raw. Peak force and plaque area are the declared primary metrics
(`METRICS` in script 20, `tier` column in every output); adhesion and extension are
exploratory. BH is applied within each metric (`q_lm`) and within each tier family
(`q_family`). Candidates come from the genome-wide BLAST annotation (226 foot / 279 gill
keyword matches, of which 201 / 251 are above the detection floor and tested). Script 21
adds the thread-level mixed ANCOVA as a sensitivity check, six pathway modules, the
leave-one-animal-out influence check and `best_hits_<T>.csv`.

**Foot: nothing.** Candidate x primary family (402 tests) minimum p = 0.0064 (PDE8B vs
force, `q_family` 0.91); no q below 0.10 in any family; the `byssal_structural` module is
unrelated to every metric (best p = 0.26): the plaque-protein genes mark thread secretion,
not strength. **Gill: leads, not results.** Heat shock 70 kDa protein 12A (LOC134718614)
vs adhesion p = 3.2e-5 (`q_lm` 0.008, `q_family` 0.016, exploratory tier) and vs force
p = 2.5e-4 (`q_lm` 0.062, `q_family` 0.12, primary), two further HSPA12A paralogs at the
top of area and extension, and in the DEG union Arp2/3 complex subunit 2 vs force
p = 3.2e-5 (`q_lm` 0.032, `q_family` 0.064); all positive slopes, all robust to the most
influential animal, all in the systemic tissue rather than the foot. The byssal secretion
state (`differential-expression/01_code/01_7-secretion_state.Rmd`) is joined into the
manifest and can be switched on as a covariate (`USE_SECRETION_STATE`, default FALSE).
Details in `01_code/gene-mechanics-pipeline_DOC.md` §4.
