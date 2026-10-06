# 09_gene-mechanics-correlation

Cross-cutting analysis linking gene expression to byssal thread mechanics, the basis for the
manuscript gene/mechanics interaction results (Section 3.4). Within each arm, and with the
animal's own baseline held constant, it asks whether animals that express a gene, a DEG set,
the genes of an enriched GO term, or their mitochondrial genes more strongly make stronger or
weaker plaques. Candidate gene families: HIF, HSP, peroxidase, foot/byssus proteins.

## Layout

```
09_gene-mechanics-correlation/
├── gene-mechanics-correlation.Rproj
├── 01_code/
│   ├── 00_run_gene_mechanics_by_tissue.Rmd  batch runner: scripts 01-06 for foot and gill
│   ├── 01_gene_mechanics_correlation.Rmd    paired table, VST, candidate set, per-gene ANCOVA (the reported test)
│   ├── 02_gene_mechanics_expanded.Rmd       modules, influence, best hits (same ANCOVA)
│   ├── 03_rna_thread_manifest_and_expression_tables.Rmd
│   ├── 04_byssus_foot_gene_list_expression.Rmd
│   ├── 05_go_term_mechanics.Rmd             mechanics against DEG sets, enriched GO terms, mitochondrial share
│   ├── 06_expression_suites.Rmd             suites: programs, co-expression modules, expression components between and within arms; multi-gene prediction
│   ├── gene-mechanics-pipeline_DOC.md       how the six chain together; config; checks
│   └── _superseded/11-byssal_thread_by_sample.Rmd   Grace's legacy per-sample joining (not in the chain)
├── 02_data/
│   ├── expected_animals.csv   every day-3 animal: in the foot / gill fits, or the reason it is out
│   └── _superseded/           gene-family count matrices written by the legacy script; not read
└── 03_analyses/
    ├── gene_mechanics/     scripts 01 and 02, tissue-suffixed
    │   └── candidate_scatter_F/, candidate_scatter_G/   script 01: one figure per heatmap gene
    ├── expr_tables/        script 03
    ├── byssus_genes/       script 04
    ├── go_mechanics/       script 05
    ├── expression_suites/  script 06 (README inside)
    ├── knit_html/          reports, logs and run_log.csv from runner 00 (git-ignored)
    └── _superseded/        the legacy script's two PDFs
```

Scripts 01 to 06 take a knit parameter `tissue` (`"F"` default, or `"G"`). Knit
`00_run_gene_mechanics_by_tissue.Rmd` to run the whole chain for both tissues, after the
`02_thread-strength`, `05_differential-expression` and `07_enrichment` runners (the
repository-level `00_run_pipeline.Rmd` runs them in order). Each script reads the previous one's CSV
handoffs; none shares an R session with another.

Foot means the phenol gland to the tip of the foot, the region sampled in every animal. The
twelve day-0 animals also have a library of the rest of the foot (IDs ending `FX`); no script
here uses them (see `05_differential-expression/README.md`).

## Inputs (cross-folder)

All paths resolve from a `repo_root` found by walking up from `here::here()` to the first
folder that contains both `02_thread-strength/` and `05_differential-expression/`.

| input | from | read by |
|---|---|---|
| `02_thread-strength/03_analyses/03_assemble-thread-summary/thread-summary.xlsx` | thread summary, 02 script 03 | 01, 03, 04 |
| `02_thread-strength/03_analyses/05_decompose-adhesion/mussel_response_classification.csv` | per-animal response, 02 script 05 | 01 |
| `02_thread-strength/03_analyses/0{4,5}_*/DATA_ancova_animals.csv` | per-animal ANCOVA rows, 02 scripts 04 and 05 | 01 (agreement check) |
| `02_thread-strength/03_analyses/02_extract-tensometer-data/thread-summary-raw-output.xlsx` | every extracted trace, 02 script 02 | 03 |
| `05_differential-expression/03_analyses/count_matrix/gene_count_matrix_clean.csv` | counts, 05 script 01 | 01, 03, 04 |
| `05_differential-expression/03_analyses/count_matrix/treatmentinfo_clean.csv` | treatment, day and region per library, 05 script 01 | 01, 03 |
| `05_differential-expression/03_analyses/DEG_lists/<Foot or Gill>/<T><X>_TC_siggene.csv` | TC DEG lists, 05 script 04 | 01, 03, 04 |
| `05_differential-expression/03_analyses/DEG_lists/GOterms_genome/*_sigs_ID.csv` | annotated TC DEGs, 05 script 06 | 01 |
| `05_differential-expression/03_analyses/count_matrix/mitochondrial_loci.csv` | the mitochondrial loci, left out of every gene universe here | 01 |
| `05_differential-expression/03_analyses/mitochondrial/mt_share_by_sample.csv` | mitochondrial protein reads as a share of each library's nuclear reads, 05 script 13 | 05 |
| `07_enrichment/03_analyses/02_topgo/topgo_enriched.csv` | enriched GO terms of each TC run and the DEGs in each, 07 script 02 | 05 |
| `07_enrichment/03_analyses/06_method-comparison/consensus_terms_TC_<ont>.csv` | terms an FDR-controlled method also enriched, 07 script 06 | 05 |
| `03_blast/03_analyses/genome-foot-sprot2026_03-noseg/LOC_GO_list.txt` | genome-wide BLAST (Swiss-Prot 2026_03 plus 241 foot and byssal proteins, and `03_blast` step 05's six genes found with the low-complexity filter off) / UniProt (release 2026_03 records) / GO annotation | 01, 04 |
| `tools/pipeline_checks.R` | `warn_unless()` checks and `RUN_provenance_<T>.txt` | 01, 02, 05, 06 |
| `tools/gene_ids.R` | `gene_key()`: count-matrix gene name to LOC key | 01, 03, 04, 05 |

## Analysis design

The reported test is one baseline-adjusted regression (ANCOVA) per gene and metric on the
per-animal values, `level_day3 ~ expression + treatment + level_baseline`, on the log scale.
The metrics are those of `02_thread-strength`: `mean_force` (the geometric mean of the
animal's thread peak forces) and `pad_area` are the declared primary metrics; `max_force`
(the animal's strongest thread, one thread per animal, so noisier) and `adhesion_kpa` are
exploratory (`METRICS` in script 01, `tier` column in every output). Extension is not
analysed (threads were cut near the plaque, so the length under test was not controlled).
Animals without baseline threads are not in the fits. BH is applied within each metric
(`q_lm`) and within each tier family (`q_family`). Candidates are the keyword matches in the
genome-wide BLAST annotation that are above the detection floor. The mitochondrial loci are
in no gene universe. Script 02 adds six pathway modules, the leave-one-animal-out influence
check and `best_hits_<T>.csv`, all through the same ANCOVA.

In the current run two genes pass q < 0.1 within their metric, and one of them within its tier
family too (`best_hits_<T>.csv`):

- **Gill, candidate, plaque area (primary).** Animals that express LOC134718612 more strongly
  have smaller day-3 plaques (partial r -0.56, 95% interval -0.74 to -0.31, p 1.1e-4, q 0.036
  over the 324 gill candidates tested, 0.072 over the primary family). It is annotated as heat shock 70 kDa
  protein 12A, but only from a 32% identity BLAST hit, so it is a distant HSP70-family gene
  rather than a confident ortholog; it is expressed at a low level in gill (DESeq2 base mean
  about 10) and is not a TC DEG. Leaving out its most influential animal (T136) gives p 1.6e-5.
  On the previous count matrix it was already the top plaque-area candidate (q 0.09, family q
  0.18).
- **Foot, DEG union, strongest thread (exploratory).** Animals that express LOC134698975
  (dynein axonemal assembly factor 8, 36% BLAST identity) more strongly make stronger threads (partial r 0.57, 95%
  interval 0.32 to 0.75, p 9.8e-5, q 0.091 over the 932 foot DEG-union genes, family q 0.18; mean
  peak force, the primary metric, r 0.56, p 1.6e-4, q 0.15). It is a Foot DO down-regulated
  DEG (log2 fold change -0.52) in the Foot DO down-regulated cilium and dynein-arm assembly terms (`07`); within arms,
  the animals with less of it make weaker threads. Robust to leaving out one animal.

The next candidates are glutathione peroxidase 7 with the strongest thread in gill
(exploratory, q 0.12); the lowest foot candidate q is 0.49 (262 foot candidates since 2026-10-06, when
`03_blast` step 05 annotated mefp-5 and two expressed mfp-3 genes: none comes near, p 0.39 or more), and
no module comes near (lowest q 0.61).
Until 2026-10-05 the gill's next result was aminoacylase-1 with adhesion (q 0.14), but it was
a candidate only because the keyword `aminoacyl`, meant for tRNA synthetases, matched its name;
the keyword is now `aminoacyl[- ]tRNA` (`01_code/gene-mechanics-pipeline_DOC.md`, caveats).

Script 05 puts the sets the expression analysis found through the same ANCOVA: the up- and
the down-regulated DEGs of each stressor, the genes behind every enriched topGO term of
those runs (all three ontologies, marked where an FDR-controlled method agrees), and the
mitochondrial protein reads as a share of the library. A set's score is the first principal
component of its genes, as for script 02's modules. BH runs in two families, the
response-level scores (DEG sets and mitochondrial share) and the GO terms. The sets were
chosen because they differ between arms and the model removes the arm differences, so a
slope is an association among animals of the same arm; with 10 to 12 animals per arm it is
exploratory. In the current run no set passes q < 0.1 (the lowest q is 0.19). The strongest
signals are in gill: animals with a higher mitochondrial share of reads make weaker threads
(strongest thread, partial r -0.34, p 0.028, q 0.19; mean peak force, r -0.33, p 0.034, q 0.24),
and animals scoring higher on the Gill OA up-regulated terms "respiratory chain complex" and
"proton motive force-driven ATP synthesis" have smaller plaques (pad area, partial r -0.38, p
0.012, and r -0.37, p 0.015; q 0.94 among the GO terms). Before the annotation and GO release
of 2026 (2026-10-04), the Gill OA up term "inner mitochondrial membrane protein complex" was
enriched and gave the lowest q (pad area, r -0.48, p 0.001, q 0.16); with the new records it is
no longer enriched, and five of its seven DEGs are among the six of "respiratory chain
complex". The mitochondrial share is counted on
the mitochondrial genome alone since 2026-10-02 (`05` step 13); with the previous genome count
the same two associations were r -0.33 and -0.35. Nothing in foot comes near (lowest q 0.58).

Script 06 asks whether combinations of genes ("suites") go with attachment, and separates the
two questions an "expression profile" plot mixes up: across arms (an association that can be
the treatment itself, since warming and hypoxia change both expression and thread strength)
and within arms (the question of scripts 01, 02 and 05). Its axes are script 05's DEG sets and
GO terms with at least 10 genes (the DEG sets have tens to hundreds of genes), co-expression
modules of the 4,000 genes with the largest within-arm variance (arm removed, signed
correlation, dynamic tree cut) and the first ten principal components of all expressed genes;
each gets the within-arm ANCOVA, the same model without the treatment term, the correlation of
its arm means with the arm means of the baseline-adjusted level and a test of differing slopes.
An elastic net then predicts held-out animals from the candidate genes, the 2,000 most variable
genes and all the axes, across and within arms, against label permutations, with the treatment
alone as a reference. It needs `glmnet` and `dynamicTreeCut` (in `renv.lock`) and runs its
prediction test on a socket cluster of `cores` workers (Windows too), about 7 minutes per tissue
on 4 workers with 100 permutations. Results: `03_analyses/expression_suites/README.md`.

Details in `01_code/gene-mechanics-pipeline_DOC.md`.
