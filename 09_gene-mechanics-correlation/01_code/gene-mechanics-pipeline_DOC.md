# Gene-mechanics correlation pipeline, scripts 01 to 05

Links foot or gill gene expression at day 3 to the same animal's byssal thread mechanics.
Five chained scripts, each reading the previous one's CSV handoffs rather than sharing an R
session, parameterised by tissue. Run in order **01 → 02 → 03 → 04 → 05**, after the
`02_thread-strength`, `05_differential-expression` and `07_enrichment` runners, or knit
**00_run_gene_mechanics_by_tissue.Rmd**, which renders all five for foot and gill. Script 05
(DEG sets, enriched GO terms and mitochondrial expression against mechanics) is described in
its own header and in the folder README.

## Running for a tissue

Each of 01 to 05 has a knit parameter in its YAML header:

```yaml
params:
  tissue: "F"   # "F" foot (default) or "G" gill
```

`TISSUE` is read from `params$tissue`, falling back to `"F"` when the script is run chunk by
chunk outside a knit. Runner 00 renders each script **in its own R process** (nested
`rmarkdown::render()` collides on knitr's chunk-label registry), through the shared
`tools/run_steps.R`, and writes the HTML reports, a log per script and a `run_log.csv` to
`03_analyses/knit_html/`, which is git-ignored.

Every output is tissue-suffixed (`_F` / `_G`) except `annotation_map.csv`, a tissue-independent
LOC-to-protein map. Treatment, day and sampled region per library come from the sample table
written by `05_differential-expression` (`treatmentinfo_clean.csv`), one row per library in the
clean count matrix. Foot means the phenol gland to the tip of the foot; the day-0 libraries of
the rest of the foot (IDs ending `FX`) are left out. An animal with a sample-table row but no
count-matrix column is dropped with a message rather than a hard stop.

```
02_thread-strength/03_analyses/03_assemble-thread-summary/
    thread-summary.xlsx                                             thread summary (script 03)
02_thread-strength/03_analyses/05_decompose-adhesion/
    mussel_response_classification.csv                              per-animal response (script 05)
02_thread-strength/03_analyses/02_extract-tensometer-data/
    thread-summary-raw-output.xlsx                                  every extracted trace (script 02)
05_differential-expression/03_analyses/count_matrix/
    gene_count_matrix_clean.csv                                     counts (script 01)
    treatmentinfo_clean.csv                                         arm, day, region per library (script 01)
05_differential-expression/03_analyses/DEG_lists/<Foot|Gill>/
    <T><X>_TC_siggene.csv                                           TC DEG lists (script 04)
        |
        v
01  paired table, VST, candidate set,
    per-gene ANCOVA (the reported test)             ->  03_analyses/gene_mechanics/
02  modules, diagnostics                            ->  03_analyses/gene_mechanics/
03  RNA x thread manifest, top-25 expression tables ->  03_analyses/expr_tables/
04  byssus/foot gene list + expression              ->  03_analyses/byssus_genes/
```

---

## 1. How the chain is set up

### Thread input and labels

Scripts 01, 03 and 04 read the thread summary `02_thread-strength/03_analyses/03_assemble-thread-summary/thread-summary.xlsx`. Pre-exposure threads are `phase == "pre"`, day-3 threads
`phase == "post"`, and the arm an animal was assigned to is `mussel_trt`; a
`read_thread_summary()` helper in each script accepts the older column spellings. Script 03
reconciles against 02_thread-strength script 02's extraction (`thread-summary-raw-output.xlsx`) and flags
animals sequenced at day 3 but never pulled (`no_post_trace`).

### The day-3 control arm is in the paired set

`INCLUDE_CONTROL_ARM <- TRUE` in script 01. The day-3 control animals (T126-T137) have
foot RNA, day-3 threads and their own baselines, and enter as a fourth treatment level;
the control animals' own trajectory is the null an expression signal must beat. An
association that holds within the control arm too is about attachment biology, not the
stress response. Set FALSE for stressor arms only.

### Metrics, tiers and the reported test

Peak force and plaque area are analysed separately as well as through their ratio
(adhesion), since a response in either component can be diluted in the ratio
(`02_thread-strength/01_code/05_decompose_adhesion_DOC.md` decomposes adhesion the same way).
Script 01 tests, per gene and metric, one baseline-adjusted regression (ANCOVA) on the
per-animal values:

    level_day3 ~ expression + treatment + level_baseline

| metric | scale | tier | per-animal value |
|---|---|---|---|
| `mean_force` | log | primary | geometric mean of the peak forces of the animal's day-3 threads, N |
| `pad_area` | log | primary | geometric mean, mm² |
| `max_force` | log | exploratory | the largest thread peak force, N (one thread per animal, so noisier; a maximum also grows with the number of threads) |
| `adhesion_kpa` | log | exploratory | geometric mean of force / area × 1000 (recomputed per plaque) |

Every metric enters on the log scale on both sides of the model, so `slope` is a log-unit
change per VST unit and exp(slope) a multiplicative one. Extension is not analysed: threads
were cut near the junction of the plaque and the distal region, so the length of distal
thread under test differed between pulls (`02_thread-strength/README.md`). `level_baseline` is the same summary of the animal's own pre-exposure threads, on the
same scale. Animals without baseline threads are not in the fits. `baseline_slope` is
reported beside `slope`; it is the same quantity as the baseline coefficient of the
02_thread-strength ANCOVA (`STATS_ancova_coefficients.csv` in 02_thread-strength scripts 3 and 4),
fitted here with expression added.

`METRICS` in script 01 fixes the scale and the tier; `tier = primary` (mean force, area) is the
declared confirmatory family and every output carries the column. Multiplicity: `q_lm` is BH
within a gene set x metric, `q_family` BH within a gene set x tier, so the candidate x
primary family (tested candidates x two primary metrics) has its own search-corrected q.
`partial_r` expresses the same test as a partial correlation, t / sqrt(t^2 + residual df),
with the sign and p of `slope`. `metrics_config_<T>.csv` is
written by 01 and read by 02, so the metric list, scales, tiers, arm levels and covariates
cannot drift between the two.

Why an ANCOVA and not a change score: the change score `log(day3) - log(baseline)` imposes a
baseline coefficient of 1; when the fitted coefficient is well below 1, the change score
adds most of the baseline's measurement noise to the response and loses power. The ANCOVA
also absorbs any between-animal baseline differences the arm assignment did not balance
(baseline balance by future arm is checked in 02_thread-strength script 3, 3b).

### Annotation map and candidate universe

`CANDIDATE_ANNOTATION = "genome"`: every gene in the count matrix is annotated with its best
UniProt hit (highest bitscore) from `03_blast/03_analyses/genome-foot-sprot2026_03-noseg/LOC_GO_list.txt`, with
`blast_pident` and `blast_evalue` carried along, and any expressed gene whose name matches
`CANDIDATE_KEYWORDS` (byssal / collagen / plaque-curing / HSP / hypoxia / tRNA-synthetase /
oxidative-stress terms) and passes the BLAST floor (`CANDIDATE_MAX_EVALUE = 1e-10`,
`CANDIDATE_MIN_PIDENT = 0`) is a candidate. The alternative `"TC_DEG"` restricts candidates
to genes named in the treatment-vs-control DEG tables plus the byssal structural genes; it
made "already a DEG in some contrast" a hidden entry condition and is kept only for comparison. `in_TC_DEG_annotation` marks the
overlap in every table.

Several best hits of the BLAST search of 2026 are *M. coruscus* byssus proteins (Qin et al.
2016, J Proteomics 144:87-98) whose names carry no byssal word. The keywords name them
explicitly: "YGH-rich protein" as byssal structural (the genes whose 2024 best hit was Foot
protein 12), and "protease inhibitor-like protein-1", "C1q-domain-containing protein-1" and
"TSP_1 domain containing protein-1" as byssal accessory (candidates in script 01, the
`byssal_collagen` module in script 02, `byssal_accessory` in script 04). These names occur
only on *M. coruscus* entries in the BLAST tables.

### Gene keys and the mitochondrial loci

Count-matrix gene names (`gene-LOC134696364|LOC134696364`; `STRG.10|LOC...` in the previous matrix) become LOC keys
through `gene_key()` in `tools/gene_ids.R`, the same function 05 and 07 use; script 01 stops
if two tested genes share a key. The 331 mitochondrial loci of
`05_differential-expression/03_analyses/count_matrix/mitochondrial_loci.csv` (the
mitochondrial genome's genes and their copies on unplaced scaffolds) are removed after the
expression filter, so they enter neither the candidate set nor the DEG union; script 05
tests their summed share of the library as one score.

### Detection floor

A gene at the VST floor (zero counts) in many paired animals gives an association driven by
presence/absence. Script 01 computes `frac_at_floor` for every tested gene, flags `exclude`
(> 0.40) and `caution` (> 0.20), writes `detection_floor_flags_<T>.csv` for the candidate
set and the DEG union, and with `EXCLUDE_FLOOR_GENES = TRUE` (default) drops the `exclude`
genes from both families before testing; script 01 reports how many candidate and
DEG-union genes are tested. Excluded genes still contribute to the module scores in
script 02.

### Script 02: modules, diagnostics

Script 02 reads script 01's association tables and does not re-derive the reported test.
Every model it fits is that same per-animal ANCOVA. It adds:

- **Block A, module eigengenes** (`module_associations_<T>.csv`, `module_members_<T>.csv`):
  PC1 of each module's members (genes passing the BLAST floor, `blast_ok`), oriented so a
  higher score is higher expression, through the same ANCOVA (`p_lm`, `q_lm` across the six
  modules within a metric, `q_family` within a tier).
  `byssal_structural` is a sixth module (foot proteins, preCols including preCOL-NG, named
  "Nongradient byssal" in UniProt, the thread matrix proteins, the *M. coruscus* YGH-rich
  proteins (Qin et al. 2016; the best hits, in the search of 2026, of the genes whose 2024 best
  hit was Foot protein 12), byssal EP/ACDC, the
  plaque-curing tyrosinase: the structural proteins of the plaque and thread), separate from the broad `byssal_collagen` regex, so the test "these genes track
  thread building, not strength" has its own row (`BYSSAL_STRUCTURAL_REGEX` in script 01
  flags the same genes in the candidate table).
- **Block B, diagnostics**, below.

### Script 02 diagnostics (block B)

- `detection_floor_flags_<T>.csv` (written by script 01, echoed here): one row per gene in
  the candidate set and the DEG union, with `frac_at_floor`, `floor_flag` and `tested`.
- `influence_top_hits_<T>.csv`: the three best candidate hits and the best DEG-union hit per
  metric (`N_INFLUENCE_HITS`, ranked by `p_lm`), each refitted as the reported ANCOVA with
  and without its most influential animal. Columns: `most_influential_mussel`,
  `max_cooks_D`, `cook_threshold_4n`, `cooks_flag` (`D>1`, `D>4/n`, `ok`), `slope`,
  `slope_without_mussel`, `p_without_mussel`, `slope_change_frac`, and `influence_flag` =
  `fragile` when the hit loses p < 0.05 without that animal or its slope moves by more than
  half, else `robust`. The maximum of ~45 Cook's distances is expected to exceed 4/n in
  most fits, so `influence_flag` is the column to read.
- `best_hits_<T>.csv`: the single best hit per metric for each gene set (candidate, module,
  DEG union) with `p_lm`, `q_lm`, `q_family`, the floor flag and the influence flag on one
  row. This is the table to quote from.

### Bioconductor masking

`S4Vectors` and `IRanges`, loaded by DESeq2, mask `dplyr::rename`, `count`, `first` and
`desc`. Script 01 loads DESeq2 first and tidyverse last, and uses `dplyr::` prefixes.
Loading tidyverse before DESeq2 fails with `object 'treatment' not found`.

---

## 2. Inputs from 02_thread-strength

| file | produced by | used by |
|---|---|---|
| `03_analyses/03_assemble-thread-summary/thread-summary.xlsx` | script 03 | 01, 03, 04 |
| `03_analyses/05_decompose-adhesion/mussel_response_classification.csv` | script 05 | 01 (joined into the paired manifest) |
| `03_analyses/02_extract-tensometer-data/thread-summary-raw-output.xlsx` | script 02 | 03 |

The response classification carries, per animal and per metric, the pre and post means, the
log-ratio, the % change, the raw direction (`decreased` / `increased`), the change relative
to the control arm's mean change, and a composite `response_class` (`weaker` if both force
and adhesion fell, `stronger` if both rose, else `mixed`) and `response_score` (mean
standardised log-ratio across force, area and adhesion). Script 01 joins the class and the
score into `paired_sample_manifest_<T>.csv` for inspection; they enter no model.

---

## 3. Outputs

### `03_analyses/gene_mechanics/` (scripts 01 and 02)

| file | contents |
|---|---|
| `paired_sample_manifest_<T>.csv` | the paired animals: arm, plaque counts, day-3 and baseline per-animal values (geometric means for `mean_force`, area and adhesion; the largest thread for `max_force`), response class and score |
| `metrics_config_<T>.csv` | metric list with `scale` (log / raw), `tier` (primary / exploratory) and label, arm levels, covariates, the `EXCLUDE_FLOOR_GENES` setting |
| `vst_paired_<T>.csv` | handoff: VST expression of the paired samples |
| `annotation_map.csv` | genome-wide best UniProt hit per LOC with `blast_pident`, `blast_evalue`, `blast_ok`, `in_TC_DEG_annotation` |
| `candidate_genes_<T>.csv` | the candidate set with `byssal_structural`, `frac_at_floor`, `floor_flag`, `tested` |
| `detection_floor_flags_<T>.csv` | every candidate and DEG-union gene: fraction of paired samples at the VST floor, `floor_flag`, `tested` |
| `assoc_candidate_<T>.csv`, `assoc_DEGunion_<T>.csv` | **the reported test**: script 01 ANCOVA per gene x metric with `scale`, `tier`, `n`, `slope`, `se`, `p_lm`, `baseline_slope`, `partial_r` with its marginal 95% interval (`partial_r_lo`, `partial_r_hi`), `q_lm`, `q_family`, floor flag |
| `module_associations_<T>.csv`, `module_members_<T>.csv` | six pathway modules x four metrics (ANCOVA); the member genes |
| `best_hits_<T>.csv` | best hit per metric and gene set with `p_lm`, `q_lm`, `q_family`, floor and influence flags |
| `influence_top_hits_<T>.csv` | top three candidate hits and the best DEG-union hit per metric with leave-one-out slope and p, `influence_flag`, `q_lm`, `q_family` |
| `RUN_provenance_<T>.txt` | settings of scripts 01 and 02 (arms, covariates, model, metrics with scale and tier, modules, family sizes), the code commit that ran and whether tracked files differed from it, R and package versions, and an MD5 of every input, all with repository-relative paths |
| `animal_reconciliation_<T>.csv` | animals whose presence in the fits differs from `02_data/expected_animals.csv` (empty when they agree) |
| `candidate_heatmap_<T>.png`, `top_candidate_scatter_<T>.png`, `best_hit_per_metric_scatter_<T>.png` | figures of genes selected by smallest p (effects biased away from zero). The two scatter files are added-variable plots: day-3 level and expression each residualised on arm and baseline, with the tested slope drawn through the origin, points coloured by arm |
| `candidate_scatter_<T>/NN_<LOC>_<name>.png` | one figure per heatmap gene (60 per tissue), numbered in the heatmap's row order: an added-variable panel per metric, each with partial r and its 95% interval, p, `q_lm` and `q_family`. Rewritten on every run (the folder's figures are removed first) |

### `03_analyses/expr_tables/` (script 03) and `03_analyses/byssus_genes/` (script 04)

`rna_thread_manifest_<T>.csv` is tissue-suffixed. `byssus_gene_expression_<T>.csv` lists the
byssal and foot genes (by category) with more than 5 reads in at least a third of the
thread-having animals, and every mussel foot protein gene with reads in the tissue even below
that filter (`mfp` TRUE, `expressed` FALSE; 12 in the foot, among them mfp-6, one mfp-3 and three
mfp-1 copies), so that no mfp gene is left out of the table; `byssus_category_scores_<T>.csv`
averages only the genes that pass the filter. The companion `sample_metadata_<T>.csv`
files carry all four arms and the per-animal thread values (arithmetic mean of the thread
peak forces as `mean_force`, the largest as `max_force`, mean area and adhesion), a
description rather than a model input.

---

## Checks

Script 01 checks, before any model is fitted:

- the animals in the fits equal `02_data/expected_animals.csv` for the tissue (every other
  day-3 animal is listed there with the reason it is out); a difference is written to
  `animal_reconciliation_<T>.csv`;
- each animal's day-3 and baseline values equal the ones 02_thread-strength's ANCOVA used
  (`DATA_ancova_animals.csv`), so the two pipelines cannot drift apart;
- every expected 05_differential-expression input exists and reads (one `<T><X>_TC_siggene.csv`
  file per stressor, six `*_sigs_ID.csv` files), and the arm in the thread key agrees with the
  Tag-seq sample table.

A failed check does not stop the run (`warn_unless()` in `tools/pipeline_checks.R`). It
prints `CHECK FAILED: <what failed>` in the knitted report and in the render log (runner 00's
`.log` files), and `RUN_provenance_<T>.txt` records how many checks ran and lists each one that
failed. Read that line before using a run.
Update `expected_animals.csv` only with a documented reason. Runner 00 stops with an error,
after writing `run_log.csv`, if any step failed.

## 4. Known data quirks

- **T047** has day-3 threads and a gill library (`T047G`) but no foot library: there is no
  `T047F` column in the count matrix and no foot row in `treatmentinfo_clean.csv`. It can enter
  the gill paired set, not the foot one.
- **Foot region.** Every foot library used here is the phenol gland to the tip of the foot.
  In day-3 animals several byssal plaque genes are at or near the VST floor (for example
  LOC134711106, foot protein-4 variant-1, and LOC134692428, byssal peroxidase-like 4: zero
  counts in 28 and 27 of the 46 day-3 foot libraries, against medians of 159 and 215 counts in
  the 12 day-0 ones), so the detection-floor filter removes part of the byssal structural
  family: in the foot run 2 of the 20 `byssal_structural` candidates are excluded and 5 more
  are flagged `caution` (`candidate_genes_F.csv`; 2 of 18 and 5 with the 2024 BLAST search;
  with the previous StringTie + prepDE counts, 6 of 16 and 5). A null result for those genes is not evidence of no association.
- The candidate keywords are regexes on UniProt names, matched without regard to case, so a
  short keyword can match inside an unrelated name. `aminoacyl` matched aminoacylase-1 and
  acylaminoacyl-peptidase until 2026-10-05; it is now `aminoacyl[- ]tRNA` (scripts 01 and 02).
  `Hsp` still matches abbreviations inside names ("HSPG", "HsPDE8B", "hSPL", "HSPK 21",
  "CRHSP-24"), which brings in 6 foot and 9 gill candidates with no heat-shock role (perlecan,
  phosphodiesterase 8B, sphingosine-1-phosphate lyase and phosphatase, Nek2, and others) and
  PERK, an ER-stress kinase, by its abbreviation "HsPEK"; and `chaperone` pulls in histone
  and assembly chaperones, so the `HSP_proteostasis` module is broad. Tighten
  `CANDIDATE_KEYWORDS` or raise `CANDIDATE_MIN_PIDENT` if a narrower family is wanted.

---

## 5. Configuration summary

| script | option | default | effect |
|---|---|---|---|
| 01 | `INCLUDE_CONTROL_ARM` | TRUE | day-3 control animals as a fourth arm |
| 01 | `METRICS` | mean force, area (primary); maximum force, adhesion (exploratory); all log | metric, thread-level column, per-animal summary (`agg`: geometric mean or largest thread), model scale and tier; reporting order = priority; `PRIMARY_METRICS` is derived from it |
| 01 | `CANDIDATE_ANNOTATION` | "genome" | best genome-wide BLAST hit per LOC ("TC_DEG": DEG-table names only) |
| 01 | `CANDIDATE_MAX_EVALUE`, `CANDIDATE_MIN_PIDENT` | 1e-10, 0 | BLAST-quality floor for candidates and module members |
| 01 | `EXCLUDE_FLOOR_GENES`, `FLOOR_EXCLUDE`, `FLOOR_CAUTION` | TRUE, 0.40, 0.20 | detection-floor filter applied before testing |
| 01 | `BYSSAL_STRUCTURAL_REGEX` | foot protein, preCol, ACDC, ... | flags the byssal structural genes |
| 01, 02 | `FDR_ALPHA` | 0.10 | BH threshold for flagging (`q_lm`, `q_family`) |
| 02 | `N_INFLUENCE_HITS` | 3 | candidate hits per metric that get the leave-one-animal-out refit |
| 02 | `modules` | six regexes | includes `byssal_structural` |
| 03 | `USE_RAW_THREAD_SET` | TRUE | any extracted trace vs thread-summary mussels only |
| 04 | `USE_RAW_THREAD_SET` | FALSE | |
| 01-04 | `params$tissue` | "F" | foot or gill |
| 00 | `params$tissues`, `params$steps`, `params$stop_on_fail` | both, "all", TRUE | what the runner renders |
