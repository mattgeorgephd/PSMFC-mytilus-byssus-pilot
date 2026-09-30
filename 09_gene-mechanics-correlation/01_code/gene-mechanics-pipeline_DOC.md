# Gene-mechanics correlation pipeline, scripts 20 to 23

Links foot or gill gene expression at day 3 to the same animal's byssal thread mechanics.
Four chained scripts, each reading the previous one's CSV handoffs rather than sharing an R
session, parameterised by tissue. Run in order **20 → 21 → 22 → 23**, after thread-strength
scripts 1 to 4, or knit **24-run_gene_mechanics_by_tissue.Rmd**, which renders all four for
foot and gill.

## Running for a tissue

Each of 20 to 23 has a knit parameter in its YAML header:

```yaml
params:
  tissue: "F"   # "F" foot (default) or "G" gill
```

`TISSUE` is read from `params$tissue`, falling back to `"F"` when the script is run chunk by
chunk outside a knit. Driver 24 renders each script **in its own R process** (nested
`rmarkdown::render()` collides on knitr's chunk-label registry) and writes the HTML reports
and a `run_log.csv` to `03_analyses/knit_html/`, which is git-ignored.

Every output is tissue-suffixed (`_F` / `_G`) except `annotation_map.csv`, a tissue-independent
LOC-to-protein map. The gill treatmentinfo codes the day-3 control arm as `control_3` and has no
leading index column; scripts 20 and 22 normalise both. An animal with a treatmentinfo row but no
count-matrix column is dropped with a message rather than a hard stop.

```
thread-strength/03_analyses/thread-summary.xlsx                     curated threads (scripts 1-2)
thread-strength/03_analyses/04_decompose-adhesion/
    mussel_response_classification.csv                              per-animal response (script 4)
thread-strength/03_analyses/01_extract-tensometer-data/
    thread-summary-raw-output.xlsx                                  every extracted trace (script 1)
differential-expression/02_data/gene_count_matrix_clean.csv         counts
differential-expression/03_analyses/DEG_lists/Foot/F_treatmentinfo.csv   Tag-seq arm per sample
        |
        v
20  paired table, VST, candidate set,
    per-gene ANCOVA (the reported test)             ->  03_analyses/gene_mechanics/
21  modules, diagnostics                          ->  03_analyses/gene_mechanics/
22  RNA x thread manifest, top-25 expression tables ->  03_analyses/expr_tables/
23  byssus/foot gene list + expression              ->  03_analyses/byssus_genes/
```

---

## 1. How the chain is set up

### Thread input and labels

Scripts 20, 22 and 23 read the curated thread table `thread-strength/03_analyses/thread-summary.xlsx`. Pre-exposure threads are `phase == "pre"`, day-3 threads
`phase == "post"`, and the arm an animal was assigned to is `mussel_trt`; a
`read_thread_summary()` helper in each script accepts the older column spellings. Script 22
reconciles against script 1's extraction (`thread-summary-raw-output.xlsx`) and flags
animals sequenced at day 3 but never pulled (`no_post_trace`).

### The day-3 control arm is in the paired set

`INCLUDE_CONTROL_ARM <- TRUE` in script 20. The day-3 control animals (T126-T137) have
foot RNA, day-3 threads and their own baselines, and enter as a fourth treatment level;
the control animals' own trajectory is the null an expression signal must beat. An
association that holds within the control arm too is about attachment biology, not the
stress response. Set FALSE for stressor arms only.

### Metrics, tiers and the reported test

Peak force and plaque area are analysed separately as well as through their ratio
(adhesion), since a response in either component can be diluted in the ratio
(`thread-strength/01_code/4_decompose_adhesion_DOC.md` decomposes adhesion the same way).
Script 20 tests, per gene and metric, one baseline-adjusted regression (ANCOVA) on the
per-animal values:

    level_day3 ~ expression + treatment + level_baseline

| metric | scale | tier | per-animal value |
|---|---|---|---|
| `max_force` | log | primary | geometric mean of the animal's day-3 plaques, N |
| `pad_area` | log | primary | geometric mean, mm² |
| `adhesion_kpa` | log | exploratory | geometric mean of force / area × 1000 (recomputed per plaque) |
| `max_displacement` | raw | exploratory | arithmetic mean of extension at break, mm |

Force, area and adhesion enter as log(geometric mean) on both sides of the model, so
`slope` is a log-unit change per VST unit and exp(slope) a multiplicative one; extension is
raw. `level_baseline` is the same summary of the animal's own pre-exposure threads, on the
same scale. Animals without baseline threads are not in the fits. `baseline_slope` is
reported beside `slope`; it is the same quantity as the baseline coefficient of the
thread-strength ANCOVA (`STATS_ancova_coefficients.csv` in thread-strength scripts 3 and 4),
fitted here with expression added.

`METRICS` in script 20 fixes the scale and the tier; `tier = primary` (force, area) is the
declared confirmatory family and every output carries the column. Multiplicity: `q_lm` is BH
within a gene set x metric, `q_family` BH within a gene set x tier, so the candidate x
primary family (tested candidates x two primary metrics) has its own search-corrected q.
`partial_r` expresses the same test as a partial correlation, t / sqrt(t^2 + residual df),
with the sign and p of `slope`. `metrics_config_<T>.csv` is
written by 20 and read by 21, so the metric list, scales, tiers, arm levels and covariates
cannot drift between the two.

Why an ANCOVA and not a change score: the change score `log(day3) - log(baseline)` imposes a
baseline coefficient of 1; when the fitted coefficient is well below 1, the change score
adds most of the baseline's measurement noise to the response and loses power. The ANCOVA
also absorbs any between-animal baseline differences the arm assignment did not balance
(baseline balance by future arm is checked in thread-strength script 3, 3b).

### Annotation map and candidate universe

`CANDIDATE_ANNOTATION = "genome"`: every gene in the count matrix is annotated with its best
UniProt hit (highest bitscore) from `blast/03_analyses/genome-foot/LOC_GO_list.txt`, with
`blast_pident` and `blast_evalue` carried along, and any expressed gene whose name matches
`CANDIDATE_KEYWORDS` (byssal / collagen / plaque-curing / HSP / hypoxia / tRNA-synthetase /
oxidative-stress terms) and passes the BLAST floor (`CANDIDATE_MAX_EVALUE = 1e-10`,
`CANDIDATE_MIN_PIDENT = 0`) is a candidate. The alternative `"TC_DEG"` restricts candidates
to genes named in the treatment-vs-control DEG tables plus the byssal structural genes; it
made "already a DEG in some contrast" a hidden entry condition and is kept only for comparison. `in_TC_DEG_annotation` marks the
overlap in every table.

### Detection floor

A gene at the VST floor (zero counts) in many paired animals gives an association driven by
presence/absence. Script 20 computes `frac_at_floor` for every tested gene, flags `exclude`
(> 0.40) and `caution` (> 0.20), writes `detection_floor_flags_<T>.csv` for the candidate
set and the DEG union, and with `EXCLUDE_FLOOR_GENES = TRUE` (default) drops the `exclude`
genes from both families before testing; script 20 reports how many candidate and
DEG-union genes are tested. Excluded genes still contribute to the module scores in
script 21.

### Script 21: modules, diagnostics

Script 21 reads script 20's association tables and does not re-derive the reported test.
Every model it fits is that same per-animal ANCOVA. It adds:

- **Block A, module eigengenes** (`module_associations_<T>.csv`, `module_members_<T>.csv`):
  PC1 of each module's members (genes passing the BLAST floor, `blast_ok`), oriented so a
  higher score is higher expression, through the same ANCOVA (`p_lm`, `q_lm` across the six
  modules within a metric, `q_family` within a tier).
  `byssal_structural` is a sixth module (foot proteins, preCols, byssal EP/ACDC, the
  plaque-curing tyrosinase: the structural proteins of the plaque and thread), separate from the broad `byssal_collagen` regex, so the test "these genes track
  thread building, not strength" has its own row (`BYSSAL_STRUCTURAL_REGEX` in script 20
  flags the same genes in the candidate table).
- **Block B, diagnostics**, below.

### Script 21 diagnostics (block B)

- `detection_floor_flags_<T>.csv` (written by script 20, echoed here): one row per gene in
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
`desc`. Script 20 loads DESeq2 first and tidyverse last, and uses `dplyr::` prefixes.
Loading tidyverse before DESeq2 fails with `object 'treatment' not found`.

---

## 2. Inputs from thread-strength

| file | produced by | used by |
|---|---|---|
| `03_analyses/thread-summary.xlsx` | scripts 1–2 | 20, 22, 23 |
| `03_analyses/04_decompose-adhesion/mussel_response_classification.csv` | script 4 | 20 (joined into the paired manifest) |
| `03_analyses/01_extract-tensometer-data/thread-summary-raw-output.xlsx` | script 1 | 22 |

The response classification carries, per animal and per metric, the pre and post means, the
log-ratio, the % change, the raw direction (`decreased` / `increased`), the change relative
to the control arm's mean change, and a composite `response_class` (`weaker` if both force
and adhesion fell, `stronger` if both rose, else `mixed`) and `response_score` (mean
standardised log-ratio across force, area and adhesion). Script 20 joins the class and the
score into `paired_sample_manifest_<T>.csv` for inspection; they enter no model.

---

## 3. Outputs

### `03_analyses/gene_mechanics/` (scripts 20 and 21)

| file | contents |
|---|---|
| `paired_sample_manifest_<T>.csv` | the paired animals: arm, plaque counts, day-3 and baseline per-animal values (geometric means for force, area, adhesion), response class and score |
| `metrics_config_<T>.csv` | metric list with `scale` (log / raw), `tier` (primary / exploratory) and label, arm levels, covariates, the `EXCLUDE_FLOOR_GENES` setting |
| `vst_paired_<T>.csv` | handoff: VST expression of the paired samples |
| `annotation_map.csv` | genome-wide best UniProt hit per LOC with `blast_pident`, `blast_evalue`, `blast_ok`, `in_TC_DEG_annotation` |
| `candidate_genes_<T>.csv` | the candidate set with `byssal_structural`, `frac_at_floor`, `floor_flag`, `tested` |
| `detection_floor_flags_<T>.csv` | every candidate and DEG-union gene: fraction of paired samples at the VST floor, `floor_flag`, `tested` |
| `assoc_candidate_<T>.csv`, `assoc_DEGunion_<T>.csv` | **the reported test**: script 20 ANCOVA per gene x metric with `scale`, `tier`, `n`, `slope`, `se`, `p_lm`, `baseline_slope`, `partial_r` with its marginal 95% interval (`partial_r_lo`, `partial_r_hi`), `q_lm`, `q_family`, floor flag |
| `module_associations_<T>.csv`, `module_members_<T>.csv` | six pathway modules x four metrics (ANCOVA); the member genes |
| `best_hits_<T>.csv` | best hit per metric and gene set with `p_lm`, `q_lm`, `q_family`, floor and influence flags |
| `influence_top_hits_<T>.csv` | top three candidate hits and the best DEG-union hit per metric with leave-one-out slope and p, `influence_flag`, `q_lm`, `q_family` |
| `RUN_provenance_<T>.txt` | settings of scripts 20 and 21 (arms, covariates, model, metrics with scale and tier, modules, family sizes), the code commit that ran and whether tracked files differed from it, R and package versions, and an MD5 of every input, all with repository-relative paths |
| `animal_reconciliation_<T>.csv` | animals whose presence in the fits differs from `02_data/expected_animals.csv` (empty when they agree) |
| `candidate_heatmap_<T>.png`, `top_candidate_scatter_<T>.png`, `best_hit_per_metric_scatter_<T>.png` | figures of genes selected by smallest p (effects biased away from zero). The two scatter files are added-variable plots: day-3 level and expression each residualised on arm and baseline, with the tested slope drawn through the origin, points coloured by arm |

### `03_analyses/expr_tables/` (script 22) and `03_analyses/byssus_genes/` (script 23)

`rna_thread_manifest_<T>.csv` is tissue-suffixed. The companion `sample_metadata_<T>.csv`
files carry all four arms and `max_displacement`.

---

## Checks

Script 20 checks, before any model is fitted:

- the animals in the fits equal `02_data/expected_animals.csv` for the tissue (every other
  day-3 animal is listed there with the reason it is out); a difference is written to
  `animal_reconciliation_<T>.csv`;
- each animal's day-3 and baseline values equal the ones thread-strength's ANCOVA used
  (`DATA_ancova_animals.csv`), so the two pipelines cannot drift apart;
- every expected differential-expression input exists and reads (one `*_TC_siggene*` file per
  stressor, six `*_sigs_ID.csv` files), and the arm in the thread key agrees with the Tag-seq
  treatment table.

A failed check does not stop the run (`warn_unless()` in `tools/pipeline_checks.R`). It
prints `CHECK FAILED: <what failed>` in the knitted report and in the render log (driver 24's
`.log` files), and `RUN_provenance_<T>.txt` records how many checks ran and lists each one that
failed. Read that line before using a run.
Update `expected_animals.csv` only with a documented reason. Driver 24 stops with an error,
after writing `run_log.csv`, if any step failed.

## 4. Known data quirks

- **T047** has day-3 threads and a gill library (`T047G`) but no foot library: there is no
  `T047F` column in the count matrix and no row in `F_treatmentinfo.csv`. It can enter the
  gill paired set, not the foot one.
- The candidate keywords are regexes on UniProt names; `Hsp` and `chaperone` in particular
  pull in co-chaperones and assembly factors, so the `HSP_proteostasis` module is broad.
  Tighten `CANDIDATE_KEYWORDS` or raise `CANDIDATE_MIN_PIDENT` if a narrower family is
  wanted.

---

## 5. Configuration summary

| script | option | default | effect |
|---|---|---|---|
| 20 | `INCLUDE_CONTROL_ARM` | TRUE | day-3 control animals as a fourth arm |
| 20 | `METRICS` | force, area (log, primary); adhesion (log), extension (raw), exploratory | metric, model scale and tier; reporting order = priority; `PRIMARY_METRICS` is derived from it |
| 20 | `CANDIDATE_ANNOTATION` | "genome" | best genome-wide BLAST hit per LOC ("TC_DEG": DEG-table names only) |
| 20 | `CANDIDATE_MAX_EVALUE`, `CANDIDATE_MIN_PIDENT` | 1e-10, 0 | BLAST-quality floor for candidates and module members |
| 20 | `EXCLUDE_FLOOR_GENES`, `FLOOR_EXCLUDE`, `FLOOR_CAUTION` | TRUE, 0.40, 0.20 | detection-floor filter applied before testing |
| 20 | `BYSSAL_STRUCTURAL_REGEX` | foot protein, preCol, ACDC, ... | flags the byssal structural genes |
| 20, 21 | `FDR_ALPHA` | 0.10 | BH threshold for flagging (`q_lm`, `q_family`) |
| 21 | `N_INFLUENCE_HITS` | 3 | candidate hits per metric that get the leave-one-animal-out refit |
| 21 | `modules` | six regexes | includes `byssal_structural` |
| 22 | `USE_RAW_THREAD_SET` | TRUE | any extracted trace vs curated only |
| 23 | `USE_RAW_THREAD_SET` | FALSE | |
| 20-23 | `params$tissue` | "F" | foot or gill |
| 24 | `params$tissues`, `params$scripts` | both, all four | what the driver renders |
