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
thread-strength/03_analyses/decompose-adhesion/
    mussel_response_classification.csv                              per-animal response (script 4)
thread-strength/03_analyses/extract-tensometer-data/
    thread-summary-raw-output.xlsx                                  every extracted trace (script 1)
differential-expression/02_data/gene_count_matrix_clean.csv         counts
differential-expression/03_analyses/DEG_lists/Foot/F_treatmentinfo.csv   Tag-seq arm per sample
        |
        v
20  paired table, VST, candidate set, per-gene lm  ->  03_analyses/gene_mechanics/
21  mixed models, weighted regression, modules,
    permutation, diagnostics                        ->  03_analyses/gene_mechanics/
22  RNA x thread manifest, top-25 expression tables ->  03_analyses/expr_tables/
23  byssus/foot gene list + expression              ->  03_analyses/byssus_genes/
```

---

## 1. What changed on 16 September 2026, and why

### The input was gone

All of 20, 22 and 23 read `thread-strength/02_data/thread-summary-trossulus-clean.xlsx`,
deleted in commit `91dea50`, and 22 scanned `tensometer_output/control/` and `treatment/`,
which no longer exist. None of the four scripts could run.

### The `treatment == "control"` overload is gone

The old file's `treatment` column meant "baseline" for pre-exposure threads and "control"
for the day-3 control arm at the same time. The scripts filtered `treatment == "control"` to
get baselines. That is now `phase == "pre"`; day-3 threads are `phase == "post"`; the arm an
animal was assigned to is `mussel_trt`. A `read_thread_summary()` helper in each script
accepts the older column spellings so an old curated file still joins.

### The day-3 control arm is in the paired set

`INCLUDE_CONTROL_ARM <- TRUE` in script 20. The ten common-garden control animals
(T126–T135) have foot RNA, day-3 threads and, since the plaque areas were completed, their
own baselines. They enter as a fourth treatment level. This does two things: n rises from
31 to **44** (42 with baselines, up from 27), and the control animals' own ~20% decline
becomes the null trajectory an expression signal must beat. An association that holds
within the control arm too is about attachment biology, not the stress response.
Set FALSE for the pre-2026-09 behaviour (stressor arms only).

### Force and area are primary; change-from-baseline metrics are first-class

On the full thread dataset the stressor effect is in peak force and plaque area separately,
not in their ratio (`thread-strength/01_code/4_decompose_adhesion_DOC.md`). Script 20 now
tests:

| metric | type | what it is |
|---|---|---|
| `max_force` | level | per-animal mean of its day-3 plaques, N |
| `pad_area` | level | mm² |
| `adhesion_kpa` | level | force / area × 1000, recomputed from the two |
| `max_displacement` | level | extension at break, mm; new |
| `dlog_max_force` | change | log(day-3 mean / baseline mean), per animal |
| `dlog_pad_area` | change | as above |
| `dlog_adhesion_kpa` | change | as above; equals `dlog_max_force − dlog_pad_area` |

`ADJUST_FOR_BASELINE` applies only to the level metrics (ANCOVA on the day-3 level with the
animal's baseline as covariate). The change metrics need no covariate; the baseline is inside
the response. The two answer related but different questions and can disagree when baseline
and change are correlated.

`metrics_config.csv` is written by 20 and read by 21, so the metric list and the arm levels
cannot drift between the two.

### Script 21: one reported test, permuted for everything (17 September 2026)

The reported statistic for **every** metric is the per-animal precision-weighted regression
(block B; weight = day-3 plaque count, or min(baseline, day-3) for a change metric), because
it is the statistic the permutation calibrates. The thread-level mixed model (block A) stays
as a sensitivity column (`p_mixed`) for the level metrics; earlier versions ranked the level
metrics by `p_mixed` and permuted the weighted regression, so a reported p and its
permutation companion referred to different tests.

The permutation (block D) now runs for all seven metrics, with `NPERM = 10000` shuffles
within arm (and within secretion state when that covariate is on), and reports, per gene
set (candidate genes, the six modules, the DEG union):

- `p_perm_metric`: chance of a minimum p this small among that metric's genes;
- `p_perm_family`: the same across genes x all seven metrics ("did the search find anything");
- `p_perm_family_primary`: the same for the pre-declared primary family (candidate genes x
  `PRIMARY_METRICS` from script 20: force, area and their change metrics);
- per gene row: `p_perm_metricwise`, `p_perm_famwise`, `p_perm_primary`, single-step min-p
  adjusted p-values (Westfall & Young 1993) so that no q is ever read without its
  search-corrected companion.

It is fast because the weighted regression is solved by Frisch-Waugh-Lovell residualisation:
weight and residualise the expression matrix on the nuisance design with one matrix
multiply per shuffle, then `t = r sqrt(df) / sqrt(1 - r^2)`. The engine is checked against
`lm()` on the observed data every run (agreement to ~1e-14); 10,000 shuffles over 73 + 6 +
523 genes x 7 metrics take under a minute (foot), about 2.5 minutes for the 1,009-gene gill
DEG union. The null minimum-p vectors are written to `permutation_null_<T>.csv.gz`.

`byssal_structural` is a sixth module (foot proteins, preCols, byssal EP/ACDC, the
plaque-curing tyrosinase: the genes that switch on together when an animal is secreting
thread), separate from the broad `byssal_collagen` regex, so the test "these genes track
thread building, not strength" has its own row. Script 20 adds those genes to the candidate
set from the genome-wide BLAST annotation (`BYSSAL_STRUCTURAL_REGEX`; the DEG-derived map
never reached them, because they are not DEGs in any treatment-vs-control contrast).

### Script 21 writes its diagnostics (block E)

Both diagnostics used to print to the console only. They are now files in
`03_analyses/gene_mechanics/`:

- `detection_floor_flags_<T>.csv`: one row per tested gene (candidate set plus the DEG union
  read back from `assoc_DEGunion_<T>.csv`), with `n_at_floor`, `frac_samples_at_floor` and
  `floor_flag` = `exclude` (> 0.40, `FLOOR_EXCLUDE`), `caution` (> 0.20, `FLOOR_CAUTION`) or
  `ok`. Drop `exclude` genes before interpreting any hit.
- `influence_top_hits_<T>.csv`: the three best candidate hits per metric (`N_INFLUENCE_HITS`;
  ranking by `p_wls` for every metric since 17 September 2026), each refitted as the
  weighted per-animal regression with and without its most influential animal, with `q_wls`
  and the permutation p-values beside them. Columns:
  `most_influential_mussel`, `max_cooks_D`, `cook_threshold_4n`, `cooks_flag` (`D>1`,
  `D>4/n`, `ok`), `slope`, `slope_without_mussel`, `p_without_mussel`, `slope_change_frac`,
  and `influence_flag` = `fragile` when the hit loses p < 0.05 without that animal or its
  slope moves by more than half, else `robust`. The 4/n screen fires for the maximum of ~44
  Cook's distances in almost every fit, so `influence_flag` is the column to read.

- `best_hits_<T>.csv`: the single best hit per metric for each gene set (candidate,
  module, DEG union) with `q_wls`, `p_mixed`, every permutation p, the floor flag and the
  influence flag on one row. This is the table to quote from.

On the current data (17 September 2026) the foot Peroxiredoxin-5 extension hit (T040) is
`fragile`; every other top-three candidate hit is `robust`. Nine foot candidates are at the
detection floor (`exclude`), all of them byssal structural genes that are silent in the
animals that were not secreting thread.

### Script 22 reconciles against script 1's extraction, not raw folders

`raw_post` / `raw_pre` counts come from `thread-summary-raw-output.xlsx`. New flag
`no_post_trace`: sequenced at day 3 but never pulled (currently T136, T137).

### Bioconductor masking

`S4Vectors` and `IRanges`, loaded by DESeq2, mask `dplyr::rename`, `count`, `first` and
`desc`. Script 20 now loads DESeq2 first and tidyverse last, and uses `dplyr::` prefixes in
the new code. Loading tidyverse before DESeq2 reproduces the original `object 'treatment'
not found` failure.

---

## 2. Inputs from thread-strength

| file | produced by | used by |
|---|---|---|
| `03_analyses/thread-summary.xlsx` | scripts 1–2, hand-curated | 20, 22, 23 |
| `03_analyses/decompose-adhesion/mussel_response_classification.csv` | script 4 | 20 (joined into the paired manifest) |
| `03_analyses/extract-tensometer-data/thread-summary-raw-output.xlsx` | script 1 | 22 |

The response classification carries, per animal and per metric, the pre and post means, the
log-ratio, the % change, the raw direction (`decreased` / `increased`), the change relative
to the control arm's mean change, and a composite `response_class` (`weaker` if both force
and adhesion fell, `stronger` if both rose, else `mixed`) and `response_score` (mean
standardised log-ratio across force, area and adhesion). Script 20 joins the class, score,
directions and control-referenced changes into `paired_sample_manifest.csv`.

---

## 3. Outputs

### `03_analyses/gene_mechanics/` (scripts 20 and 21)

| file | contents |
|---|---|
| `paired_sample_manifest_<T>.csv` | the paired animals: arm, day-3 and baseline means, `dlog_*`, response class and score, `secretion_score` / `secretion_state` (joined from `differential-expression/02_data/secretion_state.csv` when it exists) |
| `metrics_config_<T>.csv` | metric list, type, `tier` (primary / exploratory), arm levels, covariates, the `USE_SECRETION_STATE` and `ADJUST_FOR_BASELINE` settings |
| `vst_paired_<T>.csv`, `annotation_map.csv`, `candidate_genes_<T>.csv`, `thread_plaques_paired_<T>.csv` | handoffs |
| `assoc_candidate_<T>.csv`, `assoc_DEGunion_<T>.csv` | script 20 per-gene lm, all seven metrics, `tier` column |
| `assoc_candidate_MIXED_<T>.csv` | script 21 thread-level mixed model, level metrics (sensitivity) |
| `assoc_candidate_WLS_<T>.csv` | script 21 weighted per-animal regression, all metrics; **the reported test**, with `p_mixed` beside it and `p_perm_metricwise` / `p_perm_famwise` / `p_perm_primary` per row |
| `assoc_DEGunion_WLS_<T>.csv` | the same for the DEG union |
| `module_associations_<T>.csv`, `module_members_<T>.csv` | six pathway modules x seven metrics (WLS reported, mixed as sensitivity, permutation columns); the member genes |
| `permutation_summary_<T>.csv` | per gene set: per-metric, family-wide and primary-family permutation p; `permutation_best_hit_<T>.csv` keeps the old per-metric candidate layout |
| `permutation_null_<T>.csv.gz` | the null minimum-p vectors (one column per gene set x metric, plus family minima), git-ignored size aside, regenerated by every run with `PERM_SEED` |
| `best_hits_<T>.csv` | best hit per metric and gene set with q, every permutation p, floor and influence flags |
| `detection_floor_flags_<T>.csv` | every tested gene: fraction of paired samples at the VST floor, `floor_flag` |
| `influence_top_hits_<T>.csv` | top three candidate hits per metric with leave-one-out slope and p, `influence_flag`, q and permutation p |
| `RUN_provenance_<T>.txt` | settings of scripts 20 and 21 (arms, covariates, primary metrics, NPERM, seed), the git commit of the thread input |
| `assoc_candidate_BASELINEADJ_<T>.csv` | ANCOVA mixed model, level metrics |
| `candidate_heatmap_<T>.png`, `top_candidate_scatter_<T>.png`, `best_hit_per_metric_scatter_<T>.png` | figures; the last is the best candidate gene per metric, coloured by arm |

### `03_analyses/expr_tables/` (script 22) and `03_analyses/byssus_genes/` (script 23)

`rna_thread_manifest_<T>.csv` is now tissue-suffixed. The companion `sample_metadata_<T>.csv`
files carry all four arms and `max_displacement`.

---

## 4. Results on the current data (17 September 2026)

Thread input: `thread-summary.xlsx` at commit `9454963` (the corrected traces, results of
record). Outputs regenerated with the 17 September scripts; `USE_SECRETION_STATE = FALSE`.

### Foot

45 paired animals (control 11, OA 12, OW 12, DO 10), 44 with baselines; 73 candidate genes
(58 stress/collagen genes from the TC DEG annotation plus 15 byssal structural genes from
the genome annotation, 17 flagged `byssal_structural`); 523 DEG-union genes; 10,090 genes
after the expression filter.

**Nothing survives the search.** Candidate family: observed minimum p = 0.0012 (PDE8B vs
the change in force, q = 0.088), permutation p = **0.38** across all metrics and **0.26**
for the primary family; the best level-metric hits (HIF-1-alpha vs force, p = 0.011;
Byssal peroxidase-like 2 vs area, p = 0.014) have per-metric permutation p of 0.46 and
0.50. Modules: `HIF_hypoxia` vs the change in area p = 0.013 (q = 0.077) is the best row;
its per-metric permutation p among the six modules is 0.019, but family-wide 0.34. DEG
union: family-wide 0.70. The `byssal_structural` module (16 genes, PC1 63 % of variance)
is unrelated to every metric (best p = 0.18, change in force; famwise 0.99), which is the
direct test that the plaque-protein genes mark thread secretion, not strength.

### Gill

46 paired animals (control 11, OA 12, OW 12, DO 11), 45 with baselines; 66 candidates;
1,009 DEG-union genes. Candidate family: minimum p = 0.0039, permutation p = 0.77 (all)
and 0.59 (primary). Modules: nothing (best p = 0.084). DEG union: Arp2/3 complex subunit
(LOC134706590) vs the change in force reaches q = 0.038 within the metric, but the search
gives per-metric permutation p = 0.077 and family-wide 0.25: a lead for the exploratory
list, not a result.

### What the permutation added

Every q < 0.10 in either tissue has a search-corrected p above 0.05. The earlier
collagen-V lead (q = 0.034 on the pre-correction thread table) does not exist on the
results of record.

## 5. Known data quirks

- **T047** has a foot RNA column and day-3 threads but no row in `F_treatmentinfo.csv`, so it
  is excluded from the paired set. Pre-existing; not introduced here.
- **T136, T137** are sequenced day-3 control animals with no day-3 trace on disk.
- The gill DEG union rests on `GDO_TC_siggene_apeglm.csv`; its companion full table was a
  stale pre-QC fit until 16 September 2026 (rewritten; the 307-gene list was already
  correct, `differential-expression/01_code/DEG-lists-provenance_DOC.md`).
- The candidate set is drawn from the TC-DEG annotation map (plus the byssal structural
  genes), so a stress gene that is not a DEG in any treatment-vs-control contrast is not a
  candidate. That is a hidden selection step; widening the map to the genome annotation
  would raise the foot candidate set from 73 to about 230 genes and dilute the family.

---

## 6. Configuration summary

| script | option | default | effect |
|---|---|---|---|
| 20 | `INCLUDE_CONTROL_ARM` | TRUE | day-3 control animals as a fourth arm |
| 20 | `ADJUST_FOR_BASELINE` | TRUE | ANCOVA for level metrics; n = 44 |
| 20 | `USE_SECRETION_STATE` | FALSE | byssal secretion state (01_7) as a nuisance covariate in 20 and 21; off = treatment only |
| 20 | `PRIMARY_METRICS` | force, area, dlog force, dlog area | the confirmatory family; everything else `exploratory` |
| 20 | `LEVEL_METRICS`, `CHANGE_METRICS` | see above | reporting order = priority |
| 20 | `BYSSAL_STRUCTURAL_REGEX` | foot protein, preCol, ACDC, ... | byssal genes added to the candidate set from the genome annotation |
| 21 | `NPERM`, `PERM_SEED` | 10000, 1 | shuffles for every metric and gene set |
| 21 | `modules` | six regexes | includes `byssal_structural` |
| 22 | `USE_RAW_THREAD_SET` | TRUE | any extracted trace vs curated only |
| 23 | `USE_RAW_THREAD_SET` | FALSE | |
| 20-23 | `params$tissue` | "F" | foot or gill |
| 24 | `params$tissues`, `params$scripts` | both, all four | what the driver renders |
