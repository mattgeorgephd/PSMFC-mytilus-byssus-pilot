# Gene-mechanics correlation pipeline, scripts 20 to 23

Links foot (or gill) gene expression at day 3 to the same animal's byssal thread mechanics.
Four chained scripts, each reading the previous one's CSV handoffs rather than sharing an R
session. Run in order: **20 → 21 → 22 → 23**, after thread-strength scripts 1 to 4.

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

### Script 21 handles change metrics correctly

A change metric is one number per animal, so the thread-level mixed model (block A) applies
to level metrics only. For change metrics the reported test is the per-animal weighted
regression (block B), with weight = min(baseline plaques, day-3 plaques). The module
analysis (block C) uses the mixed model for level metrics and the weighted regression for
change metrics; the permutation (block D) now runs for `max_force`, `pad_area` and
`dlog_adhesion_kpa` and writes `permutation_best_hit_<T>.csv`.

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
| `paired_sample_manifest.csv` | the 44 animals: arm, day-3 and baseline means, `dlog_*`, response class and score |
| `metrics_config.csv` | metric list, type, arm levels |
| `vst_paired_<T>.csv`, `annotation_map.csv`, `candidate_genes_<T>.csv`, `thread_plaques_paired_<T>.csv` | handoffs |
| `assoc_candidate_<T>.csv`, `assoc_DEGunion_<T>.csv` | script 20 per-gene lm, all seven metrics |
| `assoc_candidate_MIXED_<T>.csv` | script 21 thread-level mixed model, level metrics |
| `assoc_candidate_WLS_<T>.csv` | script 21 weighted per-animal regression, all metrics; **the reported test for change metrics** |
| `module_associations_<T>.csv` | five pathway modules × seven metrics |
| `permutation_best_hit_<T>.csv` | search-corrected p for the best candidate hit, three metrics |
| `assoc_candidate_BASELINEADJ_<T>.csv` | ANCOVA mixed model, level metrics |
| `candidate_heatmap_<T>.png`, `top_candidate_scatter_<T>.png` | figures |

### `03_analyses/expr_tables/` (script 22) and `03_analyses/byssus_genes/` (script 23)

Unchanged in shape; the companion `sample_metadata_<T>.csv` files now carry all four arms
and `max_displacement`.

---

## 4. Results on the current data (foot, 16 September 2026)

Paired animals 44 (control 10, OA 12, OW 12, DO 10); 42 with baselines. 57 candidate genes,
523 DEG-union genes, 9,974 genes after the expression filter.

**Level metrics: nothing.** Best candidate hit on `max_force` p = 0.013, q = 0.70;
permutation p = 0.44. On `pad_area` q = 0.72, permutation p = 0.50. Nothing in the
DEG union passes q < 0.10 on any metric.

**Change metrics: one borderline signal, in the collagen machinery.**

| gene | metric | slope | p | q (WLS) |
|---|---|---|---|---|
| Collagen alpha-1(V) chain | `dlog_adhesion_kpa` | −0.38 | 0.0006 | **0.034** |
| Collagen alpha-2(IV) chain | `dlog_adhesion_kpa` | −0.51 | 0.002 | 0.064 |
| Prolyl 4-hydroxylase alpha-2 | `dlog_adhesion_kpa` | −0.42 | 0.013 | 0.25 |
| Collagen alpha-2(IV) chain | `dlog_pad_area` | +0.30 | 0.008 | 0.30 |

Higher collagen expression goes with a **larger fall in adhesion** and a **larger rise in
plaque area**; the association with the change in force is weak (q = 0.58). The reading is
that animals investing in collagen built bigger plaques that did not hold proportionally
more. Modules agree: `byssal_collagen` vs `dlog_adhesion_kpa` q = 0.052;
`HIF_hypoxia` vs `dlog_pad_area` q = 0.051.

**How much to trust it.** The permutation-corrected p for the best `dlog_adhesion_kpa` hit,
across the 57-gene search, is **0.054**. Without the control arm (n = 32 with baselines)
Collagen V is still the top hit with the same sign and a similar slope (−0.43, p = 0.004)
but q = 0.21. So: stable in direction and rank across specifications, borderline after
correction for the search, and its BH significance depends on the ten control animals.
Hypothesis-strengthening, not a result to state as established.

For comparison, the pre-change analysis (31 animals, 27 with baselines) had a minimum
candidate q of 0.29 and a minimum module q of 0.20.

---

## 5. Known data quirks

- **T047** has a foot RNA column and day-3 threads but no row in `F_treatmentinfo.csv`, so it
  is excluded from the paired set. Pre-existing; not introduced here.
- **T136, T137** are sequenced day-3 control animals with no day-3 trace on disk.
- The gill arm (`TISSUE = "G"`) has not been run.

---

## 6. Configuration summary

| script | option | default | effect |
|---|---|---|---|
| 20 | `INCLUDE_CONTROL_ARM` | TRUE | day-3 control animals as a fourth arm |
| 20 | `ADJUST_FOR_BASELINE` | TRUE | ANCOVA for level metrics; n = 42 |
| 20 | `LEVEL_METRICS`, `CHANGE_METRICS` | see above | reporting order = priority |
| 21 | `PMETRICS` | force, area, dlog adhesion | which metrics get a permutation null |
| 21 | `NPERM` | 1000 | |
| 22 | `USE_RAW_THREAD_SET` | TRUE | any extracted trace vs curated only |
| 23 | `USE_RAW_THREAD_SET` | FALSE | |
