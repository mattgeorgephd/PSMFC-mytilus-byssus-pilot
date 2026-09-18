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
20  paired table, VST, candidate set,
    per-gene ANCOVA (the reported test)             ->  03_analyses/gene_mechanics/
21  thread-level mixed ANCOVA, modules, diagnostics ->  03_analyses/gene_mechanics/
22  RNA x thread manifest, top-25 expression tables ->  03_analyses/expr_tables/
23  byssus/foot gene list + expression              ->  03_analyses/byssus_genes/
```

---

## 1. How the chain is set up

### Thread input and labels

Scripts 20, 22 and 23 read the curated thread table `thread-strength/03_analyses/thread-summary.xlsx`
(scripts 1-2, hand-curated). Pre-exposure threads are `phase == "pre"`, day-3 threads
`phase == "post"`, and the arm an animal was assigned to is `mussel_trt`; a
`read_thread_summary()` helper in each script accepts the older column spellings. Script 22
reconciles against script 1's extraction (`thread-summary-raw-output.xlsx`) and flags
animals sequenced at day 3 but never pulled (`no_post_trace`; currently T137).

### The day-3 control arm is in the paired set

`INCLUDE_CONTROL_ARM <- TRUE` in script 20. The day-3 control animals (T126-T137) have
foot RNA, day-3 threads and their own baselines, and enter as a fourth treatment level:
n = 45 foot / 46 gill (44 / 45 with baselines), and the control animals' own trajectory is
the null an expression signal must beat. An association that holds within the control arm
too is about attachment biology, not the stress response. Set FALSE for stressor arms only.

### Metrics, tiers and the reported test

On the thread dataset the stressor effect is in peak force and plaque area separately more
than in their ratio (`thread-strength/01_code/4_decompose_adhesion_DOC.md`: arm x timepoint
p = 2.9e-5 for force, < 1e-10 for area, 0.0057 for adhesion). Script 20 tests, per gene and
metric, one baseline-adjusted regression (ANCOVA) on the per-animal values:

    level_day3 ~ expression + treatment [+ secretion_state] + level_baseline

| metric | scale | tier | per-animal value |
|---|---|---|---|
| `max_force` | log | primary | geometric mean of the animal's day-3 plaques, N |
| `pad_area` | log | primary | geometric mean, mm² |
| `adhesion_kpa` | log | exploratory | geometric mean of force / area × 1000 (recomputed per plaque) |
| `max_displacement` | raw | exploratory | arithmetic mean of extension at break, mm |

Force, area and adhesion enter as log(geometric mean) on both sides of the model, so
`slope` is a log-unit change per VST unit and exp(slope) a multiplicative one; extension is
raw. `level_baseline` is the same summary of the animal's own pre-exposure threads, on the
same scale. Animals without baseline threads are not in the fits (44 of 45 foot, 45 of 46
gill). `baseline_slope` is reported beside `slope`; on the current data it is 0.28-0.37 for
force and adhesion and 0.03-0.14 for area (median over genes), consistent with the low
repeatability in `thread-strength/01_code/3_analyze_thread_strength_DOC.md` (ICC 0.13-0.38).

`METRICS` in script 20 fixes the scale and the tier; `tier = primary` (force, area) is the
declared confirmatory family and every output carries the column. Multiplicity: `q_lm` is BH
within a gene set x metric, `q_family` BH within a gene set x tier, so the candidate x
primary family (201 x 2 foot, 251 x 2 gill tests) has its own search-corrected q. Beside
the parametric slope, `rho_partial` is the Spearman correlation with arm, state and
baseline partialled out of the ranks, on the same animals. `metrics_config_<T>.csv` is
written by 20 and read by 21, so the metric list, scales, tiers, arm levels and covariates
cannot drift between the two.

Why an ANCOVA and not a change score: the change score `log(day3) - log(baseline)` imposes a
baseline coefficient of 1, while the fitted coefficient is 0.03-0.37, so the change score
adds most of the baseline's measurement noise to the response and loses power; the ANCOVA
also absorbs the between-animal baseline differences the arm assignment did not balance
(baseline area and extension differ by future arm, script 3, 3b). On the current data the
per-animal fit and the thread-level mixed fit agree to within a few percent in p because 44
of 45 foot and 45 of 46 gill day-3 animals have exactly three plaques.

### Annotation map and candidate universe

`CANDIDATE_ANNOTATION = "genome"`: every gene in the count matrix is annotated with its best
UniProt hit (highest bitscore) from `blast/03_analyses/genome-foot/LOC_GO_list.txt`, with
`blast_pident` and `blast_evalue` carried along, and any expressed gene whose name matches
`CANDIDATE_KEYWORDS` (byssal / collagen / plaque-curing / HSP / hypoxia / tRNA-synthetase /
oxidative-stress terms) and passes the BLAST floor (`CANDIDATE_MAX_EVALUE = 1e-10`,
`CANDIDATE_MIN_PIDENT = 0`) is a candidate: 226 foot, 279 gill. The alternative
`"TC_DEG"` restricts candidates to genes named in the treatment-vs-control DEG tables (57 /
62 of those) plus the byssal structural genes; it made "already a DEG in some contrast" a
hidden entry condition and is kept only for comparison. `in_TC_DEG_annotation` marks the
overlap in every table.

### Detection floor

A gene at the VST floor (zero counts) in many paired animals gives an association driven by
presence/absence. Script 20 computes `frac_at_floor` for every tested gene, flags `exclude`
(> 0.40) and `caution` (> 0.20), writes `detection_floor_flags_<T>.csv` for the candidate
set and the DEG union, and with `EXCLUDE_FLOOR_GENES = TRUE` (default) drops the `exclude`
genes from both families before testing: 201 of 226 foot candidates and 251 of 279 gill
candidates are tested, 520 of 523 and 998 of 1,009 DEG-union genes. Excluded genes still
contribute to the module scores in script 21 (the `byssal_structural` module is mostly
floor genes by construction).

### Secretion state

`differential-expression/02_data/secretion_state.csv` (DE script `01_7`) labels each foot
sample `on` / `off` for the plaque-protein module; script 20 joins it to the manifest by
animal for both tissues (15 of 45 foot animals `on`). `USE_SECRETION_STATE` (default FALSE)
adds it as a nuisance covariate in every fit of 20 and 21.

### Script 21: sensitivity, modules, diagnostics

Script 21 reads script 20's association tables and does not re-derive the reported test. It
adds:

- **Block A, thread-level mixed ANCOVA** (`assoc_candidate_MIXED_<T>.csv`, `p_mixed`): the
  same model on every day-3 plaque, `plaque ~ expression + treatment [+ secretion_state] +
  level_baseline + (1 | mussel)`, plaque values on the metric's scale, Satterthwaite df
  (lmerTest). It weights each animal by the precision of its mean instead of equally; on the
  current data (three plaques for all but one animal per tissue) it reproduces `p_lm` to
  within a few percent, which is the check that no hit is carried by plaque count.
- **Block B, module eigengenes** (`module_associations_<T>.csv`, `module_members_<T>.csv`):
  PC1 of each module's members (genes passing the BLAST floor, `blast_ok`), oriented so a
  higher score is higher expression, through the same ANCOVA (`p_lm`, `q_lm` across the six
  modules within a metric, `q_family` within a tier) and the mixed ANCOVA (`p_mixed`).
  `byssal_structural` is a sixth module (foot proteins, preCols, byssal EP/ACDC, the
  plaque-curing tyrosinase: the genes that switch on together when an animal is secreting
  thread), separate from the broad `byssal_collagen` regex, so the test "these genes track
  thread building, not strength" has its own row (`BYSSAL_STRUCTURAL_REGEX` in script 20
  flags the same genes in the candidate table).
- **Block C, diagnostics**, below.

### Script 21 diagnostics (block C)

- `detection_floor_flags_<T>.csv` (written by script 20, echoed here): one row per gene in
  the candidate set and the DEG union, with `frac_at_floor`, `floor_flag` and `tested`.
- `influence_top_hits_<T>.csv`: the three best candidate hits and the best DEG-union hit per
  metric (`N_INFLUENCE_HITS`, ranked by `p_lm`), each refitted as the reported ANCOVA with
  and without its most influential animal. Columns: `most_influential_mussel`,
  `max_cooks_D`, `cook_threshold_4n`, `cooks_flag` (`D>1`, `D>4/n`, `ok`), `slope`,
  `slope_without_mussel`, `p_without_mussel`, `slope_change_frac`, and `influence_flag` =
  `fragile` when the hit loses p < 0.05 without that animal or its slope moves by more than
  half, else `robust`. The 4/n screen fires for the maximum of ~45 Cook's distances in
  almost every fit, so `influence_flag` is the column to read.
- `best_hits_<T>.csv`: the single best hit per metric for each gene set (candidate, module,
  DEG union) with `p_lm`, `q_lm`, `q_family`, `p_mixed`, the floor flag and the influence
  flag on one row. This is the table to quote from.

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
standardised log-ratio across force, area and adhesion). Script 20 joins the class and the
score into `paired_sample_manifest_<T>.csv` for inspection; they enter no model.

---

## 3. Outputs

### `03_analyses/gene_mechanics/` (scripts 20 and 21)

| file | contents |
|---|---|
| `paired_sample_manifest_<T>.csv` | the paired animals: arm, plaque counts, day-3 and baseline per-animal values (geometric means for force, area, adhesion), response class and score, `secretion_score` / `secretion_state` (joined from `differential-expression/02_data/secretion_state.csv` when it exists) |
| `metrics_config_<T>.csv` | metric list with `scale` (log / raw), `tier` (primary / exploratory) and label, arm levels, covariates, the `USE_SECRETION_STATE` and `EXCLUDE_FLOOR_GENES` settings |
| `vst_paired_<T>.csv`, `thread_plaques_paired_<T>.csv` | handoffs |
| `annotation_map.csv` | genome-wide best UniProt hit per LOC with `blast_pident`, `blast_evalue`, `blast_ok`, `in_TC_DEG_annotation` |
| `candidate_genes_<T>.csv` | the candidate set with `byssal_structural`, `frac_at_floor`, `floor_flag`, `tested` |
| `detection_floor_flags_<T>.csv` | every candidate and DEG-union gene: fraction of paired samples at the VST floor, `floor_flag`, `tested` |
| `assoc_candidate_<T>.csv`, `assoc_DEGunion_<T>.csv` | **the reported test**: script 20 ANCOVA per gene x metric with `scale`, `tier`, `n`, `slope`, `se`, `p_lm`, `baseline_slope`, `rho_partial`, `q_lm`, `q_family`, floor flag |
| `assoc_candidate_MIXED_<T>.csv` | script 21 thread-level mixed ANCOVA, all metrics (sensitivity, `p_mixed`) |
| `module_associations_<T>.csv`, `module_members_<T>.csv` | six pathway modules x four metrics (ANCOVA reported, mixed as sensitivity); the member genes |
| `best_hits_<T>.csv` | best hit per metric and gene set with `p_lm`, `q_lm`, `q_family`, `p_mixed`, floor and influence flags |
| `influence_top_hits_<T>.csv` | top three candidate hits and the best DEG-union hit per metric with leave-one-out slope and p, `influence_flag`, `q_lm`, `q_family` |
| `RUN_provenance_<T>.txt` | settings of scripts 20 and 21 (arms, covariates, model, metrics with scale and tier, modules), the git commit of the thread input |
| `candidate_heatmap_<T>.png`, `top_candidate_scatter_<T>.png`, `best_hit_per_metric_scatter_<T>.png` | figures; the scatter y axes are baseline-adjusted day-3 levels on the model scale; the last is the best candidate gene per metric, coloured by arm |

### `03_analyses/expr_tables/` (script 22) and `03_analyses/byssus_genes/` (script 23)

`rna_thread_manifest_<T>.csv` is now tissue-suffixed. The companion `sample_metadata_<T>.csv`
files carry all four arms and `max_displacement`.

---

## 4. Results on the current data (18 September 2026, ANCOVA)

Thread input: `thread-summary.xlsx` at commit `9454963`. Genome-wide candidate map, floor
genes excluded, `USE_SECRETION_STATE = FALSE`, the ANCOVA above as the reported test. The
committed tables are the run of 18 September 2026 on the analysis machine (R 4.2.2, DESeq2
1.38.3, lme4 1.1-31); an independent run under R 4.3.3 / DESeq2 1.42.0 / lme4 1.1-35
reproduced every ANCOVA p-value to a relative 1e-13 and every mixed-model p-value to 2e-6,
so the figures below hold to the two significant figures quoted.

### Foot

45 paired animals (control 11, OA 12, OW 12, DO 10), 44 in the fits (OW 11); 201 candidate
genes tested (226 matched, 25 at the detection floor); 520 DEG-union genes tested; 10,090
genes after the expression filter.

**Nothing.** Candidate x primary family (402 tests): minimum p = 0.0064 (PDE8B vs force,
slope -0.66 log N per VST unit, `q_lm` 0.76, `q_family` 0.91); 8 of 201 force tests and 5
of 201 area tests have p < 0.05 against 10 expected by chance. Exploratory: Collagen
alpha-1(XXII) vs extension p = 0.00083 (`q_lm` 0.17, `q_family` 0.33, `caution` floor flag)
is the smallest p in the tissue; adhesion best p = 0.015 (GST A4, q = 0.82). Modules: best
row `tRNA_translation` (29 genes) vs adhesion p = 0.053 (q = 0.32), vs area p = 0.074; the
`byssal_structural` module (16 genes, PC1 63 % of variance) is unrelated to every metric
(best p = 0.26), the direct test that the plaque-protein genes mark thread secretion, not
strength. DEG union: best Zinc finger CCCH 18 vs force p = 0.00061 (`q_lm` 0.32); no q
below 0.10 anywhere; 19 of 520 force tests and 7 of 520 area tests at p < 0.05 against 26
expected.

### Gill

46 paired animals (control 11, OA 12, OW 12, DO 11), 45 in the fits; 251 candidates tested
(279 matched, 28 at the floor); 998 DEG-union genes tested; 13,233 genes after the filter.

The gill carries the only associations below q = 0.10, all with positive slopes (higher
day-3 expression, stronger attachment, net of arm and baseline) and all `robust` to the
most influential animal:

- Heat shock 70 kDa protein 12A, LOC134718614, vs adhesion: p = 3.2e-5, slope +0.97 log kPa
  per VST unit, `q_lm` 0.0079, `q_family` 0.016 (exploratory tier); vs force: p = 2.5e-4,
  slope +0.78, `q_lm` 0.062, `q_family` 0.12 (primary family). Its paralog LOC134718612 vs
  area: p = 4.6e-4, slope -0.37, `q_lm` 0.12 (`caution` floor flag); LOC134697060 vs force
  and adhesion p = 0.009 / 0.002. Three HSPA12A paralogs at the top of three metrics is
  the one pattern worth a follow-up; it is a gill (systemic-state) readout, not foot.
- DEG union: Arp2/3 complex subunit 2, LOC134706590, vs force p = 3.2e-5, slope +1.29,
  `q_lm` 0.032, `q_family` 0.064; vs adhesion p = 2.2e-4 (`q_lm` 0.11). Lactadherin,
  LOC134683070, vs adhesion p = 2.1e-4 (`q_lm` 0.11) and force p = 4.0e-4 (`q_lm` 0.20).
  100 of 998 force tests and 128 of 998 extension tests have p < 0.05 against 50 expected:
  a broad, correlated gill expression signature tracks force, which is what a systemic
  covariate (condition, handling response) looks like, not a gene-specific effect.
- Modules: nothing (best p = 0.11, `tRNA_translation` vs extension).

Reading these: `q_family` = 0.016 for HSPA12A vs adhesion is a false-discovery rate within
the candidate x exploratory family (502 correlated tests), not a family-wise error rate; in
the primary family the same gene sits at `q_family` = 0.12. These are leads for the
exploratory list and a targeted follow-up, not results. The candidate x primary family has
no q below 0.10 in either tissue.

## 5. Known data quirks

- **T047** has a foot RNA column and day-3 threads but no row in `F_treatmentinfo.csv`, so it
  is excluded from the paired set. Pre-existing; not introduced here.
- **T137** is a sequenced day-3 control animal with no day-3 trace on disk.
- The gill DEG union rests on `GDO_TC_siggene_apeglm.csv`; its companion full table was a
  stale pre-QC fit until 16 September 2026 (rewritten; the 307-gene list was already
  correct, `differential-expression/01_code/DEG-lists-provenance_DOC.md`).
- The candidate keywords are regexes on UniProt names; `Hsp` and `chaperone` in particular
  pull in co-chaperones and assembly factors, so the `HSP_proteostasis` module is broad (92
  foot / 140 gill genes). Tighten `CANDIDATE_KEYWORDS` or raise `CANDIDATE_MIN_PIDENT`
  (40 halves the candidate set) if a narrower family is wanted.

---

## 6. Configuration summary

| script | option | default | effect |
|---|---|---|---|
| 20 | `INCLUDE_CONTROL_ARM` | TRUE | day-3 control animals as a fourth arm |
| 20 | `USE_SECRETION_STATE` | FALSE | byssal secretion state (DE 01_7) as a nuisance covariate in 20 and 21; off = treatment only |
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
