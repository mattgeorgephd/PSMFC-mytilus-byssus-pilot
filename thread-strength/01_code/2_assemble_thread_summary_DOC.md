# 2_assemble_thread_summary.Rmd

Turns the trace-level extraction into a **curation-ready candidate** shaped like
`03_analyses/thread-summary.xlsx`, so the manual step is a review rather than a retype.

## Pipeline position

```
1_extract_tensometer_data.Rmd  ->  03_analyses/thread-summary-raw-output.xlsx   (375 traces)
2_assemble_thread_summary.Rmd  ->  03_analyses/assemble-thread-summary/
                                     thread-summary-candidate.xlsx
   (manual: add pad_area / failure, drop bad runs)
                               ->  03_analyses/thread-summary.xlsx
3_analyze_thread_strength.Rmd  ->  03_analyses/analyze-thread-strength/
```

## What changed

The previous version read four extraction files (`max_force_{control,treatment}.xlsx`,
`integral_{control,treatment}.xlsx`), `morphometrics/02_data/morphometrics.xlsx`, and
`02_data/GOOGLESHEET-PSMFC-mytilus-byssus-pilot-threads.xlsx`. **None of those paths exist.**

Now: one extraction file, with mussel-level metadata already joined by script 1 from the
mussel key. `group` is gone, replaced by `phase`.

## Inputs

| Path | Required | Role |
|---|---|---|
| `03_analyses/thread-summary-raw-output.xlsx` | yes | trace-level extraction |
| `03_analyses/thread-summary.xlsx` | no | source of `pad_area` and `failure` |

## The join key

`mussel` + `thread` + `thread_trt`. Verified unique on both sides, and `thread_trt` is not
optional: 45 animals were pulled before and after exposure, so `mussel` + `thread` alone
matches two different traces and would silently duplicate rows.

## Output

`03_analyses/assemble-thread-summary/thread-summary-candidate.xlsx`, sheet `data`.

Columns match `thread-summary.xlsx` except:

- `group` is replaced by `phase` (`lab` / `pre` / `post`)
- adds `note`, `max_displacement`, and the `source_folder` / `file` provenance pair

Drop `source_folder` and `file` before pasting if you want an exact column match.

## Checks worth reading before curating

1. **Peak-force disagreement.** Any row where the freshly extracted `max_force` differs from
   the value already in `thread-summary.xlsx` is tabulated. This is not rounding: it means
   the curated value came from a different file or was hand-edited, and it fed straight into
   `adhesion_kpa`. Three such rows exist today, listed in `1_extract_tensometer_data_DOC.md`.
2. **Orphan curated rows.** Rows in `thread-summary.xlsx` with no matching trace. Currently
   none.
3. **What still needs measuring.** Per folder, how many traces lack `pad_area`.
4. **Pairing by arm.** How many animals in each arm have both a `pre` and a `post` pull.
   Read this before fitting a mixed model: unbalanced pairing is what makes a per-arm
   marginal mean from a random-intercept fit unreliable, and it is the known cause of the
   OA/OW rank swap between the `lmer` contrasts and the paired t-tests.

Nothing is excluded here. Every extracted trace is carried through so bad runs can be
reviewed against their QC plots and removed deliberately.
