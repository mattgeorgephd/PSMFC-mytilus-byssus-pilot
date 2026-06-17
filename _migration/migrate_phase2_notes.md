# Phase 2 migration notes

Covers respirometry, instrument-reference, templates, morphometrics, `_to_delete`,
the GitHub README renames, and migration-artifact relocation. All moves use `git mv`.

## Apply

1. `git switch -c reorg/phase2`
2. `bash _migration/migrate_phase2.sh` from the repo root (place the script first, or run
   from wherever you keep it).
3. Copy the rewritten respirometry scripts into `respirometry/01_code/`, overwriting the
   moved originals.
4. Copy the README files into place (see the path mapping in chat). The script already
   renames the thread-strength READMEs to `README.md`; overwrite those two with the
   updated versions, which drop the instrument-reference references.
5. **Also copy the rewritten `1_extract_tensometer_data.Rmd` into `thread-strength/01_code/`.**
   The version currently on main is the stale pre-`here()` original and will not run.
6. `git status && git diff`, knit/run to confirm, commit, push.

## Respirometry rewrite

- `setwd(dirname(current_path))` (rstudioapi) and every `setwd('..')` chain replaced with
  `here::here()`. The two scripts had different numbers of `..` for the same folder depth,
  a sign the originals were hand-tuned to a starting directory and were fragile.
- Reads use `list.files(meas_dir)` + `read_excel(file.path(meas_dir, i))`, keeping the bare
  filename `i` for the downstream name parsing.
- **Constraints are no longer clobbered.** The template write is guarded by
  `if (!file.exists(constr_file))`, so a re-run does not wipe the hand-edited start/stop
  windows. This is a behavior change from the original; revert if undesired.
- Baseline still writes `SMR_slope_heat_only_T5.xlsx` (a leftover name); flagged in a code
  comment and the folder README.

## What was verified / not

- Verified: phase-2 `git mv` runs clean against current `main`; respirometry scripts have
  no `setwd`/`getwd`/`rstudioapi` left and every `here()` target resolves on disk.
- Not verified: a full run (no R/CRAN in the sandbox). Run both scripts to confirm.
