#!/usr/bin/env bash
# Phase 2: respirometry, instrument-reference, templates, morphometrics,
#          _to_delete, README renames (GitHub), migration-artifact relocation.
# Run from repo ROOT on a branch. `git mv` preserves history; nothing analytical is deleted.
set -euo pipefail
[ -d .git ] || { echo "Run from the repo root."; exit 1; }

# preflight: refuse if any target already exists (avoids half-redo)
for t in respirometry/01_code instrument-reference templates _to_delete _migration morphometrics/02_data; do
  [ -e "$t" ] && { echo "ERROR: '$t' already exists. Reset first: git reset --hard && git clean -fd"; exit 1; }
done
for s in respirometry/metabolic_rate_analysis_baseline.R tagseq/oyster-code \
         morphometrics/morphometrics.xlsx thread-strength/02_data/instrument-reference; do
  [ -e "$s" ] || { echo "ERROR: expected source '$s' missing. Current clone of main?"; exit 1; }
done

# -------------------------------------------------------------- respirometry --
R=respirometry
mkdir -p "$R/01_code" "$R/02_data" "$R/03_analyses/metabolic-rate"
git mv "$R/metabolic_rate_analysis_baseline.R" "$R/01_code/metabolic_rate_analysis_baseline.R"
git mv "$R/metabolic_rate_analysis_T3.R"       "$R/01_code/metabolic_rate_analysis_T3.R"
git mv "$R/measurements" "$R/02_data/measurements"
git mv "$R/constraints"  "$R/02_data/constraints"
git mv "$R/QC_plots"     "$R/03_analyses/metabolic-rate/QC_plots"
git mv "$R/output"       "$R/03_analyses/metabolic-rate/output"
cat > "$R/respirometry.Rproj" <<'RPROJ'
Version: 1.0

RestoreWorkspace: No
SaveWorkspace: No
AlwaysSaveHistory: No

EnableCodeIndexing: Yes
UseSpacesForTab: Yes
NumSpacesForTab: 2
Encoding: UTF-8
RPROJ

# ------------------------------------------------------- instrument-reference --
mkdir -p instrument-reference/thread-tensometer
git mv thread-strength/02_data/instrument-reference/OMEGA_DFG51-2_manual.pdf            instrument-reference/thread-tensometer/OMEGA_DFG51-2_manual.pdf
git mv thread-strength/02_data/instrument-reference/Omega_DFG51-2_Force_Gauge_Logger.vi instrument-reference/thread-tensometer/Omega_DFG51-2_Force_Gauge_Logger.vi
git mv thread-strength/02_data/instrument-reference/wire_colors.txt                     instrument-reference/thread-tensometer/wire_colors.txt
git mv thread-strength/02_data/instrument-reference/original_thread_strength_README.md  instrument-reference/thread-tensometer/original_thread_strength_README.md
rmdir thread-strength/02_data/instrument-reference 2>/dev/null || true

# -------------------------------------------------------------------- templates --
git mv tagseq/oyster-code templates

# ---------------------------------------------------------------- morphometrics --
mkdir -p morphometrics/02_data
git mv morphometrics/morphometrics.xlsx morphometrics/02_data/morphometrics.xlsx

# ------------------------------------------------------ README renames (GitHub) --
git mv thread-strength/thread-strength_README.md                        thread-strength/README.md
git mv thread-strength/02_data/thread-strength_02_data_README.md         thread-strength/02_data/README.md
git mv thread-strength/03_analyses/thread-strength_03_analyses_README.md thread-strength/03_analyses/README.md

# ------------------------------------------------ migration artifacts -> _migration --
mkdir -p _migration
git mv migrate_thread-strength.sh         _migration/migrate_thread-strength.sh
git mv thread-strength_migration_notes.md _migration/thread-strength_migration_notes.md

# ----------------------------------------------------------------- _to_delete --
mkdir -p _to_delete
git mv .Rhistory       _to_delete/root.Rhistory.bak
git mv code/.Rhistory  _to_delete/code.Rhistory.bak
git mv plots/.Rhistory _to_delete/plots.Rhistory.bak
git mv .Rproj.user     _to_delete/Rproj.user.bak

# -------------------------------------------- delete Excel lock temp files --
git rm -q 'tagseq/~$PSMFC-mytilus-byssus-pilot_samplelist_updated-mytilus.xlsx' \
          'tagseq/~$PSMFC-mytilus-byssus-pilot_samplelist_updated.xlsx' 2>/dev/null || true

# ------------------------------------------ keep junk from coming back --
for pat in '.Rhistory' '.Rproj.user/' '~$*'; do
  grep -qxF "$pat" .gitignore 2>/dev/null || echo "$pat" >> .gitignore
done
git add .gitignore

echo "Phase 2 staged."
