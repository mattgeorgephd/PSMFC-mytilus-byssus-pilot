#!/usr/bin/env bash
# =============================================================================
# Migration: thread-strength  ->  standardized 01_code / 02_data / 03_analyses
# Moves scattered thread-strength files into one self-contained analysis folder
# using `git mv` (history preserved). Nothing is deleted.
#
# RUN from the repository ROOT, ideally on a branch:
#     git switch -c reorg/thread-strength
#     bash migrate_thread-strength.sh
# Then copy the rewritten .Rmd files into thread-strength/01_code/, add the
# README files, review `git status` / `git diff`, knit, commit, push.
# =============================================================================
set -euo pipefail

# --- preflight ---------------------------------------------------------------
if [ ! -d .git ]; then
  echo "ERROR: run this from the repository root (no .git directory here)."; exit 1
fi
if [ -e "thread-strength" ]; then
  echo "ERROR: 'thread-strength/' already exists. A previous run probably failed partway."
  echo "       Reset to a clean state first, then re-run this script:"
  echo "         git reset --hard && git clean -fd thread-strength"
  exit 1
fi
for src in code/1_extract_tensometer_data.Rmd \
           thread_strength/tensometer_output \
           thread_strength/QC_plots \
           analyses/thread_strength; do
  [ -e "$src" ] || { echo "ERROR: expected '$src' not found. Is this a pristine clone at the repo root?"; exit 1; }
done

A="thread-strength"

mkdir -p "$A/01_code" \
         "$A/02_data/instrument-reference" \
         "$A/03_analyses/extract-tensometer-data/summarized_data" \
         "$A/03_analyses"

# --- 01_code -----------------------------------------------------------------
git mv code/1_extract_tensometer_data.Rmd "$A/01_code/1_extract_tensometer_data.Rmd"
git mv code/2_analyze_thread_strength.Rmd "$A/01_code/2_analyze_thread_strength.Rmd"
git mv code/thread-code                   "$A/01_code/thread-code"

# --- 02_data -----------------------------------------------------------------
git mv thread_strength/tensometer_output            "$A/02_data/tensometer_output"
git mv thread_strength/pictures                     "$A/02_data/pictures"
git mv thread_strength/summarized_data/summary.xlsx "$A/02_data/summary.xlsx"
git mv thread_strength/Omega_DFG51-2_Force_Gauge_Logger.vi "$A/02_data/instrument-reference/Omega_DFG51-2_Force_Gauge_Logger.vi"
git mv thread_strength/OMEGA_DFG51-2_manual.pdf            "$A/02_data/instrument-reference/OMEGA_DFG51-2_manual.pdf"
git mv thread_strength/wire_colors.txt                     "$A/02_data/instrument-reference/wire_colors.txt"
git mv thread_strength/README.md.txt                       "$A/02_data/instrument-reference/original_thread_strength_README.md"

# --- 03_analyses : outputs binned by the script that produced them -----------
# script 1 (extract): QC plots (whole dir) + summarized force/AUC tables
git mv thread_strength/QC_plots "$A/03_analyses/extract-tensometer-data/QC_plots"
for f in max_force_control integral_control max_force_treatment integral_treatment; do
  git mv "thread_strength/summarized_data/$f.xlsx" \
         "$A/03_analyses/extract-tensometer-data/summarized_data/$f.xlsx"
done
# script 2 (analyze): boxplots / line graphs
git mv analyses/thread_strength "$A/03_analyses/analyze-thread-strength"

# --- RStudio project file ----------------------------------------------------
cat > "$A/thread-strength.Rproj" <<'RPROJ'
Version: 1.0

RestoreWorkspace: No
SaveWorkspace: No
AlwaysSaveHistory: No

EnableCodeIndexing: Yes
UseSpacesForTab: Yes
NumSpacesForTab: 2
Encoding: UTF-8

RnwWeave: Sweave
LaTeX: pdfLaTeX
RPROJ

# --- tidy emptied source dirs ------------------------------------------------
rmdir thread_strength/summarized_data thread_strength 2>/dev/null || true

echo
echo "thread-strength migration staged."
echo "NEXT: copy the rewritten .Rmd files into thread-strength/01_code/,"
echo "      add the README.md files, then: git status && git diff"
