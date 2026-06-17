#!/usr/bin/env bash
# =============================================================================
# Migration: thread-strength  ->  standardized 01_code / 02_data / 03_analyses
# -----------------------------------------------------------------------------
# WHAT IT DOES
#   Moves the thread-strength files out of their scattered locations
#   (code/, thread_strength/, analyses/thread_strength/) into a single
#   self-contained analysis folder, using `git mv` so file history is kept.
#   Nothing is deleted.
#
# HOW TO RUN
#   1. From the ROOT of a clean, up-to-date clone:
#        bash migrate_thread-strength.sh
#   2. Overwrite the two moved scripts in thread-strength/01_code/ with the
#      rewritten versions provided alongside this script (paths converted to
#      here::here()). Then drop the provided README.md files into place.
#   3. Review with `git status` and `git diff`, then commit and push.
#
# SAFETY
#   - Run on a branch, not main:  git switch -c reorg/thread-strength
#   - `set -euo pipefail` aborts on the first error so a partial move is obvious.
# =============================================================================
set -euo pipefail

if [ ! -d .git ]; then
  echo "ERROR: run this from the repository root (no .git directory found here)."
  exit 1
fi

A="thread-strength"

mkdir -p "$A/01_code" \
         "$A/02_data/instrument-reference" \
         "$A/03_analyses/extract-tensometer-data/summarized_data" \
         "$A/03_analyses"

# --- 01_code : the two analysis scripts + the thread-code note ---------------
git mv code/1_extract_tensometer_data.Rmd "$A/01_code/1_extract_tensometer_data.Rmd"
git mv code/2_analyze_thread_strength.Rmd "$A/01_code/2_analyze_thread_strength.Rmd"
git mv code/thread-code                   "$A/01_code/thread-code"

# --- 02_data : raw inputs + curated summary + instrument/reference docs ------
git mv thread_strength/tensometer_output            "$A/02_data/tensometer_output"
git mv thread_strength/pictures                     "$A/02_data/pictures"
git mv thread_strength/summarized_data/summary.xlsx "$A/02_data/summary.xlsx"
git mv thread_strength/Omega_DFG51-2_Force_Gauge_Logger.vi "$A/02_data/instrument-reference/Omega_DFG51-2_Force_Gauge_Logger.vi"
git mv thread_strength/OMEGA_DFG51-2_manual.pdf            "$A/02_data/instrument-reference/OMEGA_DFG51-2_manual.pdf"
git mv thread_strength/wire_colors.txt                     "$A/02_data/instrument-reference/wire_colors.txt"
git mv thread_strength/README.md.txt                       "$A/02_data/instrument-reference/original_thread_strength_README.md"

# --- 03_analyses : outputs binned by the script that produced them -----------
# script 1 (extract): QC plots + summarized force/AUC tables
git mv thread_strength/QC_plots/control   "$A/03_analyses/extract-tensometer-data/QC_plots/control"
git mv thread_strength/QC_plots/treatment "$A/03_analyses/extract-tensometer-data/QC_plots/treatment"
for f in max_force_control integral_control max_force_treatment integral_treatment; do
  git mv "thread_strength/summarized_data/$f.xlsx" \
         "$A/03_analyses/extract-tensometer-data/summarized_data/$f.xlsx"
done
# script 2 (analyze): boxplots / line graphs (whole directory)
git mv analyses/thread_strength "$A/03_analyses/analyze-thread-strength"

# --- RStudio project file (sets working dir to the analysis root) ------------
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

# --- tidy up now-empty leftover directories ----------------------------------
rmdir thread_strength/QC_plots thread_strength/summarized_data thread_strength 2>/dev/null || true

echo
echo "thread-strength migration staged."
echo "NEXT: copy the rewritten 1_extract_tensometer_data.Rmd and"
echo "      2_analyze_thread_strength.Rmd into thread-strength/01_code/,"
echo "      add the README.md files, then: git status && git diff"
