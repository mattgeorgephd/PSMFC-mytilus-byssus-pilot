#!/usr/bin/env bash
# Phase 3a: summary-plots (cross-cutting; reads respirometry + morphometrics).
# Run from repo ROOT on a branch. git mv preserves history.
set -euo pipefail
[ -d .git ] || { echo "Run from the repo root."; exit 1; }
[ -e summary-plots ] && { echo "ERROR: 'summary-plots' already exists. Reset: git reset --hard && git clean -fd"; exit 1; }
[ -e code/plot_analysis_mussel.R ] || { echo "ERROR: code/plot_analysis_mussel.R missing."; exit 1; }
[ -d plots ] || { echo "ERROR: plots/ missing."; exit 1; }

S=summary-plots
mkdir -p "$S/01_code" "$S/03_analyses"
git mv code/plot_analysis_mussel.R "$S/01_code/plot_analysis_mussel.R"
git mv plots "$S/03_analyses/plot-analysis-mussel"
cat > "$S/summary-plots.Rproj" <<'RPROJ'
Version: 1.0

RestoreWorkspace: No
SaveWorkspace: No
AlwaysSaveHistory: No

EnableCodeIndexing: Yes
UseSpacesForTab: Yes
NumSpacesForTab: 2
Encoding: UTF-8
RPROJ
echo "Phase 3a staged."
