#!/usr/bin/env bash
set -euo pipefail
# Point 3_analyze_thread_strength.Rmd at the curated clean file and fix the load.lib identifier
# corruption (the old oa->DO find/replace turned every "load"/"Load" into "lDOd"/"LDOd"). These
# tokens appear only as the corruption, so the replacement is global and complete. Run from the repo
# root on a branch, review, push. Re-running is a no-op.
#   git checkout -b fix/analyze-input
#   bash apply_3analyze_downstream.sh
#   git diff thread-strength/01_code/3_analyze_thread_strength.Rmd

F="thread-strength/01_code/3_analyze_thread_strength.Rmd"
[ -f "$F" ] || { echo "ERROR: $F not found. Run from the repository root." >&2; exit 1; }

echo "Before:  lDOd/LDOd=$(grep -cE 'lDOd|LDOd' "$F" || true)  \"summary.xlsx\"=$(grep -c '"summary\.xlsx"' "$F" || true)"
perl -i -pe '
  s/lDOd/load/g;
  s/LDOd/Load/g;
  s/"summary\.xlsx"/"thread-summary-trossulus-clean.xlsx"/g;
' "$F"
echo "After:   lDOd/LDOd=$(grep -cE 'lDOd|LDOd' "$F" || true) (want 0)  load.lib=$(grep -c 'load\.lib' "$F" || true) (want 3)  clean-file ref=$(grep -c 'thread-summary-trossulus-clean\.xlsx' "$F" || true) (want 1)"
echo "Done. Review: git diff $F"
