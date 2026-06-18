#!/usr/bin/env bash
set -euo pipefail
# Renumber the analyze script so the order is extract(1) -> assemble(2) -> analyze(3).
# Run from the repository root on a branch, review, then push. Re-running is a no-op.
#   git checkout -b feature/thread-summary
#   bash apply_assemble_setup.sh
#   # place 2_assemble_thread_summary.Rmd into thread-strength/01_code/ from your download
#   git add -A && git commit -m "Add thread-summary assembly; renumber analyze to 3" && git push

OLD="thread-strength/01_code/2_analyze_thread_strength.Rmd"
NEW="thread-strength/01_code/3_analyze_thread_strength.Rmd"

if [ -f "$NEW" ]; then
  echo "Already renamed: $NEW exists."
elif [ -f "$OLD" ]; then
  git mv "$OLD" "$NEW"
  echo "Renamed $(basename "$OLD") -> $(basename "$NEW")"
else
  echo "ERROR: neither $OLD nor $NEW found. Run from the repository root." >&2
  exit 1
fi
echo "Now add 2_assemble_thread_summary.Rmd to thread-strength/01_code/ (from your download)."
