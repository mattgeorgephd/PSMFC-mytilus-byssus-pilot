#!/usr/bin/env bash
#
# Rename the tensometer output subfolders to a phase-first scheme.
#
#   00_laboratory_control  ->  00_lab_reference     day 0, never entered the experimental system
#   00_baseline            ->  01_pre_exposure      day 1, shared holding, before exposure
#   01_treatment_control   ->  02_post_control      day 3, common-garden control tank
#   02_OA_treatment        ->  03_post_OA           day 3
#   03_OW_treatment        ->  04_post_OW           day 3
#   04_DO_treatment        ->  05_post_DO           day 3
#
# Why: the old prefixes mislead. `00_baseline` and `00_laboratory_control` shared the `00`
# prefix but are different phases, and `01_treatment_control` sorted as though it were a
# phase when it is one of the four day-3 arms. The new names put the phase first, so the
# prefix order is meaningful and the phase is readable straight off the folder name.
#
# Uses `git mv`, so history follows the files and the rename is a normal staged change.
#
# DRY RUN BY DEFAULT. Nothing moves until you pass --apply.
#
#   bash thread-strength/01_code/rename_tensometer_folders.sh            # show the plan
#   bash thread-strength/01_code/rename_tensometer_folders.sh --apply    # do it
#
# Run from anywhere inside the repository. Review with `git status` before committing.

set -euo pipefail

APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1

REPO_ROOT=$(git rev-parse --show-toplevel)
BASE="$REPO_ROOT/thread-strength/02_data/tensometer_output"

[ -d "$BASE" ] || { echo "ERROR: not found: $BASE" >&2; exit 1; }

# old:new pairs, applied in order
PAIRS=(
  "00_laboratory_control:00_lab_reference"
  "00_baseline:01_pre_exposure"
  "01_treatment_control:02_post_control"
  "02_OA_treatment:03_post_OA"
  "03_OW_treatment:04_post_OW"
  "04_DO_treatment:05_post_DO"
)

echo "Repository : $REPO_ROOT"
echo "Folder     : thread-strength/02_data/tensometer_output"
[ "$APPLY" -eq 1 ] && echo "Mode       : APPLY" || echo "Mode       : DRY RUN (pass --apply to execute)"
echo

# ---- Pre-flight: check every move before making any of them -----------------------------
PROBLEMS=0
TODO=0
for pair in "${PAIRS[@]}"; do
  OLD="${pair%%:*}"; NEW="${pair##*:}"
  if [ ! -d "$BASE/$OLD" ]; then
    if [ -d "$BASE/$NEW" ]; then
      echo "  skip    $OLD  ->  $NEW   (already renamed)"
    else
      echo "  MISSING $OLD                (source folder not found)"
      PROBLEMS=$((PROBLEMS + 1))
    fi
    continue
  fi
  if [ -e "$BASE/$NEW" ]; then
    echo "  BLOCKED $OLD  ->  $NEW   (destination already exists)"
    PROBLEMS=$((PROBLEMS + 1))
    continue
  fi
  N=$(find "$BASE/$OLD" -maxdepth 1 -type f | wc -l | tr -d ' ')
  echo "  rename  $OLD  ->  $NEW   ($N files)"
  TODO=$((TODO + 1))
done

echo
if [ "$PROBLEMS" -gt 0 ]; then
  echo "$PROBLEMS problem(s) found. Nothing has been changed. Resolve them and re-run." >&2
  exit 1
fi

if [ "$TODO" -eq 0 ]; then
  echo "Nothing to do."
  exit 0
fi

if [ "$APPLY" -eq 0 ]; then
  echo "Dry run only. Re-run with --apply to perform $TODO rename(s)."
  exit 0
fi

# ---- Apply ------------------------------------------------------------------------------
for pair in "${PAIRS[@]}"; do
  OLD="${pair%%:*}"; NEW="${pair##*:}"
  [ -d "$BASE/$OLD" ] || continue
  git -C "$REPO_ROOT" mv "thread-strength/02_data/tensometer_output/$OLD" \
                         "thread-strength/02_data/tensometer_output/$NEW"
  echo "  renamed $OLD -> $NEW"
done

echo
echo "Done. Next:"
echo "  1. git status                       # review the staged renames"
echo "  2. Delete the legacy rows from \`folder_labels\` in 1_extract_tensometer_data.Rmd"
echo "     and from \`folder_arm\` in 0_build_mussel_key.Rmd."
echo "  3. Re-run 1_extract_tensometer_data.Rmd. QC plot folders under"
echo "     03_analyses/extract-tensometer-data/QC_plots/ are regenerated under the new"
echo "     names; the old ones can be deleted."
