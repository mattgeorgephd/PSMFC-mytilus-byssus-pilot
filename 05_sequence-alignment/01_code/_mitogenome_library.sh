#!/usr/bin/env bash
## One library of step 05: the trimmed Tag-seq reads aligned to the mitochondrial genome alone
## and counted per mitochondrial gene. Called by 05_mitogenome_counts.Rmd (and by any batch
## driver) with the same arguments, so every library is processed the same way. Writes into OUT,
## and OUT/done last; the alignments are deleted at the end.
##
## usage: _mitogenome_library.sh READS OUT INDEX SAF BIN_DIR THREADS
##   READS    trimmed reads (*_L099_R1_cmb.trim.fastq.gz)
##   OUT      output folder of this library
##   INDEX    HISAT2 index prefix of 02_data/mitogenome_NC_007687.1.fa
##   SAF      02_data/mitogenome_genes.saf
##   BIN_DIR  folder holding hisat2, samtools and featureCounts
##   THREADS  threads
set -euo pipefail
READS=$1; OUT=$2; INDEX=$3; SAF=$4; BIN=$5; T=$6
mkdir -p "$OUT"
TMP=$(mktemp -d "$OUT/tmp.XXXX")
trap 'rm -rf "$TMP"' EXIT

## reads in the library and their mean length
zcat "$READS" | awk 'NR % 4 == 2 {n++; s += length($0)} END {printf "%d\t%.2f\n", n, s / n}' > "$OUT/reads.tsv"

## of record: HISAT2 defaults (minimum score -0.2 x read length, about two mismatches in a
## 63-base read), no spliced alignment (the mitogenome is not spliced); and, as a sensitivity
## setting, a permissive minimum score (-0.6 x read length, about six mismatches)
for v in default permissive; do
  extra=""; [ "$v" = permissive ] && extra="--score-min L,0,-0.6"
  "$BIN/hisat2" -p "$T" --no-spliced-alignment --no-unal $extra -x "$INDEX" -U "$READS" 2> "$OUT/hisat2_$v.log" \
    | "$BIN/samtools" sort -@ 1 -T "$TMP/sort_$v" -o "$TMP/$v.bam" -
  ## one count per read, sense strand; the few reads overlapping two adjacent genes, or aligned
  ## to two places, are shared fractionally
  "$BIN/featureCounts" -T "$T" -s 1 -M -O --fraction -F SAF -a "$SAF" -o "$TMP/fc_$v.txt" "$TMP/$v.bam" \
    > /dev/null 2> "$TMP/fc_$v.log"
  grep -v '^#' "$TMP/fc_$v.txt" | tail -n +2 | cut -f1,7 > "$OUT/counts_$v.tsv"
  mv "$TMP/fc_$v.txt.summary" "$OUT/featurecounts_$v.summary"
done
touch "$OUT/done"
