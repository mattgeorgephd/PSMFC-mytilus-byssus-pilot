#!/usr/bin/env bash
## One library of step 06 (option B): align the trimmed Tag-seq reads to the genome and
## count them on the three annotations of step 05 with StringTie + prepDE and featureCounts.
## Called by 06_genome_recount.Rmd (and by any batch driver) with the same arguments, so every
## library is processed the same way. Writes into OUT, and OUT/done last; the alignment itself
## is deleted at the end (it is about 200 MB per library and nothing later reads it).
##
## usage: _recount_library.sh READS OUT INDEX SPLICE_SITES ANN_DIR BIN_DIR RSCRIPT THREADS
##   READS         trimmed reads (*_L099_R1_cmb.trim.fastq.gz)
##   OUT           output folder of this library
##   INDEX         HISAT2 index prefix
##   SPLICE_SITES  hisat2_extract_splice_sites.py output for the RefSeq GTF
##   ANN_DIR       step 05's folder (annotation_<a>.gff and annotation_<a>.saf)
##   BIN_DIR       folder holding hisat2, samtools, stringtie and featureCounts
##   RSCRIPT       Rscript binary (for _prepde_sample.R)
##   THREADS       threads
set -euo pipefail
READS=$1; OUT=$2; INDEX=$3; SS=$4; ANN=$5; BIN=$6; RSCRIPT=$7; T=$8
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"
TMP=$(mktemp -d "$OUT/tmp.XXXX")
trap 'rm -rf "$TMP"' EXIT

## HISAT2 as in the previous record (step 01: default settings and --dta), with the RefSeq splice
## sites given at alignment time instead of built into the index
"$BIN/hisat2" -p "$T" --dta --known-splicesite-infile "$SS" -x "$INDEX" -U "$READS" 2> "$OUT/hisat2.log" \
  | "$BIN/samtools" sort -@ 2 -m 1G -T "$TMP/sort" -o "$TMP/aln.bam" -
"$BIN/samtools" index "$TMP/aln.bam"

for a in refseq ext3 full; do
  ## StringTie as in the record (-e against the annotation; -B is left out, it only writes
  ## Ballgown tables), reduced at once to prepDE's gene counts (read length 75)
  "$BIN/stringtie" -e -p "$T" -G "$ANN/annotation_$a.gff" -o "$TMP/st_$a.gtf" "$TMP/aln.bam"
  "$RSCRIPT" "$HERE/_prepde_sample.R" "$TMP/st_$a.gtf" "$OUT/stringtie_$a.tsv.gz"
  rm -f "$TMP/st_$a.gtf"
  ## featureCounts: one count per uniquely aligned read, sense strand (-s 1), reads on more
  ## than one gene left out (the defaults)
  "$BIN/featureCounts" -T "$T" -s 1 -F SAF -a "$ANN/annotation_$a.saf" -o "$TMP/fc_$a.txt" "$TMP/aln.bam" > /dev/null 2> "$TMP/fc_$a.log"
  grep -v '^#' "$TMP/fc_$a.txt" | tail -n +2 | cut -f1,7 | gzip > "$OUT/featurecounts_$a.tsv.gz"
  mv "$TMP/fc_$a.txt.summary" "$OUT/featurecounts_$a.summary"
done
touch "$OUT/done"
