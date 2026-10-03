#!/bin/bash
## The genome blastx of 01_genome_blast.Rmd, in parts that can be resumed.
##   _blastx_parts.sh <blast_bin> <query.fasta> <db> <out.tab> <threads> <parts> [<jobs>]
## With parts = 1 it runs the 2024 command once. With parts > 1 it splits the query FASTA, in
## order, into that many files of about equal numbers of sequences, runs the same command on each
## part not yet done (writing a temporary file of its own and renaming it to <part>.tab only when
## blastx ends without error, so an interrupted part is redone and never half-kept), <jobs> parts
## at a time with <threads> / <jobs> threads each, and once every part is done writes <out.tab>
## as the parts in order. BLAST reports each
## query independently of the others, so the result is the table one run would write (checked on
## a 99-CDS sample, byte for byte). Rerun the same command to resume after an interruption.
## BLASTX_DEADLINE (optional, seconds since 1970): start no part after this time; the run then
## stops with exit status 3 and is resumed by running the same command again. A second run on
## the same output refuses to start (exit status 4) while one is still going.
set -euo pipefail
bin=$1; query=$2; db=$3; out=$4; threads=$5; parts=$6; jobs=${7:-1}
run_blastx() {  # query, output, threads: the 2024 parameters
  "$bin/blastx" -query "$1" -db "$db" -out "$2" -evalue 1E-20 -num_threads "$3" -max_target_seqs 1 -outfmt 6
}
if [ "$parts" -le 1 ]; then run_blastx "$query" "$out" "$threads"; exit 0; fi
dir="$out.parts"; mkdir -p "$dir"
## one run at a time (a blastx still exiting after an interruption writes only its own file)
exec 9> "$dir/.lock"
flock -n 9 || { echo "another run is still working on $dir" >&2; exit 4; }
rm -f "$dir"/part_*.tab.tmp.*      # parts an interrupted run left unfinished
if [ ! -f "$dir/split.done" ]; then
  n=$(grep -c '^>' "$query")
  awk -v n="$n" -v parts="$parts" -v dir="$dir" '
    /^>/ { i++; p = int((i - 1) * parts / n) + 1; f = sprintf("%s/part_%04d.fa", dir, p) }
    { print > f }' "$query"
  echo "$n sequences in $parts parts" > "$dir/split.done"
fi
per=$(( threads / jobs )); [ "$per" -ge 1 ] || per=1
export -f run_blastx; export bin db; export BLASTX_DEADLINE="${BLASTX_DEADLINE:-}"
ls "$dir"/part_*.fa | while read -r fa; do
  tab="${fa%.fa}.tab"; [ -f "$tab" ] || echo "$fa"
done | xargs -r -P "$jobs" -I{} bash -c 'fa={}; tab="${fa%.fa}.tab"
    if [ -n "$BLASTX_DEADLINE" ] && [ "$(date +%s)" -ge "$BLASTX_DEADLINE" ]; then exit 0; fi
    run_blastx "$fa" "$tab.tmp.$$" '"$per"' && mv "$tab.tmp.$$" "$tab"'
missing=0; for fa in "$dir"/part_*.fa; do [ -f "${fa%.fa}.tab" ] || missing=$((missing + 1)); done
if [ "$missing" -gt 0 ]; then echo "$missing of $parts parts not done; rerun to resume" >&2; exit 3; fi
cat "$dir"/part_*.tab > "$out.tmp" && mv "$out.tmp" "$out"
echo "all $parts parts done: $(wc -l < "$out") hits in $out"
