#!/bin/bash
## The search of 05_byssal_noseg_search.Rmd: every genome CDS against the 241 foot and byssal
## proteins alone, with blastx's low-complexity filter off.
##   _blastx_noseg.sh <blast_bin> <query.fasta> <foot.fasta> <out.tab> <threads> <dbsize>
## The settings are those of step 01 (-evalue 1E-20 -max_target_seqs 1 -outfmt 6) plus -seg no,
## and -dbsize set to the total residues of step 01's database, so that the e-values (and the
## 1E-20 cut-off) are those a search of step 01's whole database would give. The database is made
## next to <out.tab>. The table is written to a temporary file and renamed only when blastx ends
## without error; if <out.tab> exists the search is not run again.
set -euo pipefail
bin=$1; query=$2; foot=$3; out=$4; threads=$5; dbsize=$6
if [ -f "$out" ]; then echo "$out exists ($(wc -l < "$out") hits); not searched again"; exit 0; fi
db="${out%.tab}_footdb"
"$bin/makeblastdb" -in "$foot" -dbtype prot -out "$db" > /dev/null
"$bin/blastx" -query "$query" -db "$db" -out "$out.tmp" -evalue 1E-20 -num_threads "$threads" \
  -max_target_seqs 1 -seg no -dbsize "$dbsize" -outfmt 6
mv "$out.tmp" "$out"
echo "search done: $(wc -l < "$out") hits in $out"
