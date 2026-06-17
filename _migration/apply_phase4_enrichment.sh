#!/usr/bin/env bash
# Phase 4 (enrichment): repoint the 3 absorbed GO-enrichment scripts. All pure R.
# Reads cross-folder from differential-expression (DEG_lists, 02_data) and blast (genome-foot);
# writes its own Func_annot_DAVID. Run from repo ROOT on a branch. Uses Perl (ships with Git).
set -euo pipefail
[ -d .git ] || { echo "Run from the repo root."; exit 1; }
command -v perl >/dev/null 2>&1 || { echo "ERROR: perl not found (it normally ships with Git)."; exit 1; }
D=enrichment
[ -d "$D/01_code" ] || { echo "ERROR: $D/01_code missing. Current clone of main?"; exit 1; }
[ -e "$D/01_code/_paths.R" ] && { echo "ERROR: _paths.R already exists; already applied? git reset --hard && git clean -fd"; exit 1; }

cat > "$D/01_code/_paths.R" <<'RPATHS'
# Shared paths for the enrichment analysis.
# Sourced by each script in 01_code/; anchored to enrichment.Rproj via here::here().

library(here)

repo_root <- normalizePath(file.path(here::here(), ".."))

# Cross-folder reads (the absorbed GO-enrichment pipeline spans several analysis folders)
deg      <- file.path(repo_root, "differential-expression", "03_analyses", "DEG_lists")  # sig lists, DAVID/REVIGO input lists
dat      <- file.path(repo_root, "differential-expression", "02_data")                   # gene_count_matrix_clean
blast_go <- file.path(repo_root, "blast", "03_analyses", "genome-foot")                   # g.spid.txt, LOC_GO_list.txt

# This analysis's own outputs
func_david <- here::here("03_analyses", "Func_annot_DAVID")                               # DAVID chart .xlsx
dir.create(func_david, recursive = TRUE, showWarnings = FALSE)
RPATHS

perl - "$D/01_code" <<'PERL'
use strict; use warnings;
my $base = shift;
my @scripts = qw(07-GOenrichment_listcreationforDAVID 08-func_enrichment_DAVID_visual 09-GOenrichment_listcreation_REVIGO);
my @prefixes = ("/home/shared/8TB_HDD_02/graceleuchtenberger/Github/byssus-exp-analysis/",
                "/home/shared/8TB_HDD_02/graceleuchtenberger/byssus-exp-analysis/");
my $chunk = "\n```{r paths, include=FALSE}\nsource(here::here(\"01_code\", \"_paths.R\"))\n```\n";
for my $name (@scripts) {
  my $f = "$base/$name.Rmd";
  local $/; open(my $in,'<',$f) or die "open $f: $!"; my $t = <$in>; close $in;
  $t =~ s{\Q$_\E}{}g for @prefixes;
  $t =~ s{"\.\./output/}{"output/}g; $t =~ s{"\.\./data/}{"data/}g;
  $t =~ s{"output/gene_count_matrix_clean\.csv"}{file.path(dat, "gene_count_matrix_clean.csv")}g;
  $t =~ s{"output/gene_count_matrix_clean"}{file.path(dat, "gene_count_matrix_clean")}g;
  $t =~ s{"output/g\.spid\.txt"}{file.path(blast_go, "g.spid.txt")}g;
  $t =~ s{"output/LOC_GO_list\.txt"}{file.path(blast_go, "LOC_GO_list.txt")}g;
  $t =~ s{"output/Func_annot_DAVID/([^"]*)"}{file.path(func_david, "$1")}g;
  $t =~ s{"output/DEG_lists/([^"]*)"}{file.path(deg, "$1")}g;
  unless ($t =~ s/\A(---\r?\n.*?\n---[ \t]*\r?\n)/$1$chunk/s) { $t = $chunk."\n".$t; }
  open(my $out,'>',$f) or die "write $f: $!"; print $out $t; close $out;
  print "  repointed $name\n";
}
PERL
echo "Phase 4 (enrichment) applied."
