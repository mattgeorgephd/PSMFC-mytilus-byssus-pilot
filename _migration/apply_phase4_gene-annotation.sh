#!/usr/bin/env bash
# Phase 4 (gene-annotation): repoint 06-get_GOSlims, 17-uniprot_summaries, 18-ortholog-lists.
# Annotation.Rmd is left UNCHANGED (HPC blast archive). 06's destructive awk/sed bash chunk is
# disabled (eval=FALSE) and path-repointed via the GOTERMS env var. Run from repo ROOT on a branch.
set -euo pipefail
[ -d .git ] || { echo "Run from the repo root."; exit 1; }
command -v perl >/dev/null 2>&1 || { echo "ERROR: perl not found (it normally ships with Git)."; exit 1; }
D=gene-annotation
[ -d "$D/01_code" ] || { echo "ERROR: $D/01_code missing. Current clone of main?"; exit 1; }
[ -e "$D/01_code/_paths.R" ] && { echo "ERROR: _paths.R already exists; already applied? git reset --hard && git clean -fd"; exit 1; }

cat > "$D/01_code/_paths.R" <<'RPATHS'
# Shared paths for the gene-annotation analysis.
# Sourced by 06/17/18 in 01_code/; anchored to gene-annotation.Rproj via here::here().
# (Annotation.Rmd is an unchanged HPC blast archive and does not source this file.)

library(here)

repo_root <- normalizePath(file.path(here::here(), ".."))

# Cross-folder reads/writes into the differential-expression DEG_lists
deg <- file.path(repo_root, "differential-expression", "03_analyses", "DEG_lists")

# This analysis's own outputs
topgenes <- here::here("03_analyses", "Top_gene_summaries")
dir.create(topgenes, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(topgenes, "Top_50_genes"), recursive = TRUE, showWarnings = FALSE)

# Path for the (disabled) awk/sed bash chunk in 06-get_GOSlims. The bash subprocess inherits
# this env var. The chunk is eval=FALSE because it rewrites a committed CSV in place and its
# tab-based logic no longer matches the current comma-separated file format.
Sys.setenv(GOTERMS = file.path(deg, "GOterms_genome"))
RPATHS

perl - "$D/01_code" <<'PERL'
use strict; use warnings;
my $base = shift;
my @scripts = qw(06-get_GOSlims 17-uniprot_summaries 18-ortholog-lists);
my @prefixes = ("/home/shared/8TB_HDD_02/graceleuchtenberger/Github/byssus-exp-analysis/",
                "/home/shared/8TB_HDD_02/graceleuchtenberger/byssus-exp-analysis/");
my $chunk = "\n```{r paths, include=FALSE}\nsource(here::here(\"01_code\", \"_paths.R\"))\n```\n";
for my $name (@scripts) {
  my $f = "$base/$name.Rmd";
  local $/; open(my $in,'<',$f) or die "open $f: $!"; my $t = <$in>; close $in;
  $t =~ s{\Q$_\E}{}g for @prefixes;
  $t =~ s{"\.\./output/}{"output/}g; $t =~ s{"\.\./data/}{"data/}g;
  $t =~ s{"output/Top_gene_summaries/([^"]*)"}{file.path(topgenes, "$1")}g;
  $t =~ s{"output/DEG_lists/([^"]*)"}{file.path(deg, "$1")}g;
  if ($name eq '06-get_GOSlims') {
    $t =~ s/```\{bash rm-spaces-from-GOs\}/```{bash rm-spaces-from-GOs, eval=FALSE}/;
    $t =~ s{(```\{bash[^\n]*\n.*?```)}{ my $b=$1; $b =~ s|output/DEG_lists/GOterms_genome/|\${GOTERMS}/|g; $b }gse;
  }
  unless ($t =~ s/\A(---\r?\n.*?\n---[ \t]*\r?\n)/$1$chunk/s) { $t = $chunk."\n".$t; }
  open(my $out,'>',$f) or die "write $f: $!"; print $out $t; close $out;
  print "  repointed $name\n";
}
PERL
echo "Phase 4 (gene-annotation) applied. (Annotation.Rmd left unchanged as HPC archive.)"
