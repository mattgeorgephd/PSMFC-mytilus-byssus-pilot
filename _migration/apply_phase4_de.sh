#!/usr/bin/env bash
# Phase 4 (differential-expression): repoint the absorbed DESeq2/DEG scripts to here()/02_data/
# 03_analyses with cross-folder reads via repo_root, add a shared _paths.R, and move the
# kallisto-era 01-diff-exp-analysis to _superseded. Run from repo ROOT on a branch.
# Uses Perl (bundled with Git for Windows) - no Python needed.
set -euo pipefail
[ -d .git ] || { echo "Run from the repo root."; exit 1; }
command -v perl >/dev/null 2>&1 || { echo "ERROR: perl not found (it normally ships with Git)."; exit 1; }
D=differential-expression
[ -d "$D/01_code" ] || { echo "ERROR: $D/01_code missing. Current clone of main?"; exit 1; }
[ -e "$D/01_code/_paths.R" ] && { echo "ERROR: _paths.R already exists; already applied? git reset --hard && git clean -fd"; exit 1; }

cat > "$D/01_code/_paths.R" <<'RPATHS'
# Shared paths for the differential-expression analysis.
# Sourced by each script in 01_code/; anchored to differential-expression.Rproj via here::here().

library(here)

repo_root <- normalizePath(file.path(here::here(), ".."))

# Inputs and outputs within this analysis
dat <- here::here("02_data")                      # count matrices, treatment design, sample metadata
deg <- here::here("03_analyses", "DEG_lists")     # DEG tables (written here; some read back as intermediates)
dir.create(deg, recursive = TRUE, showWarnings = FALSE)

# Cross-folder reads/writes (absorbed pipeline spans several analysis folders)
blast_go <- file.path(repo_root, "blast", "03_analyses", "genome-foot")              # LOC_GO_list.txt, g.spid.txt
topgenes <- file.path(repo_root, "gene-annotation", "03_analyses", "Top_gene_summaries")
RPATHS

perl - "$D/01_code" <<'PERL'
use strict; use warnings;
my $base = shift;
my @scripts = qw(01_5-gene_count_matrix 02_5_DESeq_Foot_TC_genome 02_5_DESeq_Gill_LC_genome
  02_5_DESeq_Gill_TC_genome 02_5_DESeq_foot_LC_genome 02_5_DESeq_tissuevtissue_control
  03-LC_Shrinkage_filtration 03-TC_shrinkage_filtration 04-File_joining
  12-DEG_venn 12-Volcano-plots 15-number_DEGS 16-top_DEGs 19-DEG_list_cleanup);
my @prefixes = ("/home/shared/8TB_HDD_02/graceleuchtenberger/Github/byssus-exp-analysis/",
                "/home/shared/8TB_HDD_02/graceleuchtenberger/byssus-exp-analysis/");
my $chunk = "\n```{r paths, include=FALSE}\nsource(here::here(\"01_code\", \"_paths.R\"))\n```\n";
for my $name (@scripts) {
  my $f = "$base/$name.Rmd";
  local $/; open(my $in,'<',$f) or die "open $f: $!"; my $t = <$in>; close $in;
  $t =~ s{\Q$_\E}{}g for @prefixes;
  $t =~ s{"\.\./output/}{"output/}g; $t =~ s{"\.\./data/}{"data/}g;
  for my $fn ("gene_count_matrix_clean.csv","gene_count_matrix.csv","transcript_count_matrix.csv","treatmentinfo_clean.csv") {
    $t =~ s{"output/\Q$fn\E"}{file.path(dat, "$fn")}g;
  }
  $t =~ s{"data/PSMFC-mytilus-byssus-pilot-RNA-tagseq_raw\.csv"}{file.path(dat, "PSMFC-mytilus-byssus-pilot-RNA-tagseq_raw.csv")}g;
  $t =~ s{"output/LOC_GO_list\.txt"}{file.path(blast_go, "LOC_GO_list.txt")}g;
  $t =~ s{"output/g\.spid\.txt"}{file.path(blast_go, "g.spid.txt")}g;
  $t =~ s{"output/DEG_lists/([^"]*)"}{file.path(deg, "$1")}g;
  $t =~ s{"output/Top_gene_summaries/([^"]*)"}{file.path(topgenes, "$1")}g;
  unless ($t =~ s/\A(---\r?\n.*?\n---[ \t]*\r?\n)/$1$chunk/s) { $t = $chunk."\n".$t; }
  open(my $out,'>',$f) or die "write $f: $!"; print $out $t; close $out;
  print "  repointed $name\n";
}
PERL

mkdir -p "$D/01_code/_superseded"
git mv "$D/01_code/01-diff-exp-analysis.Rmd" "$D/01_code/_superseded/01-diff-exp-analysis.Rmd"
echo "Phase 4 (differential-expression) applied."
