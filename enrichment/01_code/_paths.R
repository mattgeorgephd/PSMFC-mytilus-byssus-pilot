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
