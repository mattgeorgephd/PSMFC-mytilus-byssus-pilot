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
