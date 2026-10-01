# Shared paths for the gene-annotation analysis.
# Sourced by the numbered scripts in 01_code/; anchored to gene-annotation.Rproj via here::here().
# This folder only reads from other folders and writes to its own 03_analyses/.

library(here)

repo_root <- normalizePath(file.path(here::here(), ".."))
dat       <- here::here("02_data")

# Read from 06_differential-expression (never written here)
deg      <- file.path(repo_root, "06_differential-expression", "03_analyses", "DEG_lists")
top50_in <- file.path(repo_root, "06_differential-expression", "03_analyses", "top_DEGs", "Top_50_genes")

# This analysis's own outputs
goslims  <- here::here("03_analyses", "goslims")
topgenes <- here::here("03_analyses", "Top_gene_summaries")
for (d in c(goslims, topgenes)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
