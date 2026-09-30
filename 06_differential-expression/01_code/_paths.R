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
