# Shared paths for the GO enrichment analysis.
# Sourced by the numbered scripts in 01_code/; anchored to enrichment.Rproj via here::here().
# This folder only reads from other folders and writes to its own 03_analyses/.

library(here)

repo_root <- normalizePath(file.path(here::here(), ".."))

# Read from other folders (never written here)
deg      <- file.path(repo_root, "06_differential-expression", "03_analyses", "DEG_lists")    # contrasts + apeglm tables
blast_go <- file.path(repo_root, "03_blast", "03_analyses", "genome-foot")                  # LOC_GO_list.txt
t_data   <- file.path(repo_root, "05_sequence-alignment", "03_analyses", "hisat", "t_data.ctab")  # transcript lengths

# This analysis's own outputs, one folder per script
out_inputs <- here::here("03_analyses", "01_go-inputs")
out_topgo  <- here::here("03_analyses", "02_topgo")
out_goseq  <- here::here("03_analyses", "03_goseq")
out_cp     <- here::here("03_analyses", "04_clusterprofiler")
out_rrvgo  <- here::here("03_analyses", "05_rrvgo")
out_cmp    <- here::here("03_analyses", "06_method-comparison")
for (d in c(out_inputs, out_topgo, out_goseq, out_cp, out_rrvgo, out_cmp))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
