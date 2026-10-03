# Shared paths for the differential-expression analysis.
# Sourced by each script in 01_code/; anchored to differential-expression.Rproj via here::here().
# This folder writes only to its own 03_analyses/; it reads earlier folders (03_blast,
# 05_sequence-alignment) and is read by 04_iso-seq-transcriptome (steps 02, 04, 05, 07), 07, 08
# and 09.

library(here)

repo_root <- normalizePath(file.path(here::here(), ".."))

# Inputs (not written by any script)
dat <- here::here("02_data")                                # raw sample sheet, RNA summary

# Outputs of this analysis
cnt     <- here::here("03_analyses", "count_matrix")        # clean counts + sample table (01)
deg     <- here::here("03_analyses", "DEG_lists")           # contrasts, DESeq2 results, annotated DEGs
dds_dir <- here::here("03_analyses", "dds")                 # fitted DESeq2 objects (03; git-ignored)
top     <- here::here("03_analyses", "top_DEGs")            # top-50 tables per contrast (07)
figs    <- here::here("03_analyses", "figures")             # volcano, MA, PCA, DEG counts
for (d in c(cnt, deg, dds_dir, top, figs)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Outputs of this analysis (continued)
mito    <- here::here("03_analyses", "mitochondrial")       # mitochondrial genes on their own (13)
dir.create(mito, recursive = TRUE, showWarnings = FALSE)

# Cross-folder reads
blast_go <- file.path(repo_root, "03_blast", "03_analyses", "genome-foot")   # LOC_GO_list.txt
t_data   <- file.path(repo_root, "05_sequence-alignment", "03_analyses", "hisat", "t_data.ctab")  # gene -> sequence
mt_annot <- file.path(repo_root, "05_sequence-alignment", "02_data", "annotation_mt_like_loci.csv")  # loci named after a mitochondrial protein
