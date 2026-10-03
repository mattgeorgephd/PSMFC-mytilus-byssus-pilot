# Shared paths for the Iso-Seq sensitivity analysis. Sourced by each script in 01_code/;
# anchored to iso-seq-transcriptome.Rproj via here::here(). This folder writes only to its own
# 03_analyses/ (and, when needed, downloads the Iso-Seq transcriptome into 02_data/, git-ignored).

library(here)

repo_root <- normalizePath(file.path(here::here(), ".."))
seq_dir   <- file.path(repo_root, "04_sequence-alignment")

# Inputs
## the Iso-Seq transcriptome (owl): 04_sequence-alignment step 04 downloads it; this folder uses
## that copy, or its own in 02_data/ (step 01 downloads one; step 02 needs it only to build its index)
ISOSEQ_URL <- "https://owl.fish.washington.edu/halfshell/genomic-databank/Mtros-hq_transcripts.fasta"
isoseq_fa  <- file.path(seq_dir, "02_data", "Mtros-hq_transcripts.fasta")
if (!file.exists(isoseq_fa)) isoseq_fa <- here::here("02_data", "Mtros-hq_transcripts.fasta")
READS_URL  <- paste0("https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/",
                     "byssus-exp-analysis/data/raw-trimmed/")        # the trimmed reads HISAT2 used
## outputs of 04_sequence-alignment steps 04 to 06
map_dir     <- file.path(seq_dir, "03_analyses", "isoform-gene-map")
ann_dir     <- file.path(seq_dir, "03_analyses", "augmented-annotation")
recount_dir <- file.path(seq_dir, "03_analyses", "genome-recount")
## differential expression (steps 03 and 04)
de05 <- file.path(repo_root, "05_differential-expression", "03_analyses")

# Outputs
salmon_dir <- here::here("03_analyses", "02_salmon")
de_dir     <- here::here("03_analyses", "03_isoseq-de")
aug_de_dir <- here::here("03_analyses", "04_augmented-de")
for (d in c(salmon_dir, de_dir, aug_de_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
