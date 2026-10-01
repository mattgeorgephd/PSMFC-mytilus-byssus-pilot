# Shared paths for the Iso-Seq branch. Sourced by each script in 01_code/; anchored to
# iso-seq-transcriptome.Rproj via here::here(). This folder writes only to its own
# 03_analyses/ (and downloads its large external inputs into 02_data/, git-ignored).

library(here)

repo_root <- normalizePath(file.path(here::here(), ".."))

# Inputs (downloaded, git-ignored)
dat         <- here::here("02_data")
isoseq_fa   <- file.path(dat, "Mtros-hq_transcripts.fasta")      # owl genomic-databank
cds_fa      <- file.path(dat, "cds_from_genomic.fasta")          # GCF_036588685.1 CDS (gannet copy)
ISOSEQ_URL  <- "https://owl.fish.washington.edu/halfshell/genomic-databank/Mtros-hq_transcripts.fasta"
CDS_URL     <- paste0("https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/",
                      "byssus-exp-analysis/data/ncbi_dataset/data/GCF_036588685.1/cds_from_genomic.fasta")
## the genome and its annotation (RefSeq GCF_036588685.1, annotation release RS_2024_02), from NCBI
ASSEMBLY    <- "GCF_036588685.1_PNRI_Mtr1.1.1.hap1"
NCBI_DIR    <- paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/588/685/", ASSEMBLY, "/")
genome_fa   <- file.path(dat, paste0(ASSEMBLY, "_genomic.fna.gz"))
gff_gz      <- file.path(dat, paste0(ASSEMBLY, "_genomic.gff.gz"))
READS_URL   <- paste0("https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/",
                      "byssus-exp-analysis/data/raw-trimmed/")       # the trimmed reads HISAT2 used

# Outputs
map_dir    <- here::here("03_analyses", "02_isoform-gene-map")
cds_map    <- here::here("03_analyses", "_superseded", "02_isoform-gene-map_cds", "isoform_gene_map.csv.gz")  # retired CDS-based map, the cross-check
salmon_dir <- here::here("03_analyses", "03_salmon")
de_dir     <- here::here("03_analyses", "04_isoseq-de")
for (d in c(map_dir, salmon_dir, de_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Cross-folder reads
de06 <- file.path(repo_root, "06_differential-expression", "03_analyses")
