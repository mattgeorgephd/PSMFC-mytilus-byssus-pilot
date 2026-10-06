# Shared paths for steps 04 to 06 (the isoform-gene map, the augmented annotation and the genome
# recount, which until 2026-10-03 were steps 02, 05 and 06 of 06_iso-seq-transcriptome). The other
# steps set their paths in their own paths chunk. Anchored to sequence-alignment.Rproj via
# here::here(). Writes only to this folder's 03_analyses/, and downloads the large external
# inputs into 02_data/ (git-ignored).

library(here)

repo_root <- normalizePath(file.path(here::here(), ".."))

# Inputs (downloaded, git-ignored)
dat         <- here::here("02_data")
isoseq_fa   <- file.path(dat, "Mtros-hq_transcripts.fasta")      # owl genomic-databank
cds_fa      <- file.path(dat, "cds_from_genomic.fasta")          # GCF_036588685.1 CDS (gannet copy; the retired CDS map)
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

## the annotation's mitochondrial loci (steps 04 and 05, through tools/mt_encoded.R): 03 step 05's
## table, which is 03 step 01's (the one the committed outputs of steps 04 and 05 were made from)
## plus six byssal genes, so the mitochondrial loci are the same
blast_go <- file.path(repo_root, "03_blast", "03_analyses", "genome-foot-sprot2026_03-noseg", "LOC_GO_list.txt")
t_data   <- here::here("03_analyses", "hisat", "t_data.ctab")
mt_annot <- here::here("02_data", "annotation_mt_like_loci.csv")

# Outputs
map_dir     <- here::here("03_analyses", "isoform-gene-map")
cds_map     <- here::here("03_analyses", "_superseded", "02_isoform-gene-map_cds", "isoform_gene_map.csv.gz")  # retired CDS-based map, the cross-check
ann_dir     <- here::here("03_analyses", "augmented-annotation")
recount_dir <- here::here("03_analyses", "genome-recount")
for (d in c(map_dir, ann_dir, recount_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
