## One-off provenance record: the mitochondrial genome and its genes, the reference of step 08
## (the mitochondrial reads counted on the mitogenome alone). Not part of the pipeline; it wrote
## 02_data/mitogenome_NC_007687.1.fa and 02_data/mitogenome_genes.saf.
##
## What. The genome assembly (RefSeq GCF_036588685.1) holds the female-type mitochondrial
## genome of M. trossulus as NC_007687.1 (18,652 bases). The FASTA is that one sequence, cut out
## of the assembly; the SAF lists its 38 annotated features (annotation release RS_2024_02):
## the 12 protein genes (gene records, named by symbol: ATP6, COX1-3, CYTB, ND1-6, ND4L), the 2
## rRNAs and the 24 tRNAs (which have no gene record and are named by their RNA ID,
## rna-NC_007687.1:1..1244). The names are the gene keys of the count matrix of record
## (tools/gene_ids.R).
##
## Run from the repository root (needs the network, or pass the FASTA's and the GFF's paths):
##   Rscript 04_sequence-alignment/01_code/_derive_mitogenome.R [genomic.fna.gz genomic.gff.gz]

args <- commandArgs(trailingOnly = TRUE)
ncbi <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/588/685/GCF_036588685.1_PNRI_Mtr1.1.1.hap1/"
get <- function(i, file) {
  if (length(args) >= i) return(args[i])
  f <- tempfile(fileext = sub("^.*?(\\..*)$", "\\1", file)); download.file(paste0(ncbi, file), f, mode = "wb"); f
}
fna_gz <- get(1, "GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.fna.gz")
gff_gz <- get(2, "GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.gff.gz")
MITO <- "NC_007687.1"
out  <- file.path("04_sequence-alignment", "02_data")

## the sequence: the record from its header line to the next header
con <- gzfile(fna_gz, "r"); keep <- FALSE; seq_lines <- character(0); header <- NA_character_
repeat {
  chunk <- readLines(con, n = 1e6)
  if (!length(chunk)) break
  for (ln in chunk) {
    if (startsWith(ln, ">")) {
      if (keep) break
      keep <- startsWith(ln, paste0(">", MITO, " ")); if (keep) header <- ln
    } else if (keep) seq_lines <- c(seq_lines, ln)
  }
  if (keep && length(seq_lines) && any(startsWith(chunk, ">") & !startsWith(chunk, paste0(">", MITO, " ")))) break
}
close(con)
s <- paste(seq_lines, collapse = "")
stopifnot(!is.na(header), nchar(s) == 18652)
writeLines(c(header, substring(s, seq(1, nchar(s), 80), pmin(seq(80, nchar(s) + 79, 80), nchar(s)))),
           file.path(out, "mitogenome_NC_007687.1.fa"))

## the features: protein genes by their gene record, rRNAs and tRNAs by their RNA record
g <- read.delim(gzfile(gff_gz), header = FALSE, comment.char = "#", quote = "", stringsAsFactors = FALSE,
                col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"))
g <- g[g$seqid == MITO, ]
id   <- sub("^ID=([^;]*).*$", "\\1", g$attr)
name <- ifelse(grepl("(^|;)Name=", g$attr), sub("^.*?(^|;)Name=([^;]*).*$", "\\2", g$attr, perl = TRUE), NA)
genes <- g[g$type == "gene", ]
rnas  <- g[g$type %in% c("rRNA", "tRNA"), ]
saf <- rbind(data.frame(GeneID = name[g$type == "gene"], Chr = MITO, Start = genes$start, End = genes$end, Strand = genes$strand),
             data.frame(GeneID = id[g$type %in% c("rRNA", "tRNA")], Chr = MITO, Start = rnas$start, End = rnas$end, Strand = rnas$strand))
saf <- saf[order(saf$Start), ]
stopifnot(nrow(saf) == 38, sum(g$type == "gene") == 12, !anyDuplicated(saf$GeneID))
write.table(saf, file.path(out, "mitogenome_genes.saf"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", file.path(out, "mitogenome_NC_007687.1.fa"), "and", nrow(saf), "features\n")
