## One-off provenance record: the loci of the genome annotation that NCBI names after a
## mitochondrially encoded protein. Not part of the pipeline; it wrote
## 02_data/annotation_mt_like_loci.csv, which tools/mt_encoded.R (06 step 01) reads.
##
## Why it is needed. 06 leaves the mitochondrial loci out of the differential expression and
## tests them per protein. It found the nuclear copies through each gene's best BLAST hit
## (03_blast, a BLAST of the annotation's CDS), which finds the 126 protein-coding copies but
## cannot see pseudogenes, which have no CDS. The RefSeq annotation (GCF_036588685.1, release
## RS_2024_02) also has 167 pseudogenes named after a mitochondrial protein ("cytochrome b-like",
## "ATP synthase subunit a-like", ...), all on unplaced scaffolds: copies of mitochondrial
## sequence (NUMTs) that take mitochondrial reads like the coding copies.
##
## How. Every gene and pseudogene record of the annotation's GFF whose description starts with
## the name of one of the 13 mitochondrially encoded proteins (optionally "-like"), anchored so
## that nuclear genes such as "succinate dehydrogenase [ubiquinone] cytochrome b small subunit"
## or "cytochrome b5" do not match. The mitogenome's own genes carry no description and are
## found by their sequence (NC_007687.1) in tools/mt_encoded.R.
##
## Run from the repository root (needs the network, or pass the GFF's path):
##   Rscript 05_sequence-alignment/01_code/_derive_mt_like_loci.R [path/to/genomic.gff.gz]

source("tools/gene_ids.R")
args <- commandArgs(trailingOnly = TRUE)
url  <- paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/588/685/GCF_036588685.1_PNRI_Mtr1.1.1.hap1/",
               "GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.gff.gz")
gff_gz <- if (length(args)) args[1] else { f <- tempfile(fileext = ".gff.gz"); download.file(url, f, mode = "wb"); f }

g <- read.delim(gzfile(gff_gz), header = FALSE, comment.char = "#", quote = "", stringsAsFactors = FALSE,
                col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"))
g <- g[g$type %in% c("gene", "pseudogene"), ]
attr_of <- function(a, key) {
  v <- ifelse(grepl(paste0("(^|;)", key, "="), a), sub(paste0("^.*?(^|;)", key, "=([^;]*).*$"), "\\2", a, perl = TRUE), NA)
  vapply(v, function(x) if (is.na(x)) NA_character_ else utils::URLdecode(x), "", USE.NAMES = FALSE)
}
g$ID <- attr_of(g$attr, "ID")
g$description <- attr_of(g$attr, "description")
g$gene_biotype <- attr_of(g$attr, "gene_biotype")

mt_rx <- paste0("^(cytochrome c oxidase subunit [123]|cytochrome b|NADH-ubiquinone oxidoreductase chain (1|2|3|4|4L|5|6)|",
                "ATP synthase (subunit a|F0 subunit 6|protein 8))(-like)?(,|$| )")
hit <- !is.na(g$description) & grepl(mt_rx, g$description, ignore.case = TRUE)
d <- tolower(g$description[hit])
sym <- rep(NA_character_, length(d))
sym[grepl("^cytochrome c oxidase subunit [123]", d)] <- paste0("COX", sub("^cytochrome c oxidase subunit ([123]).*$", "\\1", d[grepl("^cytochrome c oxidase subunit [123]", d)]))
sym[grepl("^cytochrome b", d)] <- "CYTB"
sym[grepl("^atp synthase (subunit a|f0 subunit 6)", d)] <- "ATP6"
sym[grepl("^atp synthase protein 8", d)] <- "ATP8"
nd <- grepl("^nadh-ubiquinone oxidoreductase chain", d)
sym[nd] <- paste0("ND", toupper(sub("^nadh-ubiquinone oxidoreductase chain (4l|[1-6]).*$", "\\1", d[nd])))
stopifnot(!anyNA(sym))

out <- data.frame(LOC_ID = gene_key(g$ID[hit]), type = g$type[hit], gene_biotype = g$gene_biotype[hit],
                  description = g$description[hit], mt_protein = sym, seqid = g$seqid[hit],
                  start = g$start[hit], end = g$end[hit], strand = g$strand[hit])
out <- out[order(out$type, out$mt_protein, out$LOC_ID), ]
stopifnot(!anyDuplicated(out$LOC_ID))
write.csv(out, "05_sequence-alignment/02_data/annotation_mt_like_loci.csv", row.names = FALSE)
print(table(out$type, out$mt_protein))
