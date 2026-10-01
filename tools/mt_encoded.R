## Mitochondrial loci in the count matrix. Base R.
##
## The alignment reference holds the female-type mitochondrial genome of M. trossulus
## (NC_007687.1): 12 protein genes (ATP6, COX1-3, CYTB, ND1-6, ND4L; mytilids have no ATP8)
## and its tRNAs and rRNAs. About 126 further LOCs, all on unplaced scaffolds, have one of
## those proteins as best UniProt hit, at the same identity and alignment length as the
## mitogenome genes themselves (about 23 identical copies each of CYTB, COX2, COX3, ND2 and
## ND6): duplicated mitochondrial sequence in the assembly. Reads from one mitochondrial
## transcript are split between the mitogenome gene and its copies, so one mitochondrial signal
## counted many times over (in Gill OA, 116 of the 711 TC DEGs were mitochondrial loci before they were
## separated).
##
## These loci are kept out of the genome differential expression (06 steps 03 and 05), the GO
## enrichment (07) and the gene-mechanics tests (09), and are analysed on their own, one
## summed count per mitochondrial protein, in 06 step 13.

MITOGENOME_SEQID <- "NC_007687.1"

MT_ENCODED_RX <- paste0("^(Cytochrome c oxidase subunit [123] |Cytochrome b \\(|",
                        "NADH-ubiquinone oxidoreductase chain (1|2|3|4|4L|5|6) |",
                        "ATP synthase subunit a \\(|ATP synthase protein 8)")
is_mt_encoded <- function(protein_name) !is.na(protein_name) & grepl(MT_ENCODED_RX, protein_name)

## gene symbol of a mitochondrially encoded protein from its UniProt protein name
mt_protein_symbol <- function(protein_name) {
  p <- ifelse(is.na(protein_name), "", protein_name)
  out <- rep(NA_character_, length(p))
  out[grepl("^Cytochrome c oxidase subunit 1 ", p)] <- "COX1"
  out[grepl("^Cytochrome c oxidase subunit 2 ", p)] <- "COX2"
  out[grepl("^Cytochrome c oxidase subunit 3 ", p)] <- "COX3"
  out[grepl("^Cytochrome b \\(", p)]                <- "CYTB"
  out[grepl("^ATP synthase subunit a \\(", p)]      <- "ATP6"
  out[grepl("^ATP synthase protein 8", p)]          <- "ATP8"
  nd <- regmatches(p, regexpr("^NADH-ubiquinone oxidoreductase chain (4L|[1-6]) ", p))
  out[grepl("^NADH-ubiquinone oxidoreductase chain (4L|[1-6]) ", p)] <-
    paste0("ND", sub("^NADH-ubiquinone oxidoreductase chain (4L|[1-6]) $", "\\1", nd))
  out
}

## Every mitochondrial row of a count matrix, with its source and protein.
##   gene_ids   : count-matrix row names, e.g. "gene-COX1|COX1", "STRG.10|LOC134702910",
##                "rna-NC_007687.1:10108..10170"
##   blast_file : 03_blast/03_analyses/genome-foot/LOC_GO_list.txt (best hit = highest bitscore)
##   t_data_file: 05_sequence-alignment/03_analyses/hisat/t_data.ctab (sequence of each gene)
## Returns one row per mitochondrial matrix row: gene, LOC_ID, source ("mitogenome" or
## "nuclear copy"), mt_protein (symbol, or "tRNA/rRNA"), protein_name, pident, seqid.
## Needs gene_key() from tools/gene_ids.R.
mitochondrial_loci <- function(gene_ids, blast_file, t_data_file) {
  if (!exists("gene_key", mode = "function"))
    stop("source tools/gene_ids.R before calling mitochondrial_loci()")
  key  <- gene_key(gene_ids)                               # LOC ID, or the mitogenome gene name
  head <- sub("[|].*", "", gene_ids)                       # gene-COX1, STRG.10, rna-NC_007687.1:...
  td <- utils::read.delim(t_data_file, stringsAsFactors = FALSE)
  mito_genes <- unique(td$gene_id[td$chr == MITOGENOME_SEQID & td$gene_id != "."])
  mito_rnas  <- unique(td$t_name[td$chr == MITOGENOME_SEQID & td$gene_id == "."])
  on_mito <- head %in% c(mito_genes, mito_rnas)

  b <- utils::read.delim(blast_file, stringsAsFactors = FALSE, check.names = FALSE)
  b <- b[order(-b$bitscore), ]
  b <- b[!duplicated(b$gene), c("gene", "Protein.names", "pident")]
  hit <- b[match(key, b$gene), ]
  mt_hit <- is_mt_encoded(hit$Protein.names)

  keep <- on_mito | mt_hit
  seqid <- td$chr[match(paste0("gene-", key), td$gene_id)]
  seqid[on_mito] <- MITOGENOME_SEQID
  out <- data.frame(gene = gene_ids, LOC_ID = key,
                    source = ifelse(on_mito, "mitogenome", "nuclear copy"),
                    mt_protein = ifelse(mt_hit, mt_protein_symbol(hit$Protein.names),
                                        ifelse(grepl("^rna-", head), "tRNA/rRNA", NA_character_)),
                    protein_name = hit$Protein.names, pident = hit$pident, seqid = seqid,
                    stringsAsFactors = FALSE)[keep, ]
  rownames(out) <- NULL
  out[order(out$source, out$mt_protein, out$gene), ]
}
