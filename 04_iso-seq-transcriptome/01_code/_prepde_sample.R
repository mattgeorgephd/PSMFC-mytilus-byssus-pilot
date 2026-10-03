## prepDE's gene counts for one StringTie sample, fast enough for the large annotations of
## 04 step 06. Same rule as the port in 05_sequence-alignment/01_code/_prepde.R (StringTie's
## prepDE.py3): for each transcript, count = ceiling(cov * transcript length / read length),
## a gene's count is the sum of its transcripts', and the gene is named "gene_id|gene_name"
## (leftmost matches, so StringTie's ref_gene_name counts as the name), else gene_id. A
## transcript with no exon length is skipped unless it is the file's last, and a repeated
## transcript ID keeps its last entry, as prepDE does. Step 06 checks this against the port.
##
## stable_strg: StringTie puts reference transcripts that have no gene record above them (about
## 300 RefSeq genes whose exons hang directly from the gene, pseudogenes for instance) under
## gene IDs of its own, STRG.N, numbered in the order it meets them, which changes from one
## library to the next. prepDE fixes them from the first sample; counting one library at a time
## cannot, so with stable_strg = TRUE those transcripts are named by their transcript ID instead
## ("<transcript_id>|<gene_name>", e.g. gene-LOC134713724|LOC134713724), which gene_key() reads
## the same way. With FALSE the names are prepDE's exactly.
##
## As a script: Rscript _prepde_sample.R <stringtie.gtf> <out.tsv.gz> [read_length]
## writes two columns, gene and count, without a header, with stable_strg = TRUE.

prepde_sample <- function(gtf, read_length = 75L, stable_strg = TRUE) {
  d <- data.table::fread(cmd = sprintf("grep -v '^#' %s | cut -f3,4,5,9", shQuote(gtf)), sep = "\t",
                         header = FALSE, quote = "", col.names = c("feat", "start", "end", "attr"),
                         colClasses = c("character", "integer", "integer", "character"))
  is_t <- d$feat == "transcript"
  block <- cumsum(is_t)
  if (any(d$feat == "exon" & block == 0)) stop("exon line before the first transcript line in ", gtf)
  ex <- d$feat == "exon"
  tlen <- numeric(sum(is_t))
  if (any(ex)) { s <- rowsum(as.numeric(d$end[ex] - d$start[ex] + 1L), block[ex]); tlen[as.integer(rownames(s))] <- s[, 1] }
  first <- function(rx, x) {                                   # leftmost capture, as re.search
    m <- regexpr(rx, x, perl = TRUE)
    st <- attr(m, "capture.start")[, 1]; ln <- attr(m, "capture.length")[, 1]
    ifelse(m > 0, substr(x, st, st + ln - 1L), NA_character_)
  }
  a <- d$attr[is_t]
  tid <- first('transcript_id "([^"]+)"', a)
  if (anyNA(tid)) stop("a transcript line without transcript_id in ", gtf)
  gid <- first('gene_id "([^"]+)"', a); gnm <- first('gene_name "([^"]+)"', a)
  gene <- ifelse(is.na(gid), tid, ifelse(is.na(gnm), gid, paste0(gid, "|", gnm)))
  if (stable_strg) {
    strg <- !is.na(gid) & startsWith(gid, "STRG.")
    gene[strg] <- ifelse(is.na(gnm[strg]), tid[strg], paste0(tid[strg], "|", gnm[strg]))
  }
  cov <- suppressWarnings(as.numeric(first('cov "([-+0-9.]+)"', a)))
  cov[is.na(cov)] <- 0; cov <- pmax(cov, 0)
  keep <- tlen > 0 | seq_along(tlen) == length(tlen)
  n <- as.integer(ceiling(cov * tlen / read_length))
  gene_of <- setNames(gene, tid)[!duplicated(tid, fromLast = TRUE)]   # a transcript's last definition
  kt <- tid[keep]; kn <- n[keep]
  lk <- !duplicated(kt, fromLast = TRUE)                       # a repeated ID keeps its last entry
  s <- rowsum(as.numeric(kn[lk]), unname(gene_of[kt[lk]]), reorder = FALSE)
  data.frame(gene = rownames(s), count = as.integer(s[, 1]))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) %in% 2:3)
  rl <- if (length(args) == 3) as.integer(args[3]) else 75L
  res <- prepde_sample(args[1], rl)
  con <- gzfile(args[2], "w")
  utils::write.table(res, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  close(con)
}
