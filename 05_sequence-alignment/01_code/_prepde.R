## R port of StringTie's prepDE.py3 (github.com/gpertea/stringtie; the script has not changed
## since 2021-10-18, so this is the version shipped with StringTie 2.2.1, which made the
## counts in this repository). Base R only.
##
## prepde() builds a gene and a transcript count matrix from the per-sample GTFs written by
## `stringtie -e -B -G <reference>`. For every transcript in every sample
##
##     count = ceiling(cov * transcript_length / read_length)
##
## where `cov` is the transcript's coverage attribute and `transcript_length` the sum of its
## exon lengths. A gene's count is the sum of its transcripts' counts. Gene IDs are
## "gene_id|gene_name" when the transcript line carries a gene name (prepDE's unanchored
## regex also matches StringTie's `ref_gene_name`, hence IDs such as
## "gene-LOC134715138|LOC134715138"), else gene_id, else the transcript ID.
##
## prepDE.py (the Python 2 script StringTie also ships, and the one used for this project's
## matrices) computes the same counts, but Python 2 dictionaries have no defined order, so its
## rows come out in an arbitrary order.
##
## The output matches prepDE.py3 byte for byte: same samples (sorted by ID in plain
## character order), same row order, CRLF line ends and minimal quoting. Three deliberate
## differences, all in cases where prepDE.py3 writes a corrupt matrix instead of stopping:
##   - a duplicated sample ID stops (prepDE's duplicate check never fires);
##   - a transcript absent from the first sample stops (prepDE labels its gene "<class 'str'>");
##   - a GTF line with fewer than three fields stops with a message (prepDE raises IndexError).
## The -c option (cluster overlapping genes with StringTie-assigned IDs) is not ported: this
## pipeline quantifies against the reference only (`stringtie -e`), where it does not apply.
##
## Arguments mirror prepDE.py3: input (-i), gene_out (-g), transcript_out (-t),
## read_length (-l, default 75), pattern (-p, directory mode only). `eol` is the line ending:
## "\r\n" (default) as prepDE.py3 writes; "\n" for files committed to this repository, which
## stores text files with LF.

prepde <- function(input = ".", gene_out = "gene_count_matrix.csv",
                   transcript_out = "transcript_count_matrix.csv",
                   read_length = 75L, pattern = ".", verbose = FALSE, eol = "\r\n") {

  ## ---- sample list: a file of "<sample ID> <path to GTF>" lines, or a directory --------
  if (file.exists(input) && !dir.exists(input)) {
    lines <- readLines(input, warn = FALSE)
    lines <- lines[!startsWith(lines, "#")]
    fields <- lapply(trimws(lines), function(l) {
      f <- strsplit(l, "[[:space:]]+")[[1]]
      if (length(f) > 2) f <- c(f[1], f[2], paste(f[-(1:2)], collapse = " "))  # split(None, 2)
      f[nzchar(f)]
    })
    bad <- which(lengths(fields) != 2)
    if (length(bad))
      stop("line should have a sample ID and a file path:\n", lines[bad[1]], call. = FALSE)
    ids   <- vapply(fields, `[`, "", 1)
    paths <- vapply(fields, `[`, "", 2)
    if (anyDuplicated(ids))
      stop("non-unique sample ID (", ids[duplicated(ids)][1], ")", call. = FALSE)
    missing <- paths[!file.exists(paths)]
    if (length(missing)) stop("GTF file not found (", missing[1], ")", call. = FALSE)
  } else {
    if (!dir.exists(input)) stop("sub-directory '", input, "' not found!", call. = FALSE)
    subdirs <- list.dirs(input, full.names = FALSE, recursive = FALSE)
    subdirs <- subdirs[grepl(pattern, subdirs, perl = TRUE)]
    ids <- character(0); paths <- character(0)
    for (d in subdirs) {
      g <- Sys.glob(file.path(input, d, "*.gtf"))
      ids <- c(ids, rep(d, length(g))); paths <- c(paths, g)
    }
  }
  if (length(ids) == 0) stop("no GTF files found under base directory ", input, " !", call. = FALSE)
  o <- order(ids, paths, method = "radix")          # Python tuple sort: plain character order
  ids <- ids[o]; paths <- paths[o]

  rx_first <- function(rx, x) {                      # first (leftmost) match, as re.search
    m <- regmatches(x, regexec(rx, x, perl = TRUE))
    vapply(m, function(z) if (length(z)) z[2] else NA_character_, "")
  }
  gene_id_of <- function(attr, tid) {
    g <- rx_first('gene_id "([^"]+)"', attr)
    n <- rx_first('gene_name "([^"]+)"', attr)
    ifelse(is.na(g), tid, ifelse(is.na(n), g, paste0(g, "|", n)))
  }
  cov_of <- function(attr) {
    v <- suppressWarnings(as.numeric(rx_first('cov "([-+0-9.]+)"', attr)))
    v[is.na(v)] <- 0
    pmax(v, 0)
  }
  split_gtf <- function(path) {
    l <- readLines(path, warn = FALSE)
    list(lines = l, f = strsplit(l, "\t", fixed = TRUE))
  }

  ## ---- transcript -> gene, from the first sample only (as prepDE) ------------------------
  first <- split_gtf(paths[1])
  is_tx <- vapply(first$f, function(v) length(v) > 2 && v[3] == "transcript", NA)
  tx_attr <- vapply(first$f[is_tx], `[`, "", 9)
  tx_id <- rx_first('transcript_id "([^"]+)"', tx_attr)
  if (anyNA(tx_id)) stop("Problem parsing file ", paths[1], " at a transcript line", call. = FALSE)
  last_def <- !duplicated(tx_id, fromLast = TRUE)          # a repeated ID keeps its last gene (dict)
  gene_of <- setNames(gene_id_of(tx_attr, tx_id)[last_def], tx_id[last_def])

  ## ---- counts per sample ------------------------------------------------------------------
  ## `tx` grows in first-seen order; NA marks a sample with no entry for that transcript,
  ## which prepDE writes as an empty cell.
  tx_ids <- character(0)
  tx_val <- matrix(NA_integer_, 0, length(ids), dimnames = list(NULL, ids))
  gene_ids <- character(0)
  gene_val <- matrix(NA_real_, 0, length(ids), dimnames = list(NULL, ids))
  guides <- NULL

  for (k in seq_along(ids)) {
    s <- ids[k]
    if (verbose) message(">processing sample ", s, " from file ", paths[k])
    g <- if (k == 1) first else split_gtf(paths[k])

    ## header checks, on the first line only and only when it is a comment (as prepDE)
    if (length(g$lines) && startsWith(g$lines[1], "#")) {
      l <- g$lines[1]
      if (!grepl("-e", l, fixed = TRUE))
        stop("sample file ", paths[k], " was not generated with -e option!", call. = FALSE)
      gf <- rx_first("-G\\s*(\\S+)", l)
      if (is.na(gf)) stop("sample ", paths[k], " was not processed with -G option!", call. = FALSE)
      if (is.null(guides)) guides <- gf
      else if (gf != guides)
        message("Warning: sample file ", paths[k], " generated with a different -G file (",
                gf, ") than the first sample (", guides, ")")
    }

    ## Vectorised form of prepDE's line loop. Each transcript line opens a block; the exon
    ## lines after it add their lengths. A block is recorded when its exon length is > 0,
    ## and the last block is recorded whatever its length.
    body <- !startsWith(g$lines, "#")
    f <- g$f[body]
    short <- which(lengths(f) < 3)
    if (length(short))
      stop("GTF line ", which(body)[short[1]], " of ", paths[k], " has fewer than 3 fields", call. = FALSE)
    feat <- vapply(f, `[`, "", 3)
    is_t <- feat == "transcript"
    rec_t <- character(0); rec_n <- integer(0)
    if (any(is_t)) {
      block <- cumsum(is_t)
      if (any(feat == "exon" & block == 0))
        stop("exon line before the first transcript line in ", paths[k], call. = FALSE)
      ex <- feat == "exon"
      ex_len <- vapply(f[ex], function(v) as.integer(v[5]) - as.integer(v[4]) + 1L, 1L)
      tlen <- numeric(sum(is_t))
      if (any(ex)) { s_len <- rowsum(as.numeric(ex_len), block[ex]); tlen[as.integer(rownames(s_len))] <- s_len[, 1] }
      tf <- f[is_t]
      t_id <- rx_first('transcript_id "([^"]+)"', vapply(tf, function(v) v[length(v)], ""))
      cov  <- cov_of(vapply(tf, function(v) if (length(v) >= 9) v[9] else "", ""))
      keep <- tlen > 0 | seq_along(tlen) == length(tlen)
      rec_t <- t_id[keep]
      rec_n <- as.integer(ceiling(cov[keep] * tlen[keep] / read_length))
    }

    ## record this sample's transcript counts (a repeated ID keeps its last value, as a dict)
    new <- setdiff(unique(rec_t), tx_ids)
    if (length(new)) {
      tx_ids <- c(tx_ids, new)
      tx_val <- rbind(tx_val, matrix(NA_integer_, length(new), length(ids)))
    }
    last <- !duplicated(rec_t, fromLast = TRUE)
    tx_val[match(rec_t[last], tx_ids), s] <- rec_n[last]

    ## gene sums over every transcript seen so far; reading a missing entry gives 0 and
    ## stores it (Python defaultdict), so later empty cells only occur for earlier samples
    unknown <- setdiff(tx_ids, names(gene_of))
    if (length(unknown))
      stop("transcript ", unknown[1], " is not in the first sample (", ids[1], "); every sample ",
           "must be quantified against the same reference with `stringtie -e -G`", call. = FALSE)
    tx_val[is.na(tx_val[, s]), s] <- 0L
    gk <- unname(gene_of[tx_ids])
    newg <- setdiff(unique(gk), gene_ids)
    if (length(newg)) {
      gene_ids <- c(gene_ids, newg)
      gene_val <- rbind(gene_val, matrix(NA_real_, length(newg), length(ids)))
    }
    sums <- rowsum(as.numeric(tx_val[, s]), gk, reorder = FALSE)
    gene_val[match(rownames(sums), gene_ids), s] <- sums[, 1]
  }

  if (verbose) message("..writing ", transcript_out)
  prepde_write_csv("transcript_id", tx_ids, tx_val, transcript_out, eol)
  if (verbose) message("..writing ", gene_out)
  storage.mode(gene_val) <- "integer"
  prepde_write_csv("gene_id", gene_ids, gene_val, gene_out, eol)
  invisible(list(genes = gene_out, transcripts = transcript_out, samples = ids))
}

## Write a count matrix as prepDE.py3 does (Python csv.DictWriter: minimal quoting, CRLF line
## ends, empty cell for a missing value). `vals` has one column per sample, named by sample.
prepde_write_csv <- function(key, rows, vals, file, eol = "\r\n") {
  q <- function(x) {
    need <- grepl('[,"\r\n]', x)
    x[need] <- paste0('"', gsub('"', '""', x[need], fixed = TRUE), '"')
    x
  }
  v <- format(vals, scientific = FALSE, trim = TRUE)
  v[is.na(vals)] <- ""
  body <- paste(q(rows), apply(matrix(v, nrow(vals)), 1, paste, collapse = ","), sep = ",")
  con <- file(file, "wb"); on.exit(close(con))
  writeLines(c(paste(q(c(key, colnames(vals))), collapse = ","), body), con, sep = eol)
  invisible(file)
}

## prepDE's gene matrix from its transcript matrix: each gene's count is the sum of its
## transcripts' counts, and genes appear in the order their first transcript does, which is
## the order prepDE writes them. `tx_to_gene` is a named vector, transcript ID -> gene ID.
prepde_genes_from_transcripts <- function(transcript_csv, tx_to_gene, gene_out, eol = "\r\n") {
  tx <- utils::read.csv(transcript_csv, check.names = FALSE, colClasses = c(transcript_id = "character"))
  miss <- setdiff(tx$transcript_id, names(tx_to_gene))
  if (length(miss)) stop(length(miss), " transcript(s) have no gene ID, e.g. ", miss[1], call. = FALSE)
  m <- as.matrix(tx[, -1, drop = FALSE])
  g <- rowsum(m, unname(tx_to_gene[tx$transcript_id]), reorder = FALSE)
  storage.mode(g) <- "integer"
  prepde_write_csv("gene_id", rownames(g), g, gene_out, eol)
}
