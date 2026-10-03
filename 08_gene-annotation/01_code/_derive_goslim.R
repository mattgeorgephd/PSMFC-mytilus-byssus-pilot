## One-off provenance record: the generic GO slim of one GO release. Not part of the pipeline;
## it wrote 02_data/goslim_generic.obo, which step 01 reads.
##
## Why. Step 01 maps each DEG's GO IDs onto the generic slim through the GO graph in GO.db, so
## the slim has to be of the same GO release as GO.db: a slim term made obsolete in GO.db's
## release has no descendants there and silently counts no genes. The official file
## (https://current.geneontology.org/ontology/subsets/goslim_generic.obo) is only served for
## the current release, and the GO release archive was not reachable from the cloud environment
## the analysis ran in, so the slim is derived from the GO editors' file instead.
##
## How. The terms tagged `subset: goslim_generic` and not obsolete in
## geneontology/go-ontology src/ontology/go-edit.obo at the last commit before the release,
## with their id, name and namespace only (the GO graph comes from GO.db). Check that the commit
## is the release: its live GO terms should be GO.db's terms (GO.db adds the root "all"). For
## GO release 2026-01-23 (GO.db 3.23.1) the commit is d52c5a1 (2026-01-22 18:38 UTC), whose
## 38,739 live terms are exactly GO.db 3.23.1's; for release 2023-07-27 (GO.db 3.18.0) it was
## 17c29bb, and this script rebuilds that slim (02_data/_superseded/) byte for byte.
##
## Run from the repository root (needs raw.githubusercontent.com, or pass go-edit.obo's path):
##   Rscript 08_gene-annotation/01_code/_derive_goslim.R [path/to/go-edit.obo]

release  <- "2026-01-23"
commit   <- "d52c5a160a478c03c9f92cb1c67bf81695ff3f6b"
when     <- "2026-01-22 18:38 UTC, the last commit before the release"
godb     <- "GO.db 3.23.1"
out      <- file.path("08_gene-annotation", "02_data", "goslim_generic.obo")

derive_goslim <- function(obo, release, commit, when, godb, out) {
  x <- readLines(obo, encoding = "UTF-8")
  start <- grep("^\\[", x)                                  # stanza headers
  stanza <- cumsum(seq_along(x) %in% start)                 # 0 = file header
  terms <- lapply(split(x, stanza)[-1], function(s) {
    if (s[1] != "[Term]") return(NULL)
    tag <- function(k) sub(paste0("^", k, ": "), "", grep(paste0("^", k, ": "), s, value = TRUE))
    subsets <- sub(" .*$", "", tag("subset"))
    if (!"goslim_generic" %in% subsets || identical(tag("is_obsolete"), "true")) return(NULL)
    c(id = tag("id"), name = tag("name"), namespace = tag("namespace"))
  })
  terms <- do.call(rbind, terms[!vapply(terms, is.null, logical(1))])
  terms <- terms[order(terms[, "id"]), , drop = FALSE]
  header <- c("format-version: 1.2",
              paste0("data-version: go/releases/", release, " (goslim_generic subset, derived)"),
              "subsetdef: goslim_generic \"Generic GO slim\"",
              "ontology: go/subsets/goslim_generic",
              paste0("remark: Terms tagged 'subset: goslim_generic' (not obsolete) in geneontology/go-ontology ",
                     "src/ontology/go-edit.obo at commit ", commit, " (", when, "), the same GO release as ",
                     godb, ". Only id, name and namespace are kept; the GO graph comes from GO.db. The official ",
                     "release file is https://current.geneontology.org/ontology/subsets/goslim_generic.obo."))
  body <- unlist(lapply(seq_len(nrow(terms)), function(i)
    c("", "[Term]", paste0("id: ", terms[i, "id"]), paste0("name: ", terms[i, "name"]),
      paste0("namespace: ", terms[i, "namespace"]), "subset: goslim_generic")))
  writeLines(c(header, body), out, useBytes = TRUE)
  invisible(terms)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  obo <- if (length(args)) args[1] else {
    f <- tempfile(fileext = ".obo")
    download.file(paste0("https://raw.githubusercontent.com/geneontology/go-ontology/", commit,
                         "/src/ontology/go-edit.obo"), f, mode = "wb")
    f
  }
  terms <- derive_goslim(obo, release, commit, when, godb, out)
  cat(nrow(terms), "slim terms,", sum(terms[, "namespace"] == "biological_process"), "biological process, written to", out, "\n")
}
