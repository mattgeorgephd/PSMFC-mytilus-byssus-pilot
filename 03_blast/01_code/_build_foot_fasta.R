## One-off provenance record: the Mytilus foot and byssal proteins of the 2026 genome BLAST.
## Not part of the pipeline; it wrote 02_data/uniprotkb_mytilus_foot_2026_03_byssal.fasta, which
## 01_genome_blast.Rmd searches with Swiss-Prot 2026_03 (foot_fasta).
##
## Why. The 2024 search added to Swiss-Prot the proteins of the UniProt free-text query
## "(mytilus foot)", which finds entries that mention both words. On 2026-10-03 (release 2026_03)
## it returned 196 proteins (02_data/uniprotkb_mytilus_foot_2026_03.fasta) and missed byssal
## proteins whose records do not say "foot": the thread matrix proteins, proximal thread matrix
## protein 1, the precollagens of M. galloprovincialis and M. californianus, mefp-5, the
## M. coruscus byssus proteome and others. A search of UniProt (Mytilidae entries naming a byssal
## protein family), NCBI Protein and the primary literature found them; the 45 confirmed in the
## literature, all sequenced from transcripts or protein and none a genome gene-model prediction,
## are listed with their papers in 02_data/byssal_additions_2026_03.tsv. They are added here.
##
## How. The query's proteins (in UniProt's order) and then the 45 (in the table's order), as
## UniProt's REST service returns them in FASTA; it stops unless UniProt serves release 2026_03,
## and checks that the query part is the 196 of uniprotkb_mytilus_foot_2026_03.fasta.
##
## Run from the repository root (needs rest.uniprot.org):
##   Rscript 03_blast/01_code/_build_foot_fasta.R

release   <- "2026_03"
dat       <- file.path("03_blast", "02_data")
additions <- read.delim(file.path(dat, "byssal_additions_2026_03.tsv"), stringsAsFactors = FALSE)
out       <- file.path(dat, "uniprotkb_mytilus_foot_2026_03_byssal.fasta")

served <- curlGetHeaders("https://rest.uniprot.org/uniprotkb/search?query=accession:P69905&size=1")
served <- trimws(sub("^[^:]*:", "", grep("^x-uniprot-release:", served, ignore.case = TRUE, value = TRUE)[1]))
if (!identical(served, release)) stop("UniProt serves release ", served, ", not ", release)

query <- readLines("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28mytilus+foot%29", warn = FALSE)
stopifnot(identical(query, readLines(file.path(dat, "uniprotkb_mytilus_foot_2026_03.fasta"))))
added <- readLines(paste0("https://rest.uniprot.org/uniprotkb/accessions?format=fasta&accessions=",
                          paste(additions$accession, collapse = ",")), warn = FALSE)
acc_of <- function(fa) sub("^>[^|]*[|]([^|]*)[|].*$", "\\1", grep("^>", fa, value = TRUE))
stopifnot(setequal(acc_of(added), additions$accession), !any(acc_of(added) %in% acc_of(query)))
## the additions in the table's order
starts <- grep("^>", added); blocks <- split(added, cumsum(seq_along(added) %in% starts))
names(blocks) <- acc_of(added)
writeLines(c(query, unlist(blocks[additions$accession], use.names = FALSE)), out)
cat(length(acc_of(query)), "query proteins +", nrow(additions), "byssal additions written to", out, "\n")
