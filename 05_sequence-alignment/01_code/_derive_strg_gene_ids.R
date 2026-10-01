## One-off provenance record: recover the gene IDs prepDE gave the 284 transcripts that the
## Ballgown table (03_analyses/hisat/t_data.ctab) lists without a gene. Not part of the
## pipeline; it wrote 02_data/strg_gene_ids.csv, which 02_prepDE.Rmd reads.
##
## Why it is needed. The per-sample StringTie GTFs that prepDE read are lost (they lived on a
## workstation). t_data.ctab still lists every reference transcript with its gene_id and
## gene_name, which give prepDE's gene ID ("gene_id|gene_name") for 67,819 of the 68,103
## transcripts. For the rest StringTie wrote its own gene ID (STRG.<n>, with ref_gene_name),
## which Ballgown records as ".". 481 of those keep their transcript ID as gene ID; the other
## 284 are named STRG.<n>|LOC<n> in the gene matrix.
##
## How. prepDE sums transcripts into genes, so each of those 284 genes' counts across all 129
## samples equal the summed counts of its transcripts. Every one of the 284 gene rows equals the
## count profile of exactly one of the 284 unassigned transcripts, with no ties and no all-zero
## profiles, so the assignment is one-to-one and unambiguous.
##
## It read the original clean gene matrix, which then sat at
## 06_differential-expression/02_data/gene_count_matrix_clean.csv; that file is now rebuilt by
## the pipeline (06_differential-expression/03_analyses/count_matrix/), partly from this
## script's output, so the original is taken from the repository's first commit.
##
## Run from the repository root:  Rscript 05_sequence-alignment/01_code/_derive_strg_gene_ids.R

td <- read.delim("05_sequence-alignment/03_analyses/hisat/t_data.ctab")
tx <- read.csv("05_sequence-alignment/03_analyses/prepDE/transcript_count_matrix.csv", check.names = FALSE)
gc <- read.csv(pipe("git show 3bee414:06_differential-expression/02_data/gene_count_matrix_clean.csv"),
               row.names = 1, check.names = FALSE)

gid <- ifelse(td$gene_name != ".", paste0(td$gene_id, "|", td$gene_name), td$gene_id)
no_gene <- td$gene_id == "."
gid[no_gene] <- td$t_name[no_gene]                          # prepDE's fallback: the transcript ID
unassigned <- td$t_name[no_gene & !(gid %in% rownames(gc))]
strg_genes <- setdiff(rownames(gc), unique(gid))
stopifnot(length(unassigned) == length(strg_genes))

## sample columns as in the cleaned gene matrix (plate well and _S<n> stripped, T051 removed)
m <- as.matrix(tx[, -1]); rownames(m) <- tx$transcript_id
colnames(m) <- sub("^[A-H][0-9]+-", "", sub("_S[0-9]+$", "", colnames(m)))
m <- m[unassigned, colnames(gc), drop = FALSE]
profile_t <- apply(m, 1, paste, collapse = ",")
profile_g <- apply(as.matrix(gc[strg_genes, ]), 1, paste, collapse = ",")
stopifnot(!anyDuplicated(profile_t), all(rowSums(m) > 0))

hit <- match(profile_g, profile_t)
stopifnot(!anyNA(hit), !anyDuplicated(hit))                 # one-to-one, every gene matched
out <- data.frame(transcript_id = unassigned[hit], gene_id = strg_genes)
out <- out[order(out$transcript_id, method = "radix"), ]
write.csv(out, "05_sequence-alignment/02_data/strg_gene_ids.csv", row.names = FALSE)
message("Wrote 05_sequence-alignment/02_data/strg_gene_ids.csv (", nrow(out), " transcripts).")
