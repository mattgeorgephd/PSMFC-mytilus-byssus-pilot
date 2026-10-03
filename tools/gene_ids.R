## Gene names in the count matrix, and the key that joins them to annotation. Base R.
##
## The count matrix of record (04_sequence-alignment step 07: featureCounts, which names a gene
## by its key) carries the names prepDE gives a gene: "gene_id|gene_name" when the reference
## has a gene name, else "gene_id", else the transcript name. Four forms:
##   gene-LOC134721619|LOC134721619   reference gene (41,271 rows)
##   gene-LOC134700001                reference gene without a gene_name attribute (739)
##   gene-COX1|COX1                   mitochondrial genes, and tRNA genes such as
##   gene-Trnaa-agc-10|Trnaa-agc      (the tRNA gene_name is shared by up to ~200 genes)
##   rna-NC_007687.1:10108..10170     mitochondrial tRNAs and rRNAs without a gene record (26)
## The previous matrix (StringTie + prepDE, 04_sequence-alignment/03_analyses/prepDE/) also
## had genes StringTie filed under its own IDs: "STRG.10|LOC134702910" (263) and "STRG.12" (21,
## the mitochondrial tRNAs and rRNAs that are named by transcript above).
##
## gene_key() returns the LOC identifier when the name holds one, otherwise the gene_id
## without its "gene-" prefix. That is the key of the BLAST / UniProt table
## (03_blast/03_analyses/genome-foot-uniprot2026_03/LOC_GO_list.txt, column `gene`: LOC IDs, and ND2, CYTB,
## ... for the mitochondrial genes) and of the reference annotation (t_data.ctab gene_id
## without "gene-"). It is unique across the matrix; taking the text after "|" is not, since
## it collapses the tRNA genes onto their shared gene_name and leaves "gene-" on the 739
## names without "|".
gene_key <- function(gene) {
  gene <- as.character(gene)
  has_loc <- grepl("LOC[0-9]+", gene)
  ifelse(has_loc, sub("^.*?(LOC[0-9]+).*$", "\\1", gene, perl = TRUE),
         sub("^gene-", "", sub("[|].*$", "", gene)))
}
