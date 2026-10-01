## Nuclear-genome gene models annotated as mitochondrially encoded proteins. Base R.
##
## About 140 LOCs in the M. trossulus genome annotation have as best UniProt hit one of the 13
## proteins encoded by mtDNA (COX1-3, CYTB, ND1-6, ND4L, ATP6, ATP8). Reads from the abundant
## mitochondrial transcripts spread over them, so one mitochondrial signal can show up as
## dozens of DEGs and dominate GO results: in Gill OA, 116 of the 711 TC DEGs are such LOCs,
## all up with log2 fold changes near 0.5. Scripts use this to flag them, not to drop them.
MT_ENCODED_RX <- paste0("^(Cytochrome c oxidase subunit [123] |Cytochrome b \\(|",
                        "NADH-ubiquinone oxidoreductase chain (1|2|3|4|4L|5|6) |",
                        "ATP synthase subunit a \\(|ATP synthase protein 8)")
is_mt_encoded <- function(protein_name) !is.na(protein_name) & grepl(MT_ENCODED_RX, protein_name)
