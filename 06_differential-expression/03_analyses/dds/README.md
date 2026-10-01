# dds

One `<code>.rds` per contrast: the fitted DESeq2 object after the count filter (at least 10
counts in at least a third of the contrast's samples), written by
`../../01_code/03_deseq_contrasts.Rmd` and read by `04_shrinkage_filtration.Rmd`. The `.rds`
files are git-ignored and rebuilt on every run (about six minutes).
