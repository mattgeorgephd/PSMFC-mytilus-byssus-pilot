# 02_data

| Item | Description | Read by |
|------|-------------|---------|
| `sample-submission/` | Tag-seq submission paperwork: sample lists, RNA summary, sequencing quote, RIN trace, shipping instructions, upload template | (record) |
| `strg_gene_ids.csv` | the StringTie gene ID of the 284 reference transcripts whose gene was named by StringTie (`STRG.854|LOC134681509`), recovered once from the original gene matrix by `../01_code/_derive_strg_gene_ids.R` | `02_prepDE.Rmd` (when it rebuilds the gene matrix from the transcript matrix) |

Sequencing reads are not stored here. The raw (untrimmed) reads of 73 of the 131 libraries
live on gannet at https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/20220405-tagseq/;
the trimmed reads of all 131, which `01_hisat_stringtie.Rmd` aligned, at
https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/byssus-exp-analysis/data/raw-trimmed/.
How they were trimmed is in `../README.md` ("The reads") and `../01_code/03_read_trimming.Rmd`.
The genome (`GCF_036588685.1`) is downloaded from NCBI by the alignment script into
`ncbi_dataset/` here (git-ignored).
