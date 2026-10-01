# 02_data

| Item | Description | Read by |
|------|-------------|---------|
| `sample-submission/` | Tag-seq submission paperwork: sample lists, RNA summary, sequencing quote, RIN trace, shipping instructions, upload template | (record) |
| `annotation_mt_like_loci.csv` | the 293 genes and pseudogenes of the RefSeq annotation (GCF_036588685.1, release RS_2024_02) that NCBI names after a mitochondrially encoded protein (126 protein-coding copies, 167 pseudogene copies, all on unplaced scaffolds), with the protein each copies; derived once from the annotation's GFF by `../01_code/_derive_mt_like_loci.R` | `tools/mt_encoded.R` via `06` step 01, which leaves them out of the genome analysis with the mitogenome genes |
| `strg_gene_ids.csv` | the StringTie gene ID of the 284 reference transcripts whose gene was named by StringTie (`STRG.854|LOC134681509`), recovered once from the original gene matrix by `../01_code/_derive_strg_gene_ids.R` | `02_prepDE.Rmd` (when it rebuilds the gene matrix from the transcript matrix) |

Sequencing reads are not stored here. The raw (untrimmed) reads of all 131 libraries live on
owl at https://owl.fish.washington.edu/nightingales/M_trossulus/ (gannet
`panopea/PSMFC-mytilus-byssus-pilot/20220405-tagseq/` has 73 of them); the trimmed reads of all 131, which `01_hisat_stringtie.Rmd` aligned, at
https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/byssus-exp-analysis/data/raw-trimmed/.
How they were trimmed is in `../README.md` ("The reads") and `../01_code/03_read_trimming.Rmd`.
The genome (`GCF_036588685.1`) is downloaded from NCBI by the alignment script into
`ncbi_dataset/` here (git-ignored).
