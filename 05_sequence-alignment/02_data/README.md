# 02_data

| Item | Description |
|------|-------------|
| `sample-submission/` | Tag-seq submission paperwork: sample lists, RNA summary, sequencing quote, RIN trace, shipping instructions, upload template |

Sequencing reads are not stored here. The raw (untrimmed) reads live on gannet at
https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/20220405-tagseq/; the
trimmed reads `13-Hisat.Rmd` reads were produced on the HPC (that folder holds raw reads
only, and the trimming step is not in this repository). The genome (`GCF_036588685.1`) is
downloaded from NCBI by the alignment script.
