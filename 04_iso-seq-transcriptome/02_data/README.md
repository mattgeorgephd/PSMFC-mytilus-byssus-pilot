# 02_data

Large external inputs, downloaded by the scripts into this folder and git-ignored:

| File | Source | Downloaded by |
|---|---|---|
| `Mtros-hq_transcripts.fasta` | owl, `halfshell/genomic-databank/` (411,251 high-quality Iso-Seq isoforms, `MT_Pool_HQ_transcript/N`; 852 MB; MD5 e71381f7b6970fd20ff6b0e72006c80a) | `01_isoseq_transcriptome_check.Rmd`, `02_isoform_gene_map.Rmd` |
| `GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.fna.gz`, `..._genomic.gff.gz` | the genome (1.30 Gb, 541 sequences) and its RefSeq annotation (release RS_2024_02), from NCBI `genomes/all/GCF/036/588/685/GCF_036588685.1_PNRI_Mtr1.1.1.hap1/`; MD5s checked against NCBI's `md5checksums.txt` | `02_isoform_gene_map.Rmd` |
| `cds_from_genomic.fasta` | the CDS of the same annotation (53,269 CDS, 34,235 genes); gannet copy at `panopea/PSMFC-mytilus-byssus-pilot/byssus-exp-analysis/data/ncbi_dataset/data/GCF_036588685.1/` | `_superseded/02_isoform_gene_map_cds.Rmd` (retired) |

The trimmed Tag-seq reads are not stored: `03_salmon_quant.Rmd` streams each of the 131 files
from gannet (`.../byssus-exp-analysis/data/raw-trimmed/`), quantifies it and deletes it.
