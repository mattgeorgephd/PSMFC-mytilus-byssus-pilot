# 02_data

Large external inputs, downloaded by the scripts into this folder and git-ignored:

| File | Source | Downloaded by |
|---|---|---|
| `Mtros-hq_transcripts.fasta` | owl, `halfshell/genomic-databank/` (411,251 high-quality Iso-Seq isoforms, `MT_Pool_HQ_transcript/N`; 852 MB; MD5 e71381f7b6970fd20ff6b0e72006c80a) | `01_isoseq_transcriptome_check.Rmd`, `02_isoform_gene_map.Rmd` |
| `cds_from_genomic.fasta` | the CDS of the genome annotation GCF_036588685.1 (53,269 CDS, 34,235 genes); gannet copy at `panopea/PSMFC-mytilus-byssus-pilot/byssus-exp-analysis/data/ncbi_dataset/data/GCF_036588685.1/` | `02_isoform_gene_map.Rmd` |

The trimmed Tag-seq reads are not stored: `03_salmon_quant.Rmd` streams each of the 131 files
from gannet (`.../byssus-exp-analysis/data/raw-trimmed/`), quantifies it and deletes it.
