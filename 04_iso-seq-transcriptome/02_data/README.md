# 02_data

Large external inputs, downloaded by the scripts into this folder and git-ignored:

| File | Source | Downloaded by |
|---|---|---|
| `Mtros-hq_transcripts.fasta` | owl, `halfshell/genomic-databank/` (411,251 high-quality Iso-Seq isoforms, `MT_Pool_HQ_transcript/N`; 852 MB; MD5 e71381f7b6970fd20ff6b0e72006c80a) | `01_isoseq_transcriptome_check.Rmd`; `02_salmon_quant.Rmd` uses `05_sequence-alignment/02_data/`'s copy (downloaded by that folder's step 04) when it is there, else this one |

The genome, its annotation and the annotation's CDS, which the isoform-to-gene map, the
augmented annotations and the genome recount download, are in `05_sequence-alignment/02_data/`
since those steps moved there (2026-10-03). The trimmed Tag-seq reads are not stored:
`02_salmon_quant.Rmd` downloads each of the 131 files from gannet
(`.../byssus-exp-analysis/data/raw-trimmed/`), processes it and deletes it.
