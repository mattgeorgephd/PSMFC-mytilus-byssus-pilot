# 01_code/_superseded

Earlier read-processing and alignment scripts, kept as records. None is run by the pipeline,
and none produced the reads or counts the analysis uses.

| Script | What it did | Replaced by |
|---|---|---|
| `1_1_process-tagseq-data-mytilus.Rmd` | Matt George's 2022 processing (copied unchanged on 2026-10-01 from gannet `seashell/bu-github/PSMFC-mytilus-byssus-pilot/`): downloads the raw reads of the first 73 libraries, FastQC, trims each lane with cutadapt (`-a A{8} -a G{8} -a AGATCGG -u 15 -m 20`; no leader clipping, no duplicate removal), concatenates the two lanes, then HISAT2 against the *M. edulis* (GCA_019925275.1) and *M. galloprovincialis* (GCA_900618805.1) genomes and a *M. trossulus* Trinity assembly (SRR6051615), with mean alignment rates of 34 to 35% to the two genomes | the `*_L099_R1_cmb.trim.fastq.gz` reads, made with a different recipe (`../03_read_trimming.Rmd`), aligned to GCF_036588685.1 by `../01_hisat_stringtie.Rmd`. `../../03_analyses/fastqc/multiqc_report_trimmed_merged_2022-08-09.html` is the QC of this script's trimmed reads |
| `07-HiSat_GL.Rmd` | an earlier HISAT2 attempt against a different assembly (GCA_032158295.1) with an augustus annotation (`mod_augustus.gtf`) | `../01_hisat_stringtie.Rmd` |
| `07-kallisto.Rmd`, `07-kallisto.md` | kallisto pseudo-alignment of the trimmed reads to the Iso-Seq transcriptome, with a fixed fragment length (`-l 65 -s 2`) and so a length correction that 3' Tag-seq does not need | `../01_hisat_stringtie.Rmd` (genome) and `04_iso-seq-transcriptome` (salmon, no length correction) |
| `02_isoform_gene_map_cds.Rmd` | the first isoform-to-gene map (2026-10-01, in `04_iso-seq-transcriptome` until 2026-10-03; paths inside are that folder's), built while the genome was not reachable: each isoform aligned to the annotation's CDS (minimap2 `-x asm20`), identity at least 0.90 over at least 100 bases. It reaches coding genes only and misses isoforms that end before the CDS. Output: `../../03_analyses/_superseded/02_isoform-gene-map_cds/` | `../04_isoform_gene_map.Rmd` (spliced alignment to the genome and exon overlap with every annotated gene), which compares its map with this one |
