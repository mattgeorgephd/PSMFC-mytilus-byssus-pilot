# 02_data

| Item | Description | Read by |
|------|-------------|---------|
| `sample-submission/` | Tag-seq submission paperwork: sample lists, RNA summary, sequencing quote, RIN trace, shipping instructions, upload template | (record) |
| `annotation_mt_like_loci.csv` | the 293 genes and pseudogenes of the RefSeq annotation (GCF_036588685.1, release RS_2024_02) that NCBI names after a mitochondrially encoded protein (126 protein-coding copies, 167 pseudogene copies, all on unplaced scaffolds), with the protein each copies; derived once from the annotation's GFF by `../01_code/_derive_mt_like_loci.R` | `tools/mt_encoded.R`, via steps 04 and 05 (which keep the annotation's mitochondrial loci as they are) and `05` step 01 (which leaves them out of the genome analysis with the mitogenome genes) |
| `mitogenome_NC_007687.1.fa` | the female-type mitochondrial genome of *M. trossulus* (NC_007687.1, 18,652 bases), cut once from the RefSeq assembly GCF_036588685.1 by `../01_code/_derive_mitogenome.R` | `08_mitogenome_counts.Rmd` (its HISAT2 index) |
| `mitogenome_genes.saf` | its 38 annotated features (annotation release RS_2024_02): the 12 protein genes by symbol (ATP6, COX1-3, CYTB, ND1-6, ND4L), the 2 rRNAs and 24 tRNAs by RNA ID (`rna-NC_007687.1:1..1244`), as featureCounts SAF; derived with the FASTA | `08_mitogenome_counts.Rmd` |
| `strg_gene_ids.csv` | the StringTie gene ID of the 284 reference transcripts whose gene was named by StringTie (`STRG.854|LOC134681509`), recovered once from the original gene matrix by `../01_code/_derive_strg_gene_ids.R` | `02_prepDE.Rmd` (when it rebuilds the gene matrix from the transcript matrix) |

Sequencing reads are not stored here. The raw (untrimmed) reads of all 131 libraries live on
owl at https://owl.fish.washington.edu/nightingales/M_trossulus/ (gannet
`panopea/PSMFC-mytilus-byssus-pilot/20220405-tagseq/` has 73 of them); the trimmed reads of all 131, which `01_hisat_stringtie.Rmd` aligned, at
https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/byssus-exp-analysis/data/raw-trimmed/.
How they were trimmed is in `../README.md` ("The reads") and `../01_code/03_read_trimming.Rmd`.
The genome (`GCF_036588685.1`) is downloaded from NCBI by the alignment script (step 01, HPC)
into `ncbi_dataset/` here (git-ignored).

Large external inputs of steps 04 to 06, downloaded by the scripts into this folder and
git-ignored (until 2026-10-03 in `06_iso-seq-transcriptome/02_data/`, with those steps):

| File | Source | Downloaded by |
|---|---|---|
| `Mtros-hq_transcripts.fasta` | owl, `halfshell/genomic-databank/` (411,251 high-quality Iso-Seq isoforms, `MT_Pool_HQ_transcript/N`; 852 MB; MD5 e71381f7b6970fd20ff6b0e72006c80a) | `04_isoform_gene_map.Rmd` (also read by `06_iso-seq-transcriptome` step 02 for its salmon index) |
| `GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.fna.gz`, `..._genomic.gff.gz` | the genome (1.30 Gb, 541 sequences) and its RefSeq annotation (release RS_2024_02), from NCBI `genomes/all/GCF/036/588/685/GCF_036588685.1_PNRI_Mtr1.1.1.hap1/`; MD5s checked against NCBI's `md5checksums.txt` | `04_isoform_gene_map.Rmd` |
| `GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.gtf.gz` | the same annotation as GTF, from the same NCBI folder (MD5 checked); its splice sites, extracted with HISAT2's `hisat2_extract_splice_sites.py` as step 01 did, guide the realignment | `06_genome_recount.Rmd` |
| `cds_from_genomic.fasta` | the CDS of the same annotation (53,269 CDS, 34,235 genes); gannet copy at `panopea/PSMFC-mytilus-byssus-pilot/byssus-exp-analysis/data/ncbi_dataset/data/GCF_036588685.1/` | `_superseded/02_isoform_gene_map_cds.Rmd` (retired) |

The trimmed reads are not stored: `06_genome_recount.Rmd` and `08_mitogenome_counts.Rmd`
(with `online: true`) download each of the 131 files from gannet, process it and delete it.
