# 05_sequence-alignment

Read QC, alignment and quantification of the Tag-seq reads against the *Mytilus trossulus*
genome (GenBank GCA_036588685.1 / RefSeq GCF_036588685.1), ending in the count matrix that
`06_differential-expression` starts from.

**The count matrix of record** (since 2026-10-02) is `03_analyses/featurecounts/gene_count_matrix.csv`,
written by step 07: the trimmed reads realigned with HISAT2 and counted with featureCounts
(uniquely aligned reads, sense strand) on the RefSeq annotation with its 3' ends extended by
the Iso-Seq isoforms. Step 04 places the Iso-Seq isoforms on the genome, step 05 builds the
extended annotation (`ext3`, option B) and step 06 realigns and counts the reads; those three
download their inputs and run only with `online: true`, so the pipeline's step 07 reads step
06's committed matrix. Until 2026-10-03 steps 04 to 06 were steps 02, 05 and 06 of
`04_iso-seq-transcriptome`, and steps 07 and 08 here were 04 and 05. Why the matrix replaced
the previous one is in `04_iso-seq-transcriptome/03_analyses/04_augmented-de/README.md`:
featureCounts counts reads (prepDE estimates read bases over 75 while the reads average about
63 bases, and its gene-wise dispersions ran about 9% higher), and the RefSeq gene models often
end before the 3' ends Tag-seq reads.

**The previous matrix** came from HISAT2 + StringTie (`01_code/01_hisat_stringtie.Rmd`), run
on a collaborator HPC workstation (paths under `/home/shared/...`) against trimmed reads and a
genome that are not stored in this repository; it is kept verbatim as the method record, with
its summary outputs committed in `03_analyses/`. Step 01 is a break in the pipeline: its
committed outputs are what later steps read. Its count-matrix step (`02_prepDE.Rmd`) runs
anywhere, and its matrix is still written, for `04_iso-seq-transcriptome` step 04's comparison.

## How to run

Open `sequence-alignment.Rproj` and knit `01_code/00_run_sequence_alignment.Rmd` (or let the
repository-level `00_run_pipeline.Rmd` do it). By default it runs steps 02, 03, 07 and 08 from
committed files (step 08 only summarises its committed counts). With `online: true` it also
runs steps 04 to 06 and step 08's alignment, which download the Iso-Seq transcriptome (owl),
the genome and its annotation (NCBI) and the trimmed reads (gannet), and need minimap2 2.31,
HISAT2 2.2.1, StringTie 2.2.1, Subread 2.1.1 and samtools (the runner's `minimap2` and `bin`
parameters; bioconda builds were used). Their large intermediate files (alignments,
annotations, the HISAT2 index, per-library counts) are git-ignored and reused when present.

| step | script | what it does |
|---|---|---|
| 01 | `01_hisat_stringtie.Rmd` | HISAT2 alignment and StringTie quantification (`-e -B` against the reference annotation); writes `03_analyses/hisat/` and its `sample_list.txt`. HPC only |
| 02 | `02_prepDE.Rmd` | gene and transcript count matrices in `03_analyses/prepDE/` (read length 75): the previous record |
| 03 | `03_read_trimming.Rmd` | how the trimmed reads were made: per-library read retention from the committed FastQC tables, and (with the script's own `online: true`, knit by hand) the trimming recipe reproduced on a sample of reads; writes `03_analyses/read_trimming/` |
| 04 | `04_isoform_gene_map.Rmd` | each Iso-Seq isoform placed on the genome (minimap2) and given the annotated gene whose exons it overlaps most, or a novel locus, a mitochondrial feature or none; writes `03_analyses/isoform-gene-map/`. Online; about 1.5 hours and 9 GB of memory |
| 05 | `05_augmented_annotation.Rmd` | option B's annotations: RefSeq unchanged (the control), `ext3` (3' ends extended) and `full` (isoforms and novel loci added); writes `03_analyses/augmented-annotation/`. Online (step 04's alignment); a few minutes |
| 06 | `06_genome_recount.Rmd` | the 131 libraries realigned with HISAT2 and counted on each annotation with StringTie + prepDE and featureCounts (six gene count matrices); writes `03_analyses/genome-recount/`. Online; about 15 minutes for the index plus 4 minutes per library |
| 07 | `07_count_matrix_of_record.Rmd` | the count matrix of record in `03_analyses/featurecounts/`: step 06's featureCounts matrix on `ext3`, its rows named as prepDE names them (`gene_id|gene_name` from `03_analyses/hisat/t_data.ctab`) |
| 08 | `08_mitogenome_counts.Rmd` | the mitochondrial genes counted on the mitochondrial genome alone, in `03_analyses/mitogenome/`: with `online: true` each library's reads are downloaded and aligned to NC_007687.1 with HISAT2 (default scoring, of record, and a permissive score) and counted per gene with featureCounts (`_mitogenome_library.sh`); read by `06_differential-expression` step 13. About 1 hour; offline, the committed counts are summarised |

## Design of the counts of record (steps 04 to 06)

1. **Isoforms to genes (step 04).** The Iso-Seq transcriptome's 411,251 isoforms were
   clustered but never collapsed to genes. Each is aligned to the genome (RefSeq
   GCF_036588685.1) with minimap2's spliced preset for accurate full-length transcripts
   (`-x splice:hq -uf`; Li 2018, *Bioinformatics* 34:3094) and takes the annotated gene
   (annotation release RS_2024_02) whose exons it overlaps most on its strand. Isoforms aligned
   well but outside every gene form novel loci; unaligned ones stay single features. Isoforms
   on mitochondrial loci are flagged; the annotation's mitochondrial loci (the mitogenome's
   genes and their 293 nuclear copies, `tools/mt_encoded.R`) are found here from the
   annotation itself, the same 331 loci `06` lists.
2. **Annotations (step 05).** Three versions of RS_2024_02: `refseq` unchanged (the control);
   `ext3`, each transcript's last exon extended to the 3' end of same-gene isoforms whose last
   block overlaps it, stopping before the next gene on the strand; and `full`, the isoforms
   added as transcripts of their genes and the novel loci as new genes (a locus overlapping
   one gene's exons joins it). Readthrough isoforms (exons on two genes), ambiguous
   placements and mitochondrial isoforms are not used, and the mitochondrial loci are kept
   unchanged.
3. **Realignment and counting (step 06).** Step 01's alignments were not kept, so the 131
   libraries are realigned with HISAT2 2.2.1 and step 01's settings. Step 01's index held the
   RefSeq splice sites and exons, which takes more memory than the 15 GB machine used here;
   here the same splice sites are given at alignment time. Each alignment is counted on the
   three annotations with StringTie + prepDE, as the previous record counted (Pertea et al.
   2015, *Nature Biotechnology* 33:290), and with featureCounts (unique alignments, sense
   strand; Liao, Smyth and Shi 2014, *Bioinformatics* 30:923), the usual counter for 3'
   Tag-seq. `04_iso-seq-transcriptome` step 04 compares the six matrices; featureCounts on
   `ext3` became the record (step 07).

`02_prepDE.Rmd` counts from the per-sample StringTie GTFs when step 01's outputs are present,
with `01_code/_prepde.R`, an R port of StringTie's `prepDE.py3` checked to give byte-identical
matrices. Without them, which is the case in this repository, it rebuilds the gene matrix from
the committed transcript matrix (the original `prepDE.py` output) through the reference
transcript table (`03_analyses/hisat/t_data.ctab`) and `02_data/strg_gene_ids.csv`, and checks
the sums. `01_code/_derive_strg_gene_ids.R` is the one-off record of how that ID map was
recovered.

## The reads

GSAF (University of Texas at Austin, job JA22078) sequenced the 131 libraries as 3' Tag-seq,
single-end 100 bp on a NovaSeq S1, over two lanes. The reads aligned in step 01 are the files
`*_L099_R1_cmb.trim.fastq.gz`: one per library, the two lanes combined (`cmb`, lane `L099`),
trimmed. Who trimmed them was not recorded, but step 03 reproduces the recipe: 99.8% of the
trimmed reads in a sample of library T001F (32,296 of 32,349) are identical, base for base and
in quality, to the output of

    tagseq_clipper.pl raw.fastq \
      | fastx_clipper -a AAAAAAAA -l 20 -Q33 \
      | fastx_clipper -a AGATCGGAAG -l 20 -Q33

that is, the Matz-lab Tag-seq clipper (https://github.com/z0on/tag-based_RNAseq; Meyer,
Aglyamova and Matz 2011, *Molecular Ecology* 20:3599), which cuts the 5' leader (degenerate
bases and the template-switching G run) and drops PCR duplicates (reads sharing the leader and
the first 20 bases), then FASTX-Toolkit 0.0.14 `fastx_clipper` for the poly-A tail and the
Illumina adapter, keeping reads of 20 bases or more, with no quality filter. The other recipes
tried match far fewer reads (`03_analyses/read_trimming/recipe_check.csv`). As a result the
analysed reads are deduplicated: in the 72 libraries with both raw lanes in the FastQC tables
(`03_analyses/fastqc/untrimmed/`), a median of 43% of raw reads remain (35 to 58%), and the 131
trimmed libraries hold 0.78 to 3.86 million reads (median 2.72 million) of 20 to 95 bases (mean
about 63).

For anyone rerunning from raw reads:

- **Where the raw reads are.** All 131 libraries, both lanes (262 files, the April 2022
  delivery), are on owl at https://owl.fish.washington.edu/nightingales/M_trossulus/
  (`T*_S*_L00[12]_R1_001.fastq.gz`; the paired `69M_1.fastq.gz`-style files there are not
  from this Tag-seq run). gannet `20220405-tagseq/` holds a partial copy: 145 files of 73 libraries (T001 to
  T030, T131G and T132 to T137; T030G lane 1 only), the set the FastQC tables cover. The one
  file compared (T030G lane 1) is byte-identical on both servers; owl keeps no checksums for
  this set.
- **The August 2022 trimming is a different one.** Matt George's 2022 script
  (`01_code/_superseded/1_1_process-tagseq-data-mytilus.Rmd`) trimmed the first 73 libraries with
  cutadapt (a fixed 15-base 5' cut, no deduplication) and aligned them to other *Mytilus*
  genomes; `03_analyses/fastqc/multiqc_report_trimmed_merged_2022-08-09.html` is the QC of those
  reads, not of the reads analysed here.

## Layout

```
05_sequence-alignment/
├── sequence-alignment.Rproj
├── 01_code/
│   ├── 00_run_sequence_alignment.Rmd   batch runner
│   ├── 01_hisat_stringtie.Rmd          HISAT2 + StringTie (HPC record)
│   ├── 02_prepDE.Rmd                   count matrices
│   ├── 03_read_trimming.Rmd            how the trimmed reads were made (retention, recipe check)
│   ├── 04_isoform_gene_map.Rmd         Iso-Seq isoforms to genome genes (online)
│   ├── 05_augmented_annotation.Rmd     option B: RefSeq augmented with the isoforms (online)
│   ├── 06_genome_recount.Rmd           option B: realign and count on each annotation (online)
│   ├── 07_count_matrix_of_record.Rmd   the count matrix of record (featureCounts on ext3, from 06)
│   ├── 08_mitogenome_counts.Rmd        mitochondrial reads on the mitogenome alone
│   ├── _paths.R                        shared paths of steps 04 to 06
│   ├── _recount_library.sh             one library of step 06 (align, count, clean up)
│   ├── _prepde_sample.R                prepDE's gene counts for one StringTie sample (step 06)
│   ├── _mitogenome_library.sh          one library of step 08 (alignment and counts)
│   ├── _derive_mitogenome.R            one-off: the mitogenome FASTA and its genes
│   ├── _prepde.R                       R port of prepDE.py3 (sourced by 02)
│   ├── _derive_strg_gene_ids.R         one-off: recovers the StringTie gene IDs
│   ├── _derive_mt_like_loci.R          one-off: loci the annotation names after a mitochondrial protein
│   └── _superseded/                    records (README inside)
│       ├── 1_1_process-tagseq-data-mytilus.Rmd   2022 cutadapt trimming + HISAT2 to other genomes
│       ├── 07-HiSat_GL.Rmd             earlier HISAT2 attempt (different assembly + augustus)
│       ├── 07-kallisto.Rmd(.md)        kallisto pseudo-alignment, superseded by HISAT2
│       └── 02_isoform_gene_map_cds.Rmd the CDS-based isoform map, superseded by step 04
├── 02_data/
│   ├── sample-submission/              Tag-seq sequencing submission paperwork
│   ├── annotation_mt_like_loci.csv     mitochondrial copies in the annotation (read by steps 04, 05 and 06 step 01)
│   ├── mitogenome_NC_007687.1.fa       the mitochondrial genome (read by step 08)
│   ├── mitogenome_genes.saf            its 38 genes, rRNAs and tRNAs (read by step 08)
│   ├── strg_gene_ids.csv               StringTie gene IDs of 284 reference transcripts
│   └── (downloads of steps 04 to 06)   Iso-Seq FASTA, genome, GFF, GTF, CDS; git-ignored (README)
└── 03_analyses/
    ├── hisat/                          reference StringTie tables + MultiQC alignment reports
    ├── isoform-gene-map/               step 04: each isoform's genome gene
    ├── augmented-annotation/           step 05: option B's annotations (GFF and SAF git-ignored)
    ├── genome-recount/                 step 06: six gene count matrices of the realignment
    ├── featurecounts/                  the count matrix of record (the DE input)
    ├── mitogenome/                     mitochondrial reads on the mitogenome alone (06 step 13's input)
    ├── prepDE/                         transcript and gene count matrices of the HPC alignment (previous record)
    ├── fastqc/{trimmed,untrimmed}/     FastQC per-sample read QC
    ├── read_trimming/                  read retention and the trimming recipe check
    ├── _superseded/kallisto/           kallisto quant per sample
    ├── _superseded/02_isoform-gene-map_cds/   the retired CDS-based isoform map
    └── knit_html/                      runner reports and logs (git-ignored)
```

## External inputs (not in repo)

| Input | Location |
|-------|----------|
| Raw (untrimmed) Tag-seq reads, all 131 libraries, two lanes each | owl: https://owl.fish.washington.edu/nightingales/M_trossulus/ (gannet `20220405-tagseq/` has 73 of them) |
| Trimmed reads (`*_L099_R1_cmb.trim.fastq.gz`), all 131 libraries, read by `01_hisat_stringtie.Rmd`, steps 06 and 08 (online) and `04_iso-seq-transcriptome` step 02 | gannet: https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/byssus-exp-analysis/data/raw-trimmed/ (recipe above) |
| Genome assembly + annotation | NCBI `GCF_036588685.1` (downloaded by step 01 into `02_data/ncbi_dataset/` on the HPC, and by step 04 into `02_data/`; git-ignored) |
| Iso-Seq transcriptome, `Mtros-hq_transcripts.fasta` | owl: https://owl.fish.washington.edu/halfshell/genomic-databank/ (downloaded by step 04 into `02_data/`, git-ignored) |

`01_hisat_stringtie.Rmd` and `_superseded/07-HiSat_GL.Rmd` invoke HISAT2 / StringTie at fixed
`/home/shared/...` paths and read inputs not committed here, so they do not run as they are
outside that HPC environment. Steps 04 to 06 and step 08 with `online: true` need the tools
above (the runner's `minimap2` and `bin` parameters). Packages for steps 02 to 08: base R,
here, rmarkdown, data.table (04 to 08), GenomicRanges and GenomicAlignments (04, 05),
rtracklayer and Biostrings (04); step 03 with `online: true` also needs curl, perl, FASTX-Toolkit 0.0.14 and
cutadapt (5.2 used).
