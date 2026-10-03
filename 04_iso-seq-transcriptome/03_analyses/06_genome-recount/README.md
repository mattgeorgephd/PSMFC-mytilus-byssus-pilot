# 06_genome-recount

Written by `../../01_code/06_genome_recount.Rmd` (option B): the 131 trimmed Tag-seq libraries
realigned to the genome and counted on each of step 05's three annotations by two counters,
six gene count matrices in all. `featurecounts_ext3_gene_counts.csv.gz` is the count matrix of
record since 2026-10-02 (`05_sequence-alignment` step 04 names its rows and passes it to `06`).

- **Alignment:** HISAT2 2.2.1 with the record's settings (defaults and `--dta`). The record's
  index held the RefSeq splice sites and exons (`hisat2-build --ss --exon`), which needs more
  memory than the 15 GB machine used here; the index here is the plain genome, with the same
  splice sites (from NCBI's GTF, `hisat2_extract_splice_sites.py`) given at alignment time
  (`--known-splicesite-infile`).
- **StringTie + prepDE:** StringTie 2.2.1 `-e -G <annotation>`, as the record, and prepDE's
  gene counts (read length 75) by `../../01_code/_prepde_sample.R`, which matches the port in
  `05_sequence-alignment` exactly (`checks.csv`). Reference transcripts StringTie files under
  its own `STRG.N` gene IDs (which change between libraries) are named by their transcript.
- **featureCounts** (Subread 2.1.1): each annotation's exons, sense strand (`-s 1`), uniquely
  aligned reads on exactly one gene.

Each library was processed by `../../01_code/_recount_library.sh` (the step's own loop, and,
to use the machine's idle cores while StringTie parsed the large annotation on one core, three
more workers running the same script with the same arguments on separate libraries).

| File | Contents |
|---|---|
| `<stringtie|featurecounts>_<refseq|ext3|full>_gene_counts.csv.gz` | genes x libraries (131; column names as in `06`'s count matrix, `T001F`); `gene` is prepDE's name for StringTie (gene_id\|gene_name) and the gene key for featureCounts |
| `mapping_summary.csv` | per library: reads, HISAT2's unique and multiple alignments and overall rate, the record's rate for the same library, featureCounts's assigned reads and StringTie's total count on each annotation |
| `checks.csv` | the per-library counter against the `05` port; the control library (T014F) against the record's StringTie output, transcript by transcript, and the record's count matrix |
| `RUN_provenance.txt` | tool versions and settings, input MD5s |
| `hisat2_index/`, `hisat2-build.log`, `refseq_splice_sites.txt`, `libraries/<library>/` | the index (2 GB), the splice sites and the per-library counts and logs (git-ignored) |

## Result (2026-10-02)

- **Alignment.** HISAT2 here aligns a median 0.45 percentage points fewer reads than the
  record for the same library (range 0.14 to 0.63; foot 69.7 against 70.1%, gill 63.3 against
  63.7%), most likely because the splice sites are given at alignment time instead of being
  built into the index (not tested further: the record's index cannot be built here).
- **The control against the record** (T014F, `checks.csv`): 60,491 of the 68,103 transcripts
  get exactly the record's coverage, many of them zero in both (Spearman 0.948 over all); the gene total is 0.993 of the
  record's, and 90% of the genes with at least 10 counts are within 10%. The 21 record genes
  named only by a StringTie `STRG` number cannot be matched (those numbers change between runs).
- **Reads counted** (median per library, relative to the same counter's RefSeq count): `ext3`
  +4.5% (StringTie + prepDE) and +4.4% (featureCounts); `full` +28.1% and +18.4%.
  featureCounts assigns a median 35.5% of the reads on RefSeq: the rest are unaligned (about
  30%), aligned to several places (left out) or outside every annotated exon.
