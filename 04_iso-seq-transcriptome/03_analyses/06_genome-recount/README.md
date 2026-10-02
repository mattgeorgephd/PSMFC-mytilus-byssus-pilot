# 06_genome-recount

Written by `../../01_code/06_genome_recount.Rmd` (option B, a parallel analysis): the 131
trimmed Tag-seq libraries realigned to the genome and counted on each of step 05's three
annotations by two counters, six gene count matrices in all.

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
