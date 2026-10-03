# augmented-annotation

Written by `../../01_code/05_augmented_annotation.Rmd` (option B): the RefSeq annotation of
GCF_036588685.1 (release RS_2024_02) augmented with the Iso-Seq isoforms, in three versions
that step 06 counts the Tag-seq reads on. `ext3` is the annotation of record since 2026-10-02
(counted by featureCounts; `06_iso-seq-transcriptome/03_analyses/04_augmented-de/README.md`).

| annotation | what it is |
|---|---|
| `refseq` | the RefSeq GFF unchanged (byte for byte), the record's annotation: the control |
| `ext3` | each RefSeq transcript's last exon extended to the furthest 3' end of a same-gene isoform whose last aligned block overlaps it (same strand), stopping before any other gene on that strand; no other change |
| `full` | RefSeq plus each isoform as an extra transcript of its gene, plus the novel loci as genes; a novel locus overlapping one gene's exons (RefSeq or added) on its strand joins that gene, one overlapping several genes or a mitochondrial locus is left out |

Isoforms used: those step 04 `assigned` to a gene with no exon overlap with a second gene, and
the isoforms of novel loci without a mitochondrial segment. Mitochondrial isoforms are never
added; the mitochondrial loci (the mitogenome's 38 features and the 293 nuclear copies) are as RefSeq has them in every version.

| File | Contents |
|---|---|
| `ext3_extensions.csv.gz` | every transcript `ext3` extends: `transcript_id`, `gene` (gene key), `seqid`, `strand`, old and new 3' end (`old3`, `new3`), `extension_bp`, `isoforms` overlapping its last exon, `capped` (stopped before another gene) |
| `full_added_transcripts.bed.gz` | the transcripts `full` adds, BED12, name `<gene>|<isoform>` (rebuilds `full` without step 04's alignment) |
| `novel_loci_fate.csv` | step 04's novel loci in `full`: joined a gene, left out, or kept as new genes |
| `annotation_summary.csv` | per annotation: genes, transcripts, exonic bases (each gene's exons merged), genes changed |
| `RUN_provenance.txt` | rules, input MD5s |
| `annotation_<refseq|ext3|full>.gff`, `annotation_<refseq|ext3|full>.saf` | StringTie's and featureCounts's inputs (git-ignored; 0.3 to 0.5 GB each, rebuilt in a few minutes) |

## Result (2026-10-02)

| annotation | genes | transcripts | exonic bases | vs RefSeq | genes changed |
|---|---|---|---|---|---|
| refseq | 47,806 | 68,103 | 89.1 Mb | 1 | 0 |
| ext3 | 47,806 | 68,103 | 93.2 Mb | 1.046 | 9,781 |
| full | 60,293 | 430,079 | 137.3 Mb | 1.541 | 31,470 |

- `ext3` extends 18,470 transcripts of 9,781 genes, by a median 76 bases (interquartile 6 to
  622; at most 6,874). No extension reached another gene on the same strand: isoforms that
  touch a second gene's exons were left out beforehand, so the cap never applied. 570
  extensions overlap a gene on the opposite strand, which matters only for the unstranded
  StringTie count.
- `full` adds 303,668 isoforms to 18,908 genes. Of step 04's 16,661 novel loci (after leaving
  out the 24 with a mitochondrial segment), 3,987 (31,102 isoforms) overlap the exons of one
  gene, mostly where that gene's own isoforms extend it, and join it; 187 overlap several
  genes or a mitochondrial locus and are left out; 12,487 stay new genes. Without this step
  the loci inside extended genes split their reads with them (in a test library,
  featureCounts left 122,596 reads unassigned as ambiguous against 14,562 on `refseq`; with
  it, 22,162).
- The isoforms' last blocks are what `ext3` uses; isoforms whose alignment continues past the
  gene across a long gap (up to hundreds of kilobases in step 04's table) extend nothing in
  `ext3` and add their blocks as exons in `full`.
