# 02_isoform-gene-map

Written by `../../01_code/02_isoform_gene_map.Rmd`: each Iso-Seq isoform assigned to a gene of
the genome annotation through its best same-strand alignment to the annotation's CDS
(minimap2 2.31, `-c -x asm20`; identity at least 0.90 over at least 100 aligned bases).

| File | Contents |
|---|---|
| `isoform_gene_map.csv.gz` | one row per isoform (411,251): `gene` (LOC key, or `isoseq:<isoform>` when unassigned), `status` (`assigned`, `ambiguous` when a different gene scores at least 98% of the best, `unassigned`), best CDS, identity, aligned bases, isoform and CDS lengths, alignment score, runner-up gene and score |
| `isoform_gene_map_summary.csv` | isoforms and genes by status |
| `RUN_provenance.txt` | settings, minimap2 version, input MD5s |
| `minimap2_isoforms_vs_cds.paf.gz`, `minimap2.log` | the alignments (git-ignored; rebuilt by the script in about 3 minutes) |

Result: 182,861 isoforms assigned and 6,910 ambiguous (189,771, 46%) reach 17,792 of the
34,235 genes that have a CDS; 221,480 (54%) stay unassigned. The identity cut-off is not what
leaves them out: only 4,816 isoforms have a best hit between 0.85 and 0.90, and about 215,000
have no same-strand CDS hit of 100 bases at all (non-coding transcripts, and 3' fragments
that end before the CDS). 610 have only an antisense hit.
