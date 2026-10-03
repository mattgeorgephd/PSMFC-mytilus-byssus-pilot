# 02_isoform-gene-map

Written by `../../01_code/02_isoform_gene_map.Rmd`: each Iso-Seq isoform placed on the genome
(GCF_036588685.1) with minimap2 2.31 (`-c -x splice:hq -uf --secondary=no`) and given the
annotated gene (RefSeq release RS_2024_02, every gene type) whose exons it shares most bases
with on its strand.

| File | Contents |
|---|---|
| `isoform_gene_map.csv.gz` | one row per isoform (411,251): `gene` (a gene key, `novel:<sequence>:<start>-<end>:<strand>`, `mito:<isoform>` or `isoseq:<isoform>`), `status` (`assigned`, `ambiguous`, `novel_locus`, `mitochondrial`, `unassigned`), the alignment (`seqid`, `strand`, `start`, `end`, `identity` as gap-compressed identity, `mapq`, `n_alignments`), `exonic_overlap` with the gene, the runner-up gene and its overlap, `mt_hit`, `mito_cov` (share of the isoform aligned to the mitogenome), `isoform_len`, `three_prime_beyond` (bases by which an assigned or ambiguous isoform's 3' end lies beyond its gene's annotated 3' end; read by step 03) |
| `isoform_gene_map_summary.csv` | isoforms and features by status; genome-branch genes reached |
| `cds_map_agreement.csv` | agreement with the retired CDS-based map (`../_superseded/02_isoform-gene-map_cds/`), nuclear genes only |
| `RUN_provenance.txt` | settings, minimap2 version, input MD5s. The step ran from a6ebbc8, before the count matrix of record changed (2026-10-02), so the MD5s it records for `06`'s `gene_count_matrix_clean.csv` and `mitochondrial_loci.csv` are those of the previous files (in the git history at 303b4b0). The change does not alter the map: with the new files the step reaches the same 19,856 genes, and the 21 loci added to the mitochondrial list are tRNAs and rRNAs of the mitogenome, to which no isoform is assigned (isoforms there are flagged by sequence) |
| `minimap2_isoforms_vs_genome.paf.gz`, `minimap2.log`, `mitogenome.fa`, `minimap2_isoforms_vs_mitogenome.paf.gz` | the alignments (git-ignored; the genome alignment takes about 45 minutes and 10 GB of memory on four cores) |

## Result (2026-10-02)

| status | isoforms | features |
|---|---|---|
| assigned to a gene | 309,378 | 19,425 genes |
| ambiguous (kept on the best gene) | 2,021 | 897 genes |
| novel locus (aligned well, outside every gene) | 62,811 | 16,685 loci |
| mitochondrial | 2,456 | left out in step 04 |
| unassigned (20,317 unaligned, the rest below 0.90 identity) | 34,585 | one each |

Isoforms reach 19,856 of the genome branch's genes. Where both this map and the CDS-based map
put an isoform on a nuclear gene (185,383 isoforms), they agree on the gene for 97.0%; the
CDS map reached 17,733 nuclear genes, this map 19,856, 17,074 of them in common.

Two choices shaped these numbers:

- **Identity is gap-compressed** (one minus minimap2's `de`). Matching bases over block length
  counted every base of a long unspliced gap as a difference and left 75,199 isoforms below
  0.90 in a first run, 60,907 of which are at 0.90 or better by gap-compressed identity.
- **Mitochondrial isoforms are set aside first.** Of the 2,456 mitochondrial isoforms, 594 align
  to the mitogenome itself; 1,858 were placed on unplaced scaffolds (mostly unannotated copies
  of mitochondrial sequence, where they would have formed novel loci, or a listed copy) but
  are at least half mitogenome sequence; 4 fall mostly on a listed mitochondrial locus. 277
  more isoforms carry a smaller mitochondrial segment (`mt_hit`): 247 sit in novel loci and 27
  are unassigned, so step 04 leaves out the features holding them; 3 keep their gene.
