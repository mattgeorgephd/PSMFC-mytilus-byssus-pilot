# genome-foot-sprot2026_03-noseg

Written by `../../01_code/05_byssal_noseg_search.Rmd` with `run: true` (cloud environment,
2026-10-06, run again 2026-10-17; R 4.6.1, BLAST+ 2.15.0). The annotation the analysis reads: the
genome blast of 2026 (`../genome-foot-sprot2026_03/`, step 01) with eight genes added that it
left without a hit: six because blastx's low-complexity filter (SEG) masked their low-complexity
adhesive proteins, and two mfp-1-like genes adopted by their RefSeq name (2026-10-17).

Every CDS of RefSeq GCF_036588685.1 (53,269) was searched against the 241 foot and byssal
proteins of `../../02_data/uniprotkb_mytilus_foot_2026_03_byssal.fasta` with step 01's settings
(`-evalue 1E-20 -max_target_seqs 1`) plus `-seg no`, and `-dbsize 209093815`, the total residues of
step 01's database, so that the e-values are those of a search of that whole database. A hit is
adopted only for a gene with no hit in step 01's search, and only at 70% identity or more, or,
whatever the identity, when RefSeq's own name for the gene names a byssal protein (`rule`
"RefSeq name"; added for completeness on 2026-10-17).

| File | Contents |
|---|---|
| `LOC_GO_list.txt` | step 01's `LOC_GO_list.txt` line for line (41,588 rows), then the 14 hits of the 8 adopted genes (one each for the six, 3 and 5 HSPs for the two mfp-1-like genes), built as step 01 builds its rows (same columns, same order; UniProt records of release 2026_03 from `../../02_data/uniprot_mytilus_foot_r2026_03.tsv`). Read by `04_sequence-alignment` steps 04 and 05 and `05_differential-expression` step 01 (the mitochondrial loci, `tools/mt_encoded.R`), `05_differential-expression` step 06, `06_iso-seq-transcriptome` step 04, `07_enrichment` step 01 and `09_gene-mechanics-correlation` steps 01 and 04 |
| `noseg_genes.csv` | one row per gene with a hit in this search (480): its best hit here (highest bitscore: accession, protein, organism, identity, alignment length, e-value, bitscore), the RefSeq name of its protein, its best hit in step 01's search if it has one (451 do), `rule` ("identity" or "RefSeq name"), `adopted` and `reason` |
| `foot_noseg_blastx.tab` | the blastx table of this search (1,784 hits, outfmt 6) |
| `foot_noseg_cds.tsv` | the FASTA header of every CDS in that table and of every CDS step 01 found on a foot protein, so that the tables can be rebuilt without the genome CDS (`run: false`) |
| `adopted_full_db_check.tsv` | the adopted genes' CDSs searched against step 01's whole database (`swissprot_n_mytilusfoot_r2026_03`) with the filter off: the same best protein, identity, e-value and bitscore for each |
| `RUN_provenance.txt` | BLAST+ and its settings, the counts, the code commit, the input MD5s |

The adopted genes (RefSeq name of the protein; best hit; identity over the alignment):

| Gene | RefSeq protein | Hit | Identity |
|---|---|---|---|
| LOC134721757 | uncharacterized protein | Q8WTE8, polyphenolic phosphoprotein mefp-5 (*M. edulis*) | 86.7% over 98 |
| LOC134721760 | neuropeptide-like protein 31 | Q2VPV4, foot protein 3 variant 4 (*M. californianus*) | 80.8% over 78 |
| LOC134721759 | neuropeptide-like protein 30 | A0A8B6DS81, foot protein-3 (*M. galloprovincialis*) | 78.5% over 79 |
| LOC134705410 | protein kreg-1-like | A0A223HCK0, foot protein 9 variant 2 (*M. californianus*) | 78.4% over 125 |
| LOC134709324 | uncharacterized protein | A5HKI4, thread matrix protein 2F (*M. galloprovincialis*) | 74.1% over 166 |
| LOC134721758 | uncharacterized protein | D3JVF0, foot protein 3 variant 5 (*M. coruscus*) | 70.5% over 78 |
| LOC134723087 | adhesive plaque matrix protein-like | Q2TCL0, foot protein 1 variant 1 (*M. californianus*) | 37.5% over 555 (RefSeq name) |
| LOC134723088 | adhesive plaque matrix protein-like | Q25434, adhesive plaque matrix protein (foot protein 1, *M. coruscus*) | 35.8% over 777 (RefSeq name) |

In step 01's search (filter on) the six genes' CDSs had best e-values of 7.5e-20 to 1.8e-9,
above the 1e-20 cut-off; with the filter off they are 5.6e-22 to 1.05e-28. The four mfp-3 and
mfp-5 genes lie together on NC_086378.1, and so do the two mfp-1-like genes (e-values 1.7e-38 and
1.5e-46; mfp-1 is made of variable decapeptide repeats, so a low identity does not argue against
it). In step 01's whole database, searched with the filter off, every adopted CDS has the same
best protein and e-value (`adopted_full_db_check.tsv`). The two mfp-1-like genes have almost no
reads in the day-3 libraries (LOC134723088's reads come mostly from the day-0 libraries), so they
pass no expression filter and enter no test; they appear in the byssal lists of
`09_gene-mechanics-correlation` step 04, flagged below the filter.

Not adopted (see `noseg_genes.csv` and the step's text): the 21 other genes with no hit in step 01
match at 29 to 45% identity, mostly through collagen-like, mfp-1 decapeptide and mfp-4
histidine-rich repeats, and RefSeq names none of them a byssal protein; and of the 451 genes with a hit in step 01,
142 score higher here (105 on the same protein, 37 on another), which is not used, since scores
with and without the filter are not comparable.
