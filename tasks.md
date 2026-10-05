# tasks.md

Current sprint: make the analysis run end to end, from thread strength and Tag-seq counts
to GO enrichment across treatments and the gene-mechanics associations. Updated 2026-10-05.
Conventions and how to run are in `AGENTS.md`.

The analysis folders were renumbered on 2026-10-03 so that their numbers are the run order:
`05_sequence-alignment` became `04_`, `06_differential-expression` became `05_` and
`04_iso-seq-transcriptome` became `06_` (root `README.md`, "What the map shows about the
order", maps the steps). Entries under Done written before then keep the numbers of their time;
the other sections use the new ones.

## Done

- **Paths and folders.** Every script points at the renumbered folders; READMEs updated;
  references to the retired candidate / worklist / pictures workflow removed.
- **Repository hygiene.** `.gitignore` for RStudio caches, rendered reports and large
  alignment files; the NCBI API key removed from the code (it now comes from `ENTREZ_KEY`);
  the original `tools/pipeline_checks.R` restored.
- **Batch runners.** `00_run_*.Rmd` in 02, 05, 06, 07, 08 and 09 and `00_run_pipeline.Rmd` at
  the root, on the shared `tools/run_steps.R`; steps renumbered `01_` onwards in every folder;
  each script writes only to its own `03_analyses/`.
- **02 thread strength.** Scripts 01-05 with numbered output folders; the mussel key moved to
  `03_analyses/`; regenerated results match the committed ones.
- **05 alignment to DE handoff.** `_prepde.R`, an R port of `prepDE.py3` (byte-identical on
  test data), and `02_prepDE.Rmd`, which rebuilds the gene matrix from committed data.
- **06 differential expression.** The contrasts defined by rule in one script (16 at first,
  reproducing the old sample lists exactly; 7 since the LC retirement); one saved DESeq2 fit per contrast, so shrinkage runs
  on its own; LC and foot-vs-gill results drawn; TC manuscript figures saved to `figures/`.
- **Foot regions.** FX libraries recorded as foot (the rest of the foot, without the phenol
  gland) with a `region` column and a crosswalk to the RNA isolation log; they enter no
  contrast.
- **07 GO enrichment.** topGO `weight01` (of record), goseq, clusterProfiler `enricher` with a
  compareCluster view, rrvgo reduction and a method comparison; each contrast's tested genes as
  universe, up and down separately. Replaces DAVID / REVIGO (kept in `_superseded/`).
- **08 GO slims.** Rewritten: the old script used only each gene's first GO ID; slim pinned in
  `02_data/goslim_generic.obo`; heatmap added.
- **09 gene mechanics.** Scripts 01-04 read 06's new outputs and reproduce the previous code's
  results exactly; runner loops over foot and gill.
- **Figures.** One colour set for the whole repository (`tools/plot_style.R`), checked for
  colour-vision deficiency; the stressor Venn labels no longer overlap; the DEG-count figure's
  x axis is titled Tissue (it read Treatment).
- **Provenance.** `07` records the GO release and package versions behind its results
  (`03_analyses/01_go-inputs/RUN_provenance.txt`), as `02` and `09` already did for theirs.
- **08 NCBI summaries.** Step 02 records a failed request as `Error: <message>`; it used to drop
  the accession silently. Two committed summary tables (FOA, GOA) turned out to come from an
  earlier top-50 list; documented in `Top_gene_summaries/README.md`.
- **Documentation.** A README in every new folder (and in the folders that were missing one);
  `AGENTS.md` and this file.
- **Mitochondrial genes on their own.** The 310 mitochondrial loci (12 protein and 5 RNA genes
  of the mitochondrial genome, 293 nuclear-scaffold copies) are out of the genome DE, GO and
  gene-mechanics universes and tested per protein in `06` step 13, with a manuscript figure
  (`figures/MT_mitochondrial_expression.png`): Gill OA raises 10 of 12 proteins, Foot OA 3,
  OW and DO none.
- **Lab controls retired.** The day-0 animals are not used as a control anywhere (different
  foot dissection); LC contrasts and their GO results moved to `_superseded/` in 06 and 07.
- **Thread metrics.** `mean_force` (mean of the animal's thread peak forces) and `max_force`
  (its strongest thread) replace the old per-animal `max_force`; extension removed from 02 and
  09 (thread length under test not controlled).
- **GO terms in gene mechanics.** `09` step 05 tests DEG sets, every enriched GO term and the
  mitochondrial share against the thread metrics.
- **All GO figure options drawn** for every family and ontology (topGO, goseq, clusterProfiler
  dot plots, rrvgo parents, method comparison).
- **Naming audit.** One gene key (`tools/gene_ids.R`); best BLAST hit by bitscore everywhere
  (06 top-50 labels, 06 Zenodo lists and 08 GO slims took the first-listed hit); name-based
  column selection in 06 step 06; GO term labels in figures cut, not edited, and kept unique;
  topGO gene lists indexed by GO ID. Checked: the BLAST table's columns are consistent row by
  row, and every gene listed under an enriched topGO term is annotated to it.
- **Figure fixes.** Failure-mode colours (greys), clipped titles, DEG-count panel spacing, the
  mitochondrial figure's estimates and intervals now from the same Wald model as its p-values.
- **Integral removed.** The area under the force-time curve is no longer extracted (it depends
  on the same unstandardized thread length as extension).
- **Top-50 labels.** Bars are labelled with the gene symbol of the best UniProt hit (locus tags
  skipped, protein name as fallback) instead of the entry-name mnemonic; the six bar plots are
  saved (`06/03_analyses/top_DEGs/Top_50_genes/<code>_top50.png`).
- **08 step 03.** A failed OrthoDB request stopped the whole step (`next` in an error handler);
  it now records the error and continues (tested offline).
- **Decisions recorded.** Outlier-replaced genes: the standard DESeq2 convention is kept. The
  manuscript GO figure is chosen later; every option stays drawn.
- **Read trimming record.** The trimmed reads every count rests on (`*_L099_R1_cmb.trim`) were
  made with the Matz-lab Tag-seq clipper (leader removal and PCR-duplicate removal) followed by
  FASTX-Toolkit `fastx_clipper` (poly-A, adapter, 20-base minimum); `05` step 03 reproduces
  99.8% of a sample read for read and tabulates read retention per library (median 43% of raw
  reads kept). Matt's August 2022 cutadapt script, found on gannet `seashell/bu-github/`, is
  kept as a record (`05/01_code/_superseded/`); it and the August 2022 MultiQC report describe
  a different trimming, not the reads analysed. The raw reads of all 131 libraries are on owl
  (`nightingales/M_trossulus/`); gannet's `20220405-tagseq/` is a partial copy (73 libraries).
- **Iso-Seq branch implemented.** Isoforms placed on the genome (minimap2 spliced alignment, exon
  overlap with RefSeq genes; mitochondrial isoforms set apart), the 131 libraries quantified
  with salmon (no length correction) and summed per gene with tximport, the six TC contrasts
  refitted as in `06` and compared: the genome DEGs largely replicate (54 to 80%, same
  direction), and the Iso-Seq reference finds about twice as many (`04` README).
- **Mitochondrial pseudogenes left out.** The mitochondrial copies had been found through a BLAST
  of the annotation's CDS, which cannot see the 167 pseudogenes NCBI names after mitochondrial
  proteins (1.17 million reads; 81 of the 543 Gill OA DEGs). They are now found by name
  (`05_sequence-alignment/02_data/annotation_mt_like_loci.csv`, which reproduces the 126
  BLAST-found copies exactly) and handled like the other copies. Gill OA has 423 DEGs, Foot OA
  75; Gill OA raises 10 of the 12 mitochondrial proteins. Summing every copy also removed the
  one outlier replacement among the proteins (T025G's ND1 and ATP6, whose excess came from how
  the reads split between copies).
- **08 NCBI and OrthoDB tables refreshed.** With NCBI and OrthoDB reachable, steps 02-03 ran
  on the current top-50 lists (rerun on 2026-10-02 after the pseudogene change altered them). Step 02 now finds each accession's gene through NCBI Protein and
  the protein-to-gene link (free-text search had returned the wrong gene for SPDEF); step 03
  could not run as written (it joined OrthoDB's headerless species table on a missing column)
  and now uses OrthoDB 12.2.
- **Figure reproducibility.** The two PCA figures with two untitled legends (06 step 03) changed
  on every run: ggplot2 4.0.3 orders untitled legends at random from one R session to the next
  (checked in six sessions). Their legend order is now fixed. 06 step 05's provenance now names
  the commit, like the others.

- **Option B as a parallel analysis.** The RefSeq annotation augmented with the Iso-Seq
  isoforms (3' extension; full models), the 131 libraries realigned (HISAT2) and counted on the
  RefSeq control and both augmented annotations with StringTie + prepDE and featureCounts, and
  the six TC contrasts fitted on each and compared with the record (`04` steps 05 to 07).
- **Option B adopted (2026-10-02, Matt's decision).** The count matrix of record is now
  featureCounts on the 3'-extended annotation (`05` step 04, from `04` step 06); `06` to `09`
  rerun on it. DEGs: Foot OA 161, Foot OW 211, Foot DO 678, Gill OA 584, Gill OW 227, Gill DO
  404 (70 to 90% of the previous DEGs kept; `04 .../07_augmented-de/record_change.csv`). `04`
  step 07 now compares every recount with the previous record and checks that its fit of the
  record reproduces `06`.
- **Mitochondrial genes counted on the mitogenome alone (2026-10-02, Matt's decision).** In the
  genome alignment the 293 nuclear copies and unannotated mitochondrial-like sequence took
  about half the mitochondrial reads (HISAT2 reports at most five placements), and 12 of the 59
  animals carry mitochondrial haplotypes that differ from the reference at fixed positions,
  which the genome alignment undercounted further (COX1 less than half). `05` step 05 now aligns
  every library to NC_007687.1 alone and `06` step 13 tests those counts, with the haplotype
  groups, a haplotype covariate and a permissive alignment as checks. Gill OA raises 10 of 12
  proteins and their sum (1.41-fold); Foot OA raises ND2, ND3 and ND5, its sum 1.19-fold (p
  0.061); OW and DO none. The 21 mitochondrial tRNAs and rRNAs the previous matrix had named by
  StringTie number are now recognised (331 mitochondrial loci).
- **Pull request merged** (2026-10-03): the work above is on `main`.
- **Loop removed** (2026-10-03). `04` steps 02 and 05 found the mitochondrial loci in `06`'s
  list, which `06` builds from the count matrix those steps lead to. They now find the
  annotation's mitochondrial loci themselves with the same function (the same 331 loci);
  rerun from their cached alignments, both gave identical outputs.
- **Step map** (2026-10-03). The root `README.md` ("How the steps connect") maps what every step
  reads and writes, traced from the code, with an overview diagram and one per folder. Folder
  READMEs and comments that named the wrong readers or inputs were corrected (06's count matrix
  source and readers, 07's step inputs, 08's annotation source, 02's ANCOVA metrics, 05's prepDE
  readers, 03's download claims, 04's step 01, `tools/` users).
- **Count-building steps moved** (2026-10-03, Matt's decision: option C, the folder numbers in
  run order). The isoform-to-gene map, the augmented annotations and the genome recount
  (`04` steps 02, 05, 06) are now `05_sequence-alignment` steps 04 to 06, with their outputs,
  caches and the retired CDS map; `05`'s count matrix of record and mitogenome counts are steps
  07 and 08; what stays in `04` is renumbered 01 to 04 (salmon 02, Iso-Seq DE 03, recount
  comparison 04). Both runners rerun offline with identical outputs.
- **Folders renumbered** (2026-10-03, option C). `05_sequence-alignment`,
  `06_differential-expression` and `04_iso-seq-transcriptome` are now `04_sequence-alignment`,
  `05_differential-expression` and `06_iso-seq-transcriptome`, with their caches; every path in
  code, the runners (the pipeline's stages now run 02, 04, 05, 06, 07, 08, 09), `.gitignore`,
  the READMEs, `AGENTS.md` and the step map follow. Every folder now reads only lower-numbered
  folders; `03_blast` and `04` step 01 (HPC) stay a break in the chain.
- **Online steps and the whole pipeline rerun under the new names** (2026-10-03).
  `04_sequence-alignment` steps 04-06 and 08, step 03's recipe check and
  `06_iso-seq-transcriptome` step 02 rerun from their caches, then `00_run_pipeline.Rmd` from a
  clean commit: every step TRUE, `09` checks 0 failed, every output table byte-identical to
  before the renumbering.
- **Genome BLAST script fixed** (2026-10-03). `03_blast/01_code/01_genome_blast.Rmd` runs again
  (the 2024 version is in `01_code/_superseded/`). From the 2024 blastx table it rebuilds
  `g.spid.txt` and, in its `as_2024` mode, `LOC_GO_list.txt` byte for byte; its database was
  Swiss-Prot release 2024_04 (not 2024_01), and its search reproduces the 2024 hits on a
  100-CDS sample.
- **Annotation and software of 2026** (2026-10-04, pull request GO-db-2026). A new genome
  blastx against Swiss-Prot 2026_03 plus the 196 "(mytilus foot)" proteins and 45 byssal
  proteins that query misses (`03_blast/02_data/byssal_additions_2026_03.tsv`), with the
  UniProt 2026_03 records; GO release 2026-01-23 (`GO.db` 3.23.1) and its generic slim; R 4.6.1
  and Bioconductor 3.23 (`renv.lock`, `tools/cloud_setup.sh`). Every consumer reads
  `genome-foot-sprot2026_03/LOC_GO_list.txt`; `03_blast` step 04 (the 2024 hits with the 2026
  records) is kept as a comparison. The byssal rules gained "thread matrix protein",
  "nongradient byssal", "YGH-rich protein" and "protease inhibitor-like protein-1" so that no
  gene leaves the byssal sets for a change of name, and (Matt's decision, 2026-10-04) the
  *M. coruscus* byssus proteins "C1q-domain-containing protein-1" and "TSP_1 domain containing
  protein-1" count as byssal accessory genes in `09`. Whole pipeline rerun from a clean commit:
  every step TRUE, `09` checks 0 failed; mitochondrial loci and DEGs unchanged; 177 of the 296
  TC GO terms of record remain among 285 (most of the change from the new UniProt records).
- **Runner checks** (2026-10-05, pull request GO-db-2026). A run on R 4.2.2 with `GO.db`
  3.16.0 stopped at stage `07` after 37 minutes, and its stage log only pointed to other logs.
  `00_run_pipeline.Rmd` now checks R, Bioconductor and the GO release against `renv.lock`
  before the first stage (`check_stack()`, parameter `check_versions`), and every runner's
  error quotes the failed step's own error (`finish_run()`). Runner code only: no output
  changes. Then (2026-10-05) the check also stops when a package the pipeline loads is not on
  R's library path (a fresh R 4.6.1 whose `renv::restore()` library was not in `R_LIBS`), and the
  tensometer QC plots are no longer committed (383 JPEGs rewritten by every run; about 8 MB of
  history per run from another machine).
- **Candidate keyword** (2026-10-05, pull request GO-db-2026). In `09` the tRNA-synthetase
  keyword `aminoacyl` also matched aminoacylase-1 and acylaminoacyl-peptidase, enzymes of
  N-acyl amino acids and peptides; scripts 01 and 02 now use `aminoacyl[- ]tRNA`. The two genes
  leave both tissues' candidate sets (259 foot and 324 gill candidates tested) and the
  `tRNA_translation` module. Slopes and p-values are unchanged; q-values move by up to 0.09 and
  family q-values by up to 0.11 (BH over two fewer genes: glutathione peroxidase 7 with the
  strongest thread in gill goes from family q 0.14 to 0.25), and the gill HSP70-family association with plaque area
  keeps q 0.036 (family q 0.072). Aminoacylase-1 had been the gill's strongest adhesion
  candidate (q 0.14). `Hsp` still matches abbreviations inside unrelated names (perlecan's
  "HSPG", "HsPDE8B" and others: 6 foot and 9 gill candidates; DOC, caveats); whether to
  tighten it is open.

## In progress

- Proofreading the regenerated figures before they go into the manuscript.

## Blocked or waiting on a decision

| item | why it is blocked | who |
|---|---|---|
| Merge GO-db-2026 | the annotation of 2026 (Swiss-Prot 2026_03 search plus 45 byssal proteins, UniProt 2026_03 records, GO release 2026-01-23) and R 4.6.1 / Bioconductor 3.23 are on the stacked pull request; the DEGs are unchanged, 177 of the 296 TC GO terms of record remain among 285, and the ER stress terms weaken (`07` README). Merge, or keep the 2024 annotation | Matt |
| Revoke the old NCBI API key | it was removed from the code but remains in the public history (commit 3bee414); only its owner can revoke it at NCBI | key owner |
| Manuscript text | numbers and wording predate this pull request: the count matrix (now featureCounts on the Iso-Seq-extended annotation) and every DEG count, the GO results (Gill OA up's ATP synthase terms are no longer FDR-supported), LC contrasts, the mitochondrial method and result (mitogenome counts; COX1 and ND3 now among the significant proteins), thread metric names (`mean_force`, `max_force`), extension | Matt |
| Manuscript GO figure | every option is drawn (`07_enrichment/03_analyses/02_topgo/` to `06_method-comparison/`); choose one | Matt |
| Library QC | three libraries carry most of DESeq2's outlier replacements: T040F (117 of 157 genes refit in Foot OW), T025G (55 of 96 in Gill OA; also the lowest alignment rate, 60.4% against a median of 67.2%) and T035G (56 of 98 in Gill OW). Check their RNA records; a sensitivity run without them would show whether any result rests on them | Matt |
| Exploratory foot-region contrast | rest of foot vs phenol gland to tip in the 12 day-0 animals (3,022 of 7,367 genes differ); add as a DE family if useful for interpreting byssal genes | Matt |
| `03_blast/` and `04` step 01 | not pipeline stages; committed outputs are used. `04` step 01 is HPC only (inputs and binaries not in the repository). `01_genome_blast.Rmd` runs again since 2026-10-03 (the 2024 version is in `01_code/_superseded/`); its search of 2026 ran in the cloud environment (BLAST+ 2.15.0, about 31 hours on 4 threads, in resumable parts), and steps 02 and 03 still need the HPC | HPC user |

## Known limitations (documented, not blocking)

- **Software environment only partly locked.** The pipeline was tested with R 4.6.1 and
  Bioconductor 3.23; `renv.lock` (2026-10-03) records all 277 R packages of that library, and each
  `RUN_provenance*.txt` the versions a step used, but the system libraries and the R build are not
  locked (no container image). System
  libraries matter too: after the cloud container was rebuilt on 2026-10-01, every figure
  re-rendered with different fonts (same data, different bytes).
- **Annotation provenance.** The analysis reads `LOC_GO_list.txt` of the genome blastx search of
  2026 (`03_blast/01_code/01_genome_blast.Rmd`, `genome-foot-sprot2026_03/`: Swiss-Prot 2026_03
  plus 241 foot and byssal proteins, UniProt 2026_03 records). The search of 2024 (Swiss-Prot
  2024_04, not 2024_01 as its file name said; its `V11` column is the e-value of the gene's
  first hit in that same search, not a second run) is kept in `genome-foot/`. `g.spid.txt`
  keeps one transcript per gene and protein, so the rows of a gene's other transcripts that hit
  the same protein carry no LOC key (16,842 of 41,588 rows in 2026); each gene keeps its first
  transcript's hits. The script keeps that rule so its tables match the 2024 ones; the CDS
  FASTA headers map every CDS to its LOC and would let those rows be keyed. 6 of the 14,855 GO
  IDs are unknown to the pinned GO release and dropped.
- **GO of record.** topGO `weight01` p < 0.01 is unadjusted; in BP only five TC runs have
  terms under FDR control (clusterProfiler; four with goseq). Headline GO claims are safest from
  `consensus_terms_TC_<ont>.csv`.
- **Repository size.** `.git` is about 765 MB, mostly committed superseded tables and BLAST
  outputs; `vst_paired_<T>.csv` (7 to 11 MB each) is rewritten on every 09 run.
- **No continuous integration.** Nothing runs the pipeline on a push; the checks are the
  runners, the provenance files and the `09` checks.
- **Numbers computed outside the pipeline.** The Iso-Seq vs genome correlation of the
  mitochondrial share (Spearman 0.91; `06` README) is computed by hand from committed files;
  no step writes it.
- **Inputs no current step writes:** `04 .../hisat/t_data.ctab` (one sample's table from an older HPC run),
  `04 .../prepDE/transcript_count_matrix.csv` (HPC `prepDE.py`) and `04 .../fastqc/` (the
  earlier `byssus-exp-analysis` repository).

## Up next

- Update the manuscript text.
- Library QC sensitivity run (T025G, T035G, T040F).
