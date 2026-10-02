# tasks.md

Current sprint: make the analysis run end to end, from thread strength and Tag-seq counts
to GO enrichment across treatments and the gene-mechanics associations. Updated 2026-10-01.
Conventions and how to run are in `AGENTS.md`.

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

## In progress

- Review and merge of the pull request carrying this work (branch `claude/jolly-clarke-7afl1v`).
- Proofreading the regenerated figures before they go into the manuscript.

## Blocked or waiting on a decision

| item | why it is blocked | who |
|---|---|---|
| Revoke the old NCBI API key | it was removed from the code but remains in the public history (commit 3bee414); only its owner can revoke it at NCBI | key owner |
| Merge the pull request | branch `claude/jolly-clarke-7afl1v`; merging deletes untracked `.Rhistory` copies from local checkouts (see the PR) | Matt |
| Adopt option B? (Iso-Seq-augmented annotation for the genome analysis) | built and run as a parallel analysis (`04` steps 05 to 07). Realignment alone moves about 5% of DEG calls. The 3' extension (StringTie + prepDE) keeps 91 to 97% of the record's DEGs, adds 21 to 67 per contrast (mostly near misses) and keeps gene keys, names and GO terms; 30 gene-contrast pairs with large fold-change shifts need a look first. The full models keep 68 to 84% and add 57 to 183, some on unnamed loci. featureCounts finds 16 to 96% more DEGs than StringTie + prepDE but drops most mitochondrial reads (multi-mapped to the nuclear copies). Adopting any version means switching `05`'s counts and rerunning `06` to `09`; evidence in `04_iso-seq-transcriptome/03_analyses/07_augmented-de/README.md` | Matt (decision: keep the record, adopt `ext3`, adopt `full`; and separately the counter) |
| Manuscript text | numbers and wording predate this pull request: Gill OA DEGs (711 to 423), the GO results, LC contrasts, the mitochondrial result, thread metric names (`mean_force`, `max_force`), extension | Matt |
| Manuscript GO figure | every option is drawn (`07_enrichment/03_analyses/02_topgo/` to `06_method-comparison/`); choose one | Matt |
| Library QC | three libraries carry most of DESeq2's outlier replacements: T040F (90 of 122 genes refit in Foot OW), T025G (48 of 77 in Gill OA; also the lowest alignment rate, 60.8% against a median of 67.7%) and T035G (43 of 84 in Gill OW). Check their RNA records; a sensitivity run without them would show whether any result rests on them | Matt |
| Exploratory foot-region contrast | rest of foot vs phenol gland to tip in the 12 day-0 animals (3,022 of 7,367 genes differ); add as a DE family if useful for interpreting byssal genes | Matt |
| `03_blast/` and `05` step 01 | HPC only (inputs and binaries not in the repository); committed outputs are used. `01_genome_blast.Rmd` would not run as written (`mytfoot_uniprot_info` read but `mytilus_foot_uniprot_info` used; `masterID` never defined) | HPC user |

## Known limitations (documented, not blocking)

- **No locked software environment.** The pipeline was tested with R 4.4.3 and Bioconductor
  3.20 (versions recorded in each `RUN_provenance*.txt`), but the repository has no `renv.lock`
  or conda environment file, so a new machine has to assemble the packages by hand. System
  libraries matter too: after the cloud container was rebuilt on 2026-10-01, every figure
  re-rendered with different fonts (same data, different bytes).
- **Annotation provenance.** `LOC_GO_list.txt` joins two BLAST runs (its second e-value column
  comes from the run against Swiss-Prot release 2024_01 plus a UniProtKB "mytilus foot" query,
  `03_blast/01_code/01_genome_blast.Rmd`) and dropped second transcripts of a gene that hit the
  same protein (40% of rows carry no LOC key; each gene keeps its other transcripts' hits). The
  CDS table on owl (`data/ncbi_dataset/data/GCF_036588685.1/cds_from_genomic.tab`) maps every
  CDS to its LOC and would let those rows be keyed. 42 of the 14,863 GO IDs are unknown to the
  pinned GO release and dropped.
- **GO of record.** topGO `weight01` p < 0.01 is unadjusted; only three TC runs have terms
  under FDR control. Headline GO claims are safest from `consensus_terms_TC_<ont>.csv`.
- **Repository size.** `.git` is about 765 MB, mostly committed superseded tables and BLAST
  outputs; `vst_paired_<T>.csv` (7 to 11 MB each) is rewritten on every 09 run.
- **No continuous integration.** Nothing runs the pipeline on a push; the checks are the
  runners, the provenance files and the `09` checks.

## Up next

- Merge this pull request, then update the manuscript text.
- Library QC sensitivity run (T025G, T035G, T040F).
- Add an environment lockfile (renv or conda) for the R 4.4.3 / Bioconductor 3.20 set.
