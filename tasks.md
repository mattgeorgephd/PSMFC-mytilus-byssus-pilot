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
- **Mitochondrial genes on their own.** The 143 mitochondrial loci (12 protein and 5 RNA genes
  of the mitochondrial genome, 126 nuclear-scaffold copies) are out of the genome DE, GO and
  gene-mechanics universes and tested per protein in `06` step 13, with a manuscript figure
  (`figures/MT_mitochondrial_expression.png`): Gill OA raises 8 of 12 proteins, Foot OA 5,
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

## In progress

- Review and merge of the pull request carrying this work (branch `claude/jolly-clarke-7afl1v`).
- Proofreading the regenerated figures before they go into the manuscript.

## Blocked or waiting on a decision

| item | why it is blocked | who |
|---|---|---|
| Revoke the old NCBI API key | it was removed from the code but remains in the public history (commit 3bee414); only its owner can revoke it at NCBI | key owner |
| Merge the pull request | branch `claude/jolly-clarke-7afl1v`; merging deletes untracked `.Rhistory` copies from local checkouts (see the PR) | Matt |
| Iso-Seq branch | inputs located (2026-10-01): the trimmed reads HISAT2 used are in Grace's `byssus-exp-analysis/data/raw-trimmed/` on owl, and the transcriptome on owl's genomic-databank; but owl, gannet, NCBI and sr320.github.io (Steven Roberts' notebooks) are not reachable from the cloud environment used here, and no isoform-to-gene table exists yet. Design in `04_iso-seq-transcriptome/README.md` | Matt (network access or copies), Steven (notebook details) |
| Manuscript text | numbers and wording predate this pull request: Gill OA DEGs (711 to 543), the GO results, LC contrasts, the mitochondrial result, thread metric names (`mean_force`, `max_force`), extension | Matt |
| Manuscript GO figure | every option is drawn (`07_enrichment/03_analyses/02_topgo/` to `06_method-comparison/`); choose one | Matt |
| Library QC | three libraries carry most of DESeq2's outlier replacements: T040F (91 of 122 genes refit in Foot OW), T025G (48 of 76 in Gill OA; also the lowest alignment rate, 60.8% against a median of 67.7%) and T035G (42 of 81 in Gill OW). Check their RNA records; a sensitivity run without them would show whether any result rests on them | Matt |
| Exploratory foot-region contrast | rest of foot vs phenol gland to tip in the 12 day-0 animals (3,174 of 7,393 genes differ); add as a DE family if useful for interpreting byssal genes | Matt |
| `08` steps 02-03 (NCBI summaries, OrthoDB orthologs) | need network access; the committed summaries predate the current top-50 lists (and the mitochondrial separation) | anyone with network |
| `03_blast/` and `05` step 01 | HPC only (inputs and binaries not in the repository); committed outputs are used. `01_genome_blast.Rmd` would not run as written (`mytfoot_uniprot_info` read but `mytilus_foot_uniprot_info` used; `masterID` never defined) | HPC user |

## Known limitations (documented, not blocking)

- **No locked software environment.** The pipeline was tested with R 4.4.3 and Bioconductor
  3.20 (versions recorded in each `RUN_provenance*.txt`), but the repository has no `renv.lock`
  or conda environment file, so a new machine has to assemble the packages by hand.
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
