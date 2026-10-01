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

## In progress

- Review and merge of the pull request carrying this work (branch `claude/jolly-clarke-7afl1v`).
- Proofreading the regenerated figures before they go into the manuscript.

## Blocked or waiting on a decision

| item | why it is blocked | who |
|---|---|---|
| Revoke the old NCBI API key | it was removed from the code but remains in the public history (commit 3bee414); only its owner can revoke it at NCBI | key owner |
| Iso-Seq branch | the trimmed reads (lab server) and the trimming command are not in the repository; owl, gannet and NCBI are not reachable from the cloud environment used here; no isoform-to-gene table exists yet. Design in `04_iso-seq-transcriptome/README.md` | Matt (inputs), anyone with access |
| Download large inputs from gannet in the scripts (curl) | needs the list of gannet file paths each script should fetch; planned as a follow-up pull request | Matt |
| Manuscript GO figure | every option is drawn (`07_enrichment/03_analyses/02_topgo/` to `06_method-comparison/`); choose one | Matt |
| Exploratory foot-region contrast | rest of foot vs phenol gland to tip in the 12 day-0 animals (3,174 of 7,393 genes differ); add as a DE family if useful for interpreting byssal genes | Matt |
| Outlier-replaced genes | DESeq2 replaces extreme counts and refits (Wald p) while apeglm fold changes use the original counts; 13 TC DEGs are affected, 1 by more than 25% (06 README). Kept as the DESeq2 convention; refit apeglm on the replaced counts if preferred | Matt |
| `08` steps 02-03 (NCBI summaries, OrthoDB orthologs) | need network access; not rerun since the restructure, and the FOA and GOA tables predate the current top-50 lists | anyone with network |
| `03_blast/` and `05` step 01 | HPC only (inputs and binaries not in the repository); committed outputs are used | HPC user |

## Up next

- Choose the manuscript GO figure.
- Gannet download follow-up pull request; then the Iso-Seq branch once its inputs are reachable.
