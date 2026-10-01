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
- **06 differential expression.** The 16 contrasts defined by rule in one script (they
  reproduce the old sample lists exactly); one saved DESeq2 fit per contrast, so shrinkage runs
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
- **Documentation.** A README in every new folder; `AGENTS.md` and this file.

## In progress

- Review and merge of the pull request carrying this work (branch `claude/jolly-clarke-7afl1v`).
- Proofreading the regenerated figures (new colours) before they go into the manuscript.

## Blocked or waiting on a decision

| item | why it is blocked | who |
|---|---|---|
| Revoke the old NCBI API key | it was removed from the code but remains in the public history (commit 3bee414); only its owner can revoke it at NCBI | key owner |
| Download large inputs from gannet in the scripts (curl) | needs the list of gannet file paths each script should fetch; planned as a follow-up pull request | Matt |
| Mitochondrially encoded protein LOCs | about 140 LOCs (116 of Gill OA's 711 TC DEGs) carry one mitochondrial signal; decide whether to keep them (now flagged), collapse them to one gene per protein, or exclude them from DEG counts and GO | Matt |
| Interpreting LC contrasts | byssal plaque genes are expressed in day-0 foot and mostly absent at day 3, controls included, so LC contrasts mix byssal secretion state into the stressor effect | Matt |
| Exploratory foot-region contrast | rest of foot vs phenol gland to tip in the 12 day-0 animals (3,174 of 7,393 genes differ); add as a DE family if useful for interpreting byssal genes | Matt |
| Iso-Seq branch | design in `04_iso-seq-transcriptome/README.md`; not implemented until it is wanted | Matt |
| `08` steps 02-03 (NCBI summaries, OrthoDB orthologs) | need network access; not rerun since the restructure | anyone with network |
| `03_blast/` and `05` step 01 | HPC only (inputs and binaries not in the repository); committed outputs are used | HPC user |

## Up next

- Decide the three analysis questions above (mitochondrial LOCs, LC interpretation, region
  contrast), then rerun `00_run_pipeline.Rmd`.
- Gannet download follow-up pull request.
- Choose the manuscript GO figure: the topGO TC dot plot, or the rrvgo parent-term view of the
  same terms (`07_enrichment/03_analyses/02_topgo/`, `05_rrvgo/`).
- Optional figure polish, left as drawn: the failure-mode figure (`02` step 05) colours
  "peeling" orange, close to the OW orange; the two per-mussel-mean box plots (`02` step 04)
  clip their titles; `panel.spacing` in the DEG-count figure (`06` step 10) is set before
  `theme_classic()` and so has no effect.
