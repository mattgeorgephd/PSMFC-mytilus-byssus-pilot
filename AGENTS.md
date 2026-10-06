# AGENTS.md

Instructions for anyone, person or coding agent, working in this repository: what the project
is for, how it is organised, how to run it, and what not to touch. `README.md` describes the
analysis folders; `tasks.md` tracks what is done, in progress and blocked.

## Project goals

- **Question.** How do ocean acidification (OA), ocean warming (OW) and hypoxia (DO) change
  byssal thread attachment in *Mytilus trossulus*, and which genes in the foot (and gill) track
  those changes? Mussels were pulled on a tensometer before and after a 3-day exposure and
  sampled for Tag-seq.
- **Analysis.** One reproducible chain, runnable from the committed data: thread mechanics
  (`02`) and the Tag-seq counts (`04`) through differential expression (`05`), GO enrichment
  across treatments (`07`, `08`) and the per-animal associations of thread mechanics with
  genes, DEG sets, enriched GO terms and mitochondrial expression (`09`). The counts of record
  are featureCounts counts on the RefSeq annotation with Iso-Seq-extended 3' ends (made by
  `04_sequence-alignment` steps 04-06, taken by its step 07; adopted 2026-10-02). The
  Iso-Seq branch (`06` steps 01-03) repeats `05`'s contrasts on the Iso-Seq transcriptome as a
  sensitivity analysis, and `06` step 04 keeps the comparison the counts of record rest on.
- **Manuscript.** Linked from `README.md`. The contrasts of record are stressor vs day-3
  treatment control (TC). The mitochondrial genes have their own analysis and figure
  (`05` step 13, `figures/MT_mitochondrial_expression.png`).

## How to run

- **Everything:** knit `00_run_pipeline.Rmd` at the repository root (open
  `PSMFC-mytilus-byssus-pilot.Rproj`). About 55 minutes on four cores (20 of them `09` step 06's prediction test); on Windows the GO steps
  in `07` run on one core (R cannot fork there), so allow longer. It knits each folder's
  runner in a fresh R process and writes `knit_html/run_log.csv`; it fails if any stage fails,
  and its error quotes the failed step's own error from that step's log. Before the first
  stage it compares this R with `renv.lock` (`check_stack()`, `tools/run_steps.R`) and stops
  within seconds if R's minor version or Bioconductor's differs, if a package the pipeline
  loads (`PIPELINE_PACKAGES`) is not on R's library path (put the `renv::restore()` library
  there with `R_LIBS=<library>` in `.Renviron`), or, when `07` or `08` is to run, if `GO.db`
  does not hold GO release 2026-01-23; its `check_versions: false` skips that check (outputs
  are then not of record).
- **One folder:** open the folder's own `.Rproj` and knit `01_code/00_run_*.Rmd`. Its
  `steps` parameter runs a subset; reports and logs go to `03_analyses/knit_html/`.
- **One script:** open the folder's `.Rproj` first, so `here::here()` resolves to the folder,
  then knit the script. Run the steps before it first.
- **Requirements:** the packages each folder's README lists (DESeq2, apeglm, ashr, topGO,
  goseq, clusterProfiler, enrichplot, rrvgo, GOSemSim, GO.db, GSEABase, org.Hs.eg.db,
  tidyverse, readxl, openxlsx, emmeans, here, rmarkdown, R.utils (for `data.table::fread` on
  `.gz` files), glmnet and dynamicTreeCut (`09` step 06), among others). Tested with R 4.6.1 and Bioconductor 3.23
  (`BiocManager::install(version = "3.23")`): DESeq2 1.52.0, apeglm 1.34.0, ashr 2.2.63, topGO
  2.64.0, goseq 1.64.0, clusterProfiler 4.20.0, rrvgo 1.24.0, GOSemSim 2.38.3, `GO.db` 3.23.1
  (GO release 2026-01-23) and `org.Hs.eg.db` 3.23.1. The previous runs used R 4.4.3 with
  Bioconductor 3.20 and, before that, R 4.2.2; on the same inputs R 4.6.1 / Bioconductor 3.23
  gave the same DEGs, GO terms, GO slims and gene-mechanics results (log2 fold changes within
  1e-6). GO results change with the GO release in `GO.db`, so `07` and `08` stop unless it is
  2026-01-23 (`check_go_release()`, `tools/pipeline_checks.R`), and `07` records it with the
  package versions in `03_analyses/01_go-inputs/RUN_provenance.txt`. `renv.lock` at the root
  records every package of the library the pipeline last ran with (282, R's recommended packages
  included: Matrix 1.7-6, MASS 7.3-66, mgcv 1.9-4, survival 3.8-12, ...): with R 4.6.1,
  `install.packages("renv"); renv::restore(lockfile = "renv.lock", library = "<library>",
  prompt = FALSE)` rebuilds it (the project does not activate renv, so nothing else changes).
  System libraries are not in it (on Ubuntu 24.04: libcurl, libssl, libxml2, libfontconfig,
  libharfbuzz, libfribidi, libfreetype, libpng, libtiff, libjpeg, libwebp, libcairo2, libglpk,
  libgmp, libicu, libuv, libnlopt, libgit2, ImageMagick and pandoc, as `-dev` packages). On a
  machine that had an older R, keep the old R's packages off the library path (on Ubuntu,
  packages in `/usr/lib/R/site-library` built for R 4.3 fail to load in R 4.6), for example
  `R_LIBS_SITE=/nonexistent R_LIBS=<library>` when installing and running, and upgrade R's
  recommended packages with R (CRAN's Ubuntu repository builds them for 4.6).
  `tools/cloud_setup.sh` does all of this on Ubuntu 24.04 (R 4.6.1, the recommended packages,
  the system libraries, BLAST+ 2.15.0 and `renv::restore()` into `/opt/R/site-library-4.6`); it
  is the setup script of the Claude Code cloud environment (environment settings, Setup
  script).
- **Not run by the pipeline:** `03_blast/` steps 01-03 and `04_sequence-alignment` step 01
  (HPC or long searches, inputs not in the repository; `03_blast` step 01 ran its 2026 search in
  the cloud environment in about 31 hours on 4 threads, in resumable parts), `03_blast` step 04
  (the 2024 hits with UniProt 2026_03 records, a comparison; it runs offline from committed
  files, and `online: true` fetches the records again, which works only while UniProt serves
  release 2026_03), `03_blast` step 05 (the foot and byssal proteins searched again with the
  low-complexity filter off; `run: true` needs BLAST+ and the genome CDS and takes about 2
  minutes, `run: false` rebuilds its tables offline), `04` steps 04-06 and step 08's mitochondrial alignment (download the
  transcriptome, genome and reads and run minimap2, HISAT2, StringTie and featureCounts; set
  the `04` runner's `online: true`; steps 04-06 make the count matrix of record, which step 07
  takes), `04` step 03's recipe check (its own `online: true`), `06_iso-seq-transcriptome`
  steps 01-02 (salmon; set `online: true`; steps 03 and 04 run as pipeline stage 06) and
  `08_gene-annotation` steps 02-03 (NCBI and OrthoDB; set `online: true`). Their committed
  outputs are what the pipeline reads.
- **After a run:** check that every `run_log.csv` row is `TRUE`; that the `RUN_provenance*.txt`
  files in `02`, `07` and `09` name the commit you ran and say `uncommitted changes : no`; and
  that every `checks` line in `09`'s says `0 failed`.

## Conventions

- **Folders.** `00_`/`01_` hold reference material and raw workbooks; `02_` to `09_` are
  analyses, each with its own `.Rproj`, `01_code/`, `02_data/`, `03_analyses/` and `README.md`.
- **Scripts.** In `01_code/`, `00_run_*.Rmd` is the batch runner and `01_`, `02_`, ... are the
  steps in run order. `_*.R` files are helpers the steps source. Retired scripts live in
  `_superseded/`.
- **Writes.** A script writes only to its own folder's `03_analyses/` (often a subfolder named
  after the step). Later folders read earlier folders' outputs and never write into them.
  `02_data/` holds inputs only; the one exception is a script downloading an external input,
  which is git-ignored.
- **Gene names and keys.** Join count-matrix genes to annotation only through `gene_key()`
  (`tools/gene_ids.R`); never rebuild the key with an ad hoc `sub()`. One BLAST hit per gene:
  the highest bitscore (ties in table order). Select and rename columns by name, never by
  position, and check a join kept one row per gene before binding columns to it.
- **Paths.** `here::here()` anchored on the folder's `.Rproj`, and `repo_root` (its parent)
  for other folders. Never absolute paths. Folder names are hard-coded in each `_paths.R`,
  the `paths` chunks of `09`, `00_run_pipeline.Rmd` and `tools/pipeline_checks.R`, so renaming
  a numbered folder means updating them.
- **Step map.** The root `README.md` ("How the steps connect") maps what every step reads and
  writes. Update its table and diagram in the same change as any change to a step's inputs or
  outputs.
- **Nothing is deleted.** Retire a script or output to a `_superseded/` folder with a README
  saying what replaced it. Every analysis folder (which lists its scripts) and its `02_data/`
  and `03_analyses/` have a README listing the files and what produces them; output subfolders
  have their own or are listed in their parent's. Update them in the same change as the code.
- **Figures.** Colours come only from `tools/plot_style.R`: control grey, OA green, OW orange,
  DO purple (thread strength adds baseline blue and lab reference light grey); red up, blue
  down; neutral greys for significance scales and set categories. The palette is checked for
  colour-vision deficiency; keep a legend or labels on every figure.
- **Statistics of record.**
  - Control: the day-3 treatment control (TC), the only one. The day-0 "lab control"
    animals are not a control anywhere: their feet were dissected differently (two pieces at
    day 0, one at day 3). The LC contrasts are retired to
    `05_differential-expression/03_analyses/_superseded/LC_contrasts/`. The foot-vs-gill
    contrast in the day-3 controls (`FG_TC`) is computed and drawn, but nothing of record
    rests on it.
  - Counts: featureCounts (Subread 2.1.1), uniquely aligned reads on the sense strand, on the
    RefSeq annotation (RS_2024_02) with each transcript's last exon extended to the 3' end of
    its Iso-Seq isoforms (`ext3`), from HISAT2 2.2.1 alignments
    (`04_sequence-alignment/03_analyses/featurecounts/gene_count_matrix.csv`). The previous
    counts (StringTie + prepDE on the HPC alignment, `04 .../prepDE/`) are kept for comparison
    (`06` step 04).
  - DEGs: apeglm-shrunk DESeq2 (Wald p from DESeq2's outlier-replaced refit, the standard
    workflow), padj < 0.05, after keeping genes with at least 10 counts in a third of the
    contrast's samples; one model per contrast, defined in
    `05_differential-expression/01_code/02_define_contrasts.Rmd`. The 331 mitochondrial loci
    (`count_matrix/mitochondrial_loci.csv`: the mitogenome's 12 protein genes, 2 rRNAs and 24
    tRNAs, and 293 nuclear-scaffold copies of the protein genes) are left out of these fits and
    of 07 to 09's gene universes.
  - Mitochondrial proteins: counted on the mitochondrial genome alone (`04_sequence-alignment`
    step 08: HISAT2 default scoring against NC_007687.1, featureCounts per gene) and tested per
    protein in `05` step 13 with the nuclear genes' size factors; the haplotype covariate and
    a permissive alignment score are sensitivity checks there.
  - Annotation: each gene's best hit (highest bitscore) in the genome blastx of 2026
    (Swiss-Prot release 2026_03 plus the 196 proteins of the UniProt query "(mytilus foot)"
    and 45 byssal proteins it misses, `03_blast/02_data/byssal_additions_2026_03.tsv`;
    `03_blast` step 01), with the UniProt records of release 2026_03; and, for genes with no hit
    there, a hit of at least 70% identity to those 241 foot and byssal proteins with blastx's
    low-complexity filter off (`03_blast` step 05: mfp-3, mfp-5, foot protein 9 and thread
    matrix protein 2F, six genes the filter had hidden)
    (`03_blast/03_analyses/genome-foot-sprot2026_03-noseg/LOC_GO_list.txt`); GO release 2026-01-23
    (`GO.db` 3.23.1) and its generic GO slim (`08_gene-annotation/02_data/goslim_generic.obo`).
    Kept for comparison: the search of 2024 (`genome-foot/`) and its hits with the 2026_03
    records (`genome-foot-uniprot2026_03/`, `03_blast` step 04).
  - Expression suites (`09` step 06, exploratory): script 05's gene sets, co-expression
    modules and expression components, each through the within-arm ANCOVA and the same model
    without treatment; elastic-net prediction of held-out animals (glmnet 5.1) against label
    permutations, seeded per run so that the result does not depend on the number of workers,
    with a power check (the treatment's columns added to all the axes, across arms).
  - GO enrichment: each contrast's tested genes (non-missing padj) are its universe; up- and
    down-regulated genes are tested separately; topGO `weight01` p < 0.01 is of record, goseq
    and clusterProfiler (BH < 0.05) are comparisons.
  - Thread strength: per-animal ANCOVA, day-3 level adjusted for the animal's own baseline,
    each arm against control, on the log scale. `max_force` in the thread table is one
    thread's peak force; per animal, `mean_force` is the mean of its threads' peak forces
    (geometric mean in the model; primary) and `max_force` its strongest thread
    (exploratory). Extension and the area under the force-time curve are not analysed:
    threads were cut near the plaque, so the length of distal thread under test was not
    standardized. Peak force is the only quantity taken from a trace.
- **Samples.** `F` libraries are the phenol gland to the tip of the foot (every animal); `FX`
  are the rest of the foot (day-0 animals only); `G` are gill. No contrast uses a day-0
  library. T051F and T051G
  were removed at QC; T047 has no foot library. `library_crosswalk.csv` (05) maps every library
  to its RNA isolation record.
- **Known data issues** (see the folder READMEs): blastx's low-complexity filter (SEG) hid the
  hits of six byssal genes in the 2026 search (three mfp-3, mfp-5, foot protein 9, thread
  matrix protein 2F), which `03_blast` step 05 adds; genes whose hit of record is a weaker
  repeat match keep it, and two genes RefSeq names "adhesive plaque matrix protein-like"
  (LOC134723087, LOC134723088, 36 to 37% to mfp-1) have no annotation; 293 loci on unplaced scaffolds (126
  protein-coding LOCs and 167 pseudogenes) are copies of the mitochondrial protein genes and
  take their reads in the genome alignment (`tools/mt_encoded.R`); 12 of the 59 animals carry
  mitochondrial haplotypes that differ from the reference at fixed positions, which the genome
  alignment undercounts (both handled by counting on the mitogenome alone, as above;
  `05` step 13); byssal plaque genes are expressed in the day-0 libraries and mostly
  absent at day 3, which the different day-0 dissection may explain; the analysed reads are
  leader-clipped and PCR-deduplicated Tag-seq reads (`04` README, "The reads"), whose raw
  reads are on owl (`nightingales/M_trossulus/`); the Iso-Seq branch (`06`, a sensitivity
  analysis) replicates most genome DEGs (58 to 71%) and finds 1.2 to 1.6 times as many; against
  the previous StringTie + prepDE counts it had found about twice as many, mostly because of
  the counter (prepDE's read-base estimates), less because the RefSeq gene models end before
  many Tag-seq 3' ends (`06` README, option B).
- **Secrets.** API keys go in environment variables (`ENTREZ_KEY` in `~/.Renviron`), never in
  a file in the repository.
- **Writing.** READMEs and comments in plain language, with commas or semicolons rather than
  dashes as punctuation; cite primary sources for methods (paper and version).
- **Git.** Work on a branch and open a pull request; do not rewrite history on shared
  branches. Generated outputs are committed so results can be read without running anything;
  regenerate them with the runner rather than editing them.

## What not to touch

- **Raw data:** `01_mussel-measurements/*.xlsx`, `02_thread-strength/02_data/tensometer_output/`,
  `05_differential-expression/02_data/` sample sheets, `00_experiment_plan/`,
  `00_treatment_conditions/`, `instrument-reference/`.
- `02_thread-strength/02_data/pad_area_measurements.xlsx` except to enter new plaque
  measurements by hand (it is the only hand-edited input).
- **Outputs that cannot be rebuilt here:** `03_blast/03_analyses/`,
  `04_sequence-alignment/03_analyses/hisat/`, `fastqc/` and
  `prepDE/transcript_count_matrix.csv` (HPC), and `08_gene-annotation/03_analyses/Top_gene_summaries/`
  (network). Replace them only by rerunning their scripts where the inputs exist.
- **`_superseded/` folders:** records; do not edit, run or delete them.
- **`template-oyster-pipeline/`:** reference code from another project, not part of this
  analysis.
- **Git-ignored files** (`knit_html/`, `dds/*.rds`, BAM/SAM files, the genome, the Iso-Seq
  FASTA and the tensometer QC plots in `02_thread-strength/03_analyses/02_extract-tensometer-data/QC_plots/`):
  never force-add them.
- **Generated tables and figures:** never hand-edit; change the code and rerun.
