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
  (`02`) and the Tag-seq counts (`05`) through differential expression (`06`), GO enrichment
  across treatments (`07`, `08`) and the per-animal associations of thread mechanics with
  genes, DEG sets, enriched GO terms and mitochondrial expression (`09`). The counts of record
  are featureCounts counts on the RefSeq annotation with Iso-Seq-extended 3' ends (made by
  `05_sequence-alignment` steps 04-06, taken by its step 07; adopted 2026-10-02). The
  Iso-Seq branch (`04` steps 01-03) repeats `06`'s contrasts on the Iso-Seq transcriptome as a
  sensitivity analysis, and `04` step 04 keeps the comparison the counts of record rest on.
- **Manuscript.** Linked from `README.md`. The contrasts of record are stressor vs day-3
  treatment control (TC). The mitochondrial genes have their own analysis and figure
  (`06` step 13, `figures/MT_mitochondrial_expression.png`).

## How to run

- **Everything:** knit `00_run_pipeline.Rmd` at the repository root (open
  `PSMFC-mytilus-byssus-pilot.Rproj`). About 35 minutes on four cores; on Windows the GO steps
  in `07` run on one core (R cannot fork there), so allow longer. It knits each folder's
  runner in a fresh R process and writes `knit_html/run_log.csv`; it fails if any stage fails.
- **One folder:** open the folder's own `.Rproj` and knit `01_code/00_run_*.Rmd`. Its
  `steps` parameter runs a subset; reports and logs go to `03_analyses/knit_html/`.
- **One script:** open the folder's `.Rproj` first, so `here::here()` resolves to the folder,
  then knit the script. Run the steps before it first.
- **Requirements:** the packages each folder's README lists (DESeq2, apeglm, ashr, topGO,
  goseq, clusterProfiler, enrichplot, rrvgo, GOSemSim, GO.db, GSEABase, org.Hs.eg.db,
  tidyverse, readxl, openxlsx, emmeans, here, rmarkdown, among others). Tested with R 4.4.3
  and Bioconductor 3.20 packages, with `GO.db` and `org.Hs.eg.db` 3.18.0 (GO release
  2023-07-27); the earlier analyses ran on R 4.2.2. GO results change with the GO release in
  `GO.db`, so `07` records it with the package versions in
  `03_analyses/01_go-inputs/RUN_provenance.txt`.
- **Not run by the pipeline:** `03_blast/` and `05_sequence-alignment` step 01 (HPC, inputs
  not in the repository), `05` steps 04-06 and step 08's mitochondrial alignment (download the
  transcriptome, genome and reads and run minimap2, HISAT2, StringTie and featureCounts; set
  the `05` runner's `online: true`; steps 04-06 make the count matrix of record, which step 07
  takes), `05` step 03's recipe check (its own `online: true`), `04_iso-seq-transcriptome`
  steps 01-02 (salmon; set `online: true`; steps 03 and 04 run as pipeline stage 04) and
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
    `06_differential-expression/03_analyses/_superseded/LC_contrasts/`. The foot-vs-gill
    contrast in the day-3 controls (`FG_TC`) is computed and drawn, but nothing of record
    rests on it.
  - Counts: featureCounts (Subread 2.1.1), uniquely aligned reads on the sense strand, on the
    RefSeq annotation (RS_2024_02) with each transcript's last exon extended to the 3' end of
    its Iso-Seq isoforms (`ext3`), from HISAT2 2.2.1 alignments
    (`05_sequence-alignment/03_analyses/featurecounts/gene_count_matrix.csv`). The previous
    counts (StringTie + prepDE on the HPC alignment, `05 .../prepDE/`) are kept for comparison
    (`04` step 04).
  - DEGs: apeglm-shrunk DESeq2 (Wald p from DESeq2's outlier-replaced refit, the standard
    workflow), padj < 0.05, after keeping genes with at least 10 counts in a third of the
    contrast's samples; one model per contrast, defined in
    `06_differential-expression/01_code/02_define_contrasts.Rmd`. The 331 mitochondrial loci
    (`count_matrix/mitochondrial_loci.csv`: the mitogenome's 12 protein genes, 2 rRNAs and 24
    tRNAs, and 293 nuclear-scaffold copies of the protein genes) are left out of these fits and
    of 07 to 09's gene universes.
  - Mitochondrial proteins: counted on the mitochondrial genome alone (`05_sequence-alignment`
    step 08: HISAT2 default scoring against NC_007687.1, featureCounts per gene) and tested per
    protein in `06` step 13 with the nuclear genes' size factors; the haplotype covariate and
    a permissive alignment score are sensitivity checks there.
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
  were removed at QC; T047 has no foot library. `library_crosswalk.csv` (06) maps every library
  to its RNA isolation record.
- **Known data issues** (see the folder READMEs): 293 loci on unplaced scaffolds (126
  protein-coding LOCs and 167 pseudogenes) are copies of the mitochondrial protein genes and
  take their reads in the genome alignment (`tools/mt_encoded.R`); 12 of the 59 animals carry
  mitochondrial haplotypes that differ from the reference at fixed positions, which the genome
  alignment undercounts (both handled by counting on the mitogenome alone, as above;
  `06` step 13); byssal plaque genes are expressed in the day-0 libraries and mostly
  absent at day 3, which the different day-0 dissection may explain; the analysed reads are
  leader-clipped and PCR-deduplicated Tag-seq reads (`05` README, "The reads"), whose raw
  reads are on owl (`nightingales/M_trossulus/`); the Iso-Seq branch (`04`, a sensitivity
  analysis) replicates most genome DEGs (58 to 71%) and finds 1.2 to 1.6 times as many; against
  the previous StringTie + prepDE counts it had found about twice as many, mostly because of
  the counter (prepDE's read-base estimates), less because the RefSeq gene models end before
  many Tag-seq 3' ends (`04` README, option B).
- **Secrets.** API keys go in environment variables (`ENTREZ_KEY` in `~/.Renviron`), never in
  a file in the repository.
- **Writing.** READMEs and comments in plain language, with commas or semicolons rather than
  dashes as punctuation; cite primary sources for methods (paper and version).
- **Git.** Work on a branch and open a pull request; do not rewrite history on shared
  branches. Generated outputs are committed so results can be read without running anything;
  regenerate them with the runner rather than editing them.

## What not to touch

- **Raw data:** `01_mussel-measurements/*.xlsx`, `02_thread-strength/02_data/tensometer_output/`,
  `06_differential-expression/02_data/` sample sheets, `00_experiment_plan/`,
  `00_treatment_conditions/`, `instrument-reference/`.
- `02_thread-strength/02_data/pad_area_measurements.xlsx` except to enter new plaque
  measurements by hand (it is the only hand-edited input).
- **Outputs that cannot be rebuilt here:** `03_blast/03_analyses/`,
  `05_sequence-alignment/03_analyses/hisat/`, `fastqc/` and
  `prepDE/transcript_count_matrix.csv` (HPC), and `08_gene-annotation/03_analyses/Top_gene_summaries/`
  (network). Replace them only by rerunning their scripts where the inputs exist.
- **`_superseded/` folders:** records; do not edit, run or delete them.
- **`template-oyster-pipeline/`:** reference code from another project, not part of this
  analysis.
- **Git-ignored files** (`knit_html/`, `dds/*.rds`, BAM/SAM files, the genome and the Iso-Seq
  FASTA): never force-add them.
- **Generated tables and figures:** never hand-edit; change the code and rerun.
