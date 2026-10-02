# PSMFC-mytilus-byssus-pilot

Byssal thread attachment of *Mytilus trossulus* under ocean acidification, warming and
hypoxia: tensometer pull tests before and after a 3-day exposure, foot and gill Tag-seq, and
the link between the two.

# How to run

Knit `00_run_pipeline.Rmd` at the repository root (inside `PSMFC-mytilus-byssus-pilot.Rproj`).
It runs, in order and each in a fresh R process, the batch runner of every folder that can run
from the committed data: thread strength (02), the count matrices (05), differential expression
(06), GO enrichment (07), GO slims (08) and the gene-mechanics associations (09). About 25
minutes. Reports and logs go to each folder's `03_analyses/knit_html/` and, one per stage, to
`knit_html/` at the root (all git-ignored). See `AGENTS.md` for the conventions and `tasks.md`
for what is done and open.

# Analysis folders

Each analysis folder (`02_` to `09_`) is self-contained with its own `.Rproj`, `01_code/`,
`02_data/`, `03_analyses/` and a README. Open the folder's own `.Rproj` (not the
repository-root one) before knitting a single script, so `here::here()` resolves to that
folder. In every `01_code/`, `00_run_*.Rmd` is the folder's batch runner and `01_` onwards are
the steps in run order; `_*.R` files are helpers they source. Each script writes only to its
own folder's `03_analyses/`; later folders read earlier ones.

| folder | what it does | run |
|---|---|---|
| `00_experiment_plan/` | experimental design slides and photos | reference only |
| `00_treatment_conditions/` | tank DO, pH, temperature and salinity record and summary table | reference only |
| `01_mussel-measurements/` | mussel size, condition and thread-production workbooks | input: `mussel-size-measurements.xlsx` feeds `02_thread-strength` script 01 |
| `02_thread-strength/` | tensometer trace extraction, thread summary, per-animal ANCOVA on adhesion, mean and maximum peak force and plaque area | `01_code/00_run_thread_strength.Rmd` |
| `03_blast/` | BLAST annotation of the genome CDS and the Iso-Seq transcriptome; the `genome-foot/` GO mapping used downstream | HPC method record; outputs committed |
| `04_iso-seq-transcriptome/` | sensitivity branch: the TC contrasts repeated with the reads quantified against the Iso-Seq transcriptome (isoforms mapped to genome genes, salmon, tximport) and compared with `06` | runner `01_code/00_run_isoseq.Rmd` (step 04 by default; steps 02-03 with `online: true`) |
| `05_sequence-alignment/` | read QC, HISAT2 + StringTie alignment (HPC record) and the count matrices | `01_code/00_run_sequence_alignment.Rmd` |
| `06_differential-expression/` | DESeq2 for 7 contrasts (each stressor vs the day-3 treatment control, and foot vs gill), DEG annotation, figures; the mitochondrial genes on their own | `01_code/00_run_differential_expression.Rmd` |
| `07_enrichment/` | GO enrichment: topGO (of record), goseq, clusterProfiler, rrvgo, method comparison | `01_code/00_run_enrichment.Rmd` |
| `08_gene-annotation/` | GO slims of the TC DEGs; NCBI summaries and orthologs for the top DEGs (network) | `01_code/00_run_gene_annotation.Rmd` |
| `09_gene-mechanics-correlation/` | per-animal ANCOVA of day-3 thread mechanics on genes, DEG sets, enriched GO terms and mitochondrial expression, foot and gill | `01_code/00_run_gene_mechanics_by_tissue.Rmd` |

Run order: `02_thread-strength` and `05` -> `06` (independent of each other), then
`07_enrichment` and `08_gene-annotation` (they read `06`), and
`09_gene-mechanics-correlation` last (it reads `02`, `03`, `06` and `07`). `00_run_pipeline.Rmd`
follows this order.

Other folders:

- `tools/`: shared helpers (README inside): `run_steps.R` (the runners), `plot_style.R` (the
  one set of figure colours: control grey, OA green, OW orange, DO purple; red up, blue
  down), `gene_ids.R` (`gene_key()`, the one way gene names are joined to annotation),
  `mt_encoded.R` (the mitochondrial loci of the count matrix) and `pipeline_checks.R` (run
  checks and `RUN_provenance*.txt`).
- `instrument-reference/`: tensometer manual, LabVIEW logger and wiring notes.
- `template-oyster-pipeline/`: Tag-seq code from the triploid oyster heatwave project, kept as
  a template; not part of this analysis.

## Samples

Tissue was foot or gill. Every animal has a library of the phenol gland to the tip of the foot
(IDs ending `F`) and of the gill (`G`); the twelve day-0 animals also have a library of the
rest of the foot (`FX`). The day-0 animals are not used as a control (their feet were
dissected differently); every contrast is against the day-3 treatment control. The sample sheets name these
inconsistently; `06_differential-expression/03_analyses/count_matrix/library_crosswalk.csv`
maps every library to its RNA isolation record.

## Cross-folder paths

Scripts refer to other analysis folders by name, so renaming a numbered folder breaks them.
The names are set in each folder's `01_code/_paths.R` (06, 07, 08), the `paths` chunk of each
`09_gene-mechanics-correlation/01_code/0*.Rmd` script and its runner,
`02_thread-strength/01_code/01_build_mussel_key.Rmd`, the stage table of `00_run_pipeline.Rmd`
and `psmfc_repo_root()` in `tools/pipeline_checks.R`.

## Large files

Files too large for GitHub, such as the raw Tag-seq reads
(`20220405-tagseq/`), are stored on gannet:
https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/


# Pertinent documents
## General
1. [Manuscript](https://docs.google.com/document/d/1fKfDU4gHPdMy9xejUA5pZ6YY2HVzsiGSGo7uFVgDln4/edit?usp=sharing)

## Manuals and Protocols
1. [Thread testing tutorial](https://monicaklopp.github.io/Thread-Testing-01-Notebook-Post/)
2. [Sam's RNA extraction notebook entries](https://robertslab.github.io/sams-notebook/2022/01/13/Project-Summary-Matt-George-PSMFC-Mytilus-Byssus-Project.html)
3. [RNA sample List](https://docs.google.com/spreadsheets/d/1PDVSGuCGeYQr6Rdl6u5M4L5vcQS1EgUQl7UjvLYDlBg/edit?usp=sharing)

## Datasets
1. [Tagseq dataset](https://docs.google.com/spreadsheets/d/1zZ6L05j-SyYJbzzQI_kBafFaReE4Ysp_9bORdbBu_r8/edit#gid=1302342348)
2. [RNA extraction results](https://docs.google.com/spreadsheets/d/1HizNOIfhSjppHDQrWLGiJhuDZKO8c-qm9JAz0Z5QIIQ/edit?usp=sharing)

## Github Issues
7. [RNA extraction github issue](https://github.com/RobertsLab/resources/issues/1352)
14. [Iso-seq analysis github issue](https://github.com/RobertsLab/resources/issues/1662)
