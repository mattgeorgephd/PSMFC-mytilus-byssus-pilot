# PSMFC-mytilus-byssus-pilot

Byssal thread attachment of *Mytilus trossulus* under ocean acidification, warming and
hypoxia: tensometer pull tests before and after a 3-day exposure, foot and gill Tag-seq, and
the link between the two.

# Analysis folders

Each analysis folder (`02_` to `09_`) is self-contained with its own `.Rproj`, `01_code/`,
`02_data/`, `03_analyses/` and a README. Open the folder's own `.Rproj` (not the
repository-root one) before knitting, so `here::here()` resolves to that folder.

| folder | what it does | run |
|---|---|---|
| `00_experiment_plan/` | experimental design slides and photos | reference only |
| `00_treatment_conditions/` | tank DO, pH, temperature and salinity record and summary table | reference only |
| `01_mussel-measurements/` | mussel size, condition and thread-production workbooks | input: `mussel-size-measurements.xlsx` feeds `02_thread-strength` script 0 |
| `02_thread-strength/` | tensometer trace extraction, thread summary, per-animal ANCOVA on adhesion, force, plaque area and extension | knit `01_code/` scripts 0 to 4 in order |
| `03_blast/` | BLAST annotation of the genome CDS and the Iso-Seq transcriptome; the `genome-foot/` GO mapping used downstream | HPC method record; outputs committed |
| `04_iso-seq-transcriptome/` | QC of the Iso-Seq transcriptome (superseded isoseq-as-reference DE kept) | knit `01_code/05-IsoSeq-transcriptome-check.Rmd` |
| `05_sequence-alignment/` | read QC and HISAT2 + StringTie alignment to the genome | HPC method record; outputs committed |
| `06_differential-expression/` | DESeq2 per tissue and contrast, DEG annotation, DEG figures and tables | knit `01_code/00-run_differential_expression_pipeline.Rmd` |
| `07_enrichment/` | DAVID and REVIGO input lists and DAVID plots | knit `01_code/` scripts 07, 08, 09 |
| `08_gene-annotation/` | GO slims, NCBI summaries and orthologs for the top DEGs | knit `01_code/` scripts 06, 17, 18 (need network access) |
| `09_gene-mechanics-correlation/` | per-gene ANCOVA of day-3 thread mechanics on expression, foot and gill | knit `01_code/00-run_gene_mechanics_by_tissue.Rmd` |

Run order: `02_thread-strength` and `06_differential-expression` first (independent of each
other), then `07_enrichment` and `08_gene-annotation` (they read `06`'s DEG tables), and
`09_gene-mechanics-correlation` last (it reads `02`, `03` and `06`).

Other folders:

- `tools/pipeline_checks.R`: shared run checks and `RUN_provenance*.txt` writer, sourced by
  `02_thread-strength` scripts 3 and 4 and `09_gene-mechanics-correlation` scripts 20 and 21.
- `instrument-reference/`: tensometer manual, LabVIEW logger and wiring notes.
- `template-oyster-pipeline/`: Tag-seq code from the triploid oyster heatwave project, kept as
  a template; not part of this analysis.

## Cross-folder paths

Scripts refer to other analysis folders by name, so renaming a numbered folder breaks them.
The names are set in `06_differential-expression/01_code/_paths.R`,
`07_enrichment/01_code/_paths.R`, `08_gene-annotation/01_code/_paths.R`, the `paths` chunk of
each `09_gene-mechanics-correlation/01_code/2*.Rmd` script and its driver,
`02_thread-strength/01_code/0_build_mussel_key.Rmd`, and `PSMFC_ROOT_MARKERS` in
`tools/pipeline_checks.R`.

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
