# Phase 4: differential-expression rewrite

Makes the absorbed DESeq2/DEG scripts runnable from the committed CSVs by repointing every
path to `here()` / `02_data` / `03_analyses`, with cross-folder reads through a `repo_root`
pointer. Delivered as a single apply script (not 14 file copies) so nothing is missed.

## Apply

```
git switch -c reorg/phase4-de
bash apply_phase4_de.sh          # from repo root (uses Perl, bundled with Git; no Python needed)
git status && git diff
git push -u origin reorg/phase4-de
```

The script: creates `differential-expression/01_code/_paths.R`; repoints the 14 genome-based
scripts in place; and `git mv`s `01-diff-exp-analysis.Rmd` into `01_code/_superseded/`. It
refuses to run twice (it aborts if `_paths.R` already exists).

## What changed

- **Shared `_paths.R`** defines `repo_root`, `dat` (= `02_data`), `deg` (= `03_analyses/DEG_lists`),
  and cross-folder `blast_go` and `topgenes`. Every script sources it via a small chunk inserted
  after the YAML. Reuse this same pattern for the other phase-4 folders.
- **Path mapping applied to all 14 scripts:**
  - `output/{gene_count_matrix_clean,gene_count_matrix,transcript_count_matrix,treatmentinfo_clean}.csv`
    and `data/PSMFC-...RNA-tagseq_raw.csv` -> `file.path(dat, ...)`
  - `output/DEG_lists/...` -> `file.path(deg, ...)`
  - `output/LOC_GO_list.txt`, `output/g.spid.txt` -> `file.path(blast_go, ...)` (cross-folder, from blast)
  - `output/Top_gene_summaries/...` -> `file.path(topgenes, ...)` (cross-folder, gene-annotation)
  - Grace's absolute `/home/shared/.../byssus-exp-analysis/` prefixes stripped first.

## Decision flagged: 01-diff-exp-analysis -> _superseded

`01-diff-exp-analysis.Rmd` is the kallisto-era pipeline: 11 bash chunks, 21 kallisto/FastQC
references, and it never reads the genome count matrix. The genome-based DESeq2 scripts
(`02_5_DESeq_*_genome`) supersede it. I moved it to `_superseded/` rather than rewrite it, since
its inputs (kallisto index, raw reads) are not committed and it is not the analysis behind the
manuscript. If you consider it still authoritative, move it back and tell me; I will rewrite it
as an HPC archive instead.

## Runnability notes

- The scripts now resolve to the committed CSVs in `02_data/` and the DEG tables in
  `03_analyses/DEG_lists/`. Run them in pipeline order (`01_5` -> `02_5_DESeq_*` ->
  `03-*_Shrinkage_filtration` -> `04-File_joining` -> `12`/`15`/`16`/`19`), because several read
  back intermediate files that earlier scripts write into `DEG_lists/`.
- **Cross-folder writes:** `16-top_DEGs` writes into `gene-annotation/03_analyses/Top_gene_summaries/`,
  a consequence of where that shared product was placed in 3c. If you would rather DE not write
  outside its own folder, I can move `Top_gene_summaries/` into `differential-expression/03_analyses/`
  instead.
- A few scripts carry Grace's `BiocManager::install(...)` lines; those are her dependency setup,
  left as-is.
- Not verified by a live run (no R/CRAN in the build environment). Knit them to confirm.

## Next in phase 4

`enrichment` and `gene-annotation` (they read DE's `DEG_lists/` and the blast GO mapping), then
the cross-cutting `gene-mechanics-correlation`. Each will reuse the `_paths.R` + `repo_root`
pattern established here.
