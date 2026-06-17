# Phase 4 rewrite: `enrichment` and `gene-annotation`

Two apply scripts, same model as `apply_phase4_de.sh`: Perl-based (ships with Git for Windows, no
Python), CRLF-safe, with an idempotency guard that aborts if `_paths.R` already exists. Each adds a
shared `01_code/_paths.R`, inserts a `source()` chunk after each script's YAML, and repoints paths.

Tested in a clean checkout and against Windows CRLF line endings: 0 residual bad paths in the
rewritten scripts, parentheses balanced in every file.

---

## How to run

Run **both** on the same branch, review the combined diff, push once:

```
git switch -c reorg/phase4-enrich-annot      # skip if you already made the branch
bash apply_phase4_enrichment.sh
bash apply_phase4_gene-annotation.sh
git status
git diff
git push -u origin reorg/phase4-enrich-annot
```

If a script aborts with "`_paths.R already exists`" (e.g. a half-run), reset first:
`git checkout . && git clean -fd`, then re-run.

---

## `enrichment` (3 scripts, all pure R)

`_paths.R` defines:
- `deg`  -> `differential-expression/03_analyses/DEG_lists` (cross-folder: sig lists in `GOterms_genome/`, plus the DAVID/REVIGO input lists this analysis writes under `DAVID_lists/` and `REVIGO_lists/`)
- `dat`  -> `differential-expression/02_data` (`gene_count_matrix_clean.csv`)
- `blast_go` -> `blast/03_analyses/genome-foot` (`g.spid.txt`, `LOC_GO_list.txt`)
- `func_david` -> `enrichment/03_analyses/Func_annot_DAVID` (own output, `dir.create`d)

Pipeline order (DAVID and REVIGO themselves are external web tools, run by hand between scripts):
1. `07-GOenrichment_listcreationforDAVID` -> writes the gene-ID lists for DAVID
2. (run DAVID on the web, save chart `.xlsx` into `Func_annot_DAVID/`)
3. `08-func_enrichment_DAVID_visual` -> reads those `.xlsx`, makes figures
4. `09-GOenrichment_listcreation_REVIGO` -> writes REVIGO input lists (independent of 07/08)

## `gene-annotation` (06/17/18 rewritten; `Annotation.Rmd` left unchanged)

`_paths.R` defines:
- `deg` -> `differential-expression/03_analyses/DEG_lists` (reads `GOterms_genome/`; writes `goslims_genome/`)
- `topgenes` -> `gene-annotation/03_analyses/Top_gene_summaries` (own output; `Top_50_genes/` subdir `dir.create`d)
- `Sys.setenv(GOTERMS=...)` for the disabled bash chunk in 06 (see decision 2)

Pipeline order (note the cross-folder dependency on DE's script 16):
1. DE `16-top_DEGs` must have run first; it writes `Top_gene_summaries/Top_50_genes/*_topgenes.csv`
2. `17-uniprot_summaries` -> reads `Top_50_genes/`, writes `*_topgene_summs.csv`
3. `18-ortholog-lists` -> reads `*_topgene_summs.csv`, writes `*_ortho.csv`
4. `06-get_GOSlims` -> reads `GOterms_genome/`, writes `goslims_genome/` (independent of 17/18)

---

## Decisions baked in (reversible; flagging so you can override on review)

1. **`Annotation.Rmd` left byte-for-byte unchanged.** It is an HPC blast script (7 bash chunks, 6
   `makeblastdb`/`blastx`, reads an uncommitted `Mtros-genome-uniprot_blastx.tab`), so it is not
   locally runnable and is preserved as a faithful archive, exactly as the blast scripts were in
   phase 3b. It still contains its original `/home/shared/...` paths by design.
   - Open question I did **not** decide: it overlaps `blast/Mtros-genome-blast.Rmd` (both blastx the
     genome against UniProt). If one supersedes the other, say which and I will move the loser to a
     `_superseded/`.

2. **`06-get_GOSlims` bash chunk set to `eval=FALSE` and path-repointed via `$GOTERMS`.** That chunk
   runs `awk -F"\t" '{print $2,"\t",$1}' ... > FDO_sigs_ID.csv` (in-place). The committed
   `FDO_sigs_ID.csv` is now a 39-column **comma-separated** table, so tab-field-splitting would blank
   `$2` and overwrite the file with garbage. The R parts of `06` read these CSVs directly with
   `read.csv()` and do **not** consume the bash chunk's output, so disabling it changes no results.
   - Options if you want to act on it: (a) keep as-is (disabled, preserved); (b) delete the chunk
     entirely since it is stale; (c) if it still serves a purpose, tell me the intended input format
     and I will rewrite it correctly and re-enable.

3. **Still open from phase 3c: `Top_gene_summaries` placement.** It currently lives in
   `gene-annotation/03_analyses/`, but DE's `16-top_DEGs` writes into it (a cross-folder write), and
   `17`/`18` here read/write it. If you would rather avoid the cross-folder write, I can relocate
   `Top_gene_summaries` into `differential-expression/03_analyses/` and repoint both folders. Your call.
