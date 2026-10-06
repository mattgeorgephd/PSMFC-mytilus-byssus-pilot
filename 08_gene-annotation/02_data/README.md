# 02_data

| Item | Description | Read by |
|------|-------------|---------|
| `goslim_generic.obo` | The generic GO slim terms of GO release 2026-01-23 (141 terms; 72 biological process), extracted from the `subset: goslim_generic` tags of `geneontology/go-ontology` `src/ontology/go-edit.obo` at the last commit before that release (see the file's header), by `../01_code/_derive_goslim.R`. Pinned so the slim does not change between runs; it must be of the release in `GO.db` (step 01 checks). To use another release, set the release and commit in `_derive_goslim.R` and rerun it, or use the official `goslim_generic.obo` of that release | `01_go_slims.Rmd` |
| `_superseded/` | the slim of GO release 2023-07-27, used with `GO.db` 3.18.0 until 2026-10-03 (README inside) | no current script |
| `Foot_proteins.txt` | FASTA of byssal foot-protein coding sequences from GenBank (nucleotide), a reference for annotation | no current script |

The annotated DEG tables and the top-50 lists are read cross-folder from
`../../05_differential-expression/03_analyses/` (see `../README.md`).
