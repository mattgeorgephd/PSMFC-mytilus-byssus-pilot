# 02_data

| Item | Description | Read by |
|------|-------------|---------|
| `goslim_generic.obo` | The generic GO slim terms of GO release 2023-07-27 (141 terms; 72 biological process), extracted from the `subset: goslim_generic` tags of `geneontology/go-ontology` `src/ontology/go-edit.obo` at that release (see the file's header). Pinned so the slim does not change between runs; to use another release, replace it with the official `goslim_generic.obo` from current.geneontology.org | `01_go_slims.Rmd` |
| `Foot_proteins.txt` | FASTA of byssal foot-protein coding sequences from GenBank (nucleotide), a reference for annotation | no current script |

The annotated DEG tables and the top-50 lists are read cross-folder from
`../../05_differential-expression/03_analyses/` (see `../README.md`).
