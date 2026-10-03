# 03_analyses

Outputs of the Iso-Seq branch (`../01_code/`); each subfolder has its own README.

| Folder | Produced by | Contents |
|---|---|---|
| `02_salmon/` | `02_salmon_quant.Rmd` | salmon mapping summary per library and the gene counts (isoforms summed per feature with tximport); the index and per-library salmon output are git-ignored |
| `03_isoseq-de/` | `03_isoseq_de_comparison.Rmd` | the six TC contrasts fitted on the Iso-Seq gene counts and their comparison with the genome results of `05` |
| `04_augmented-de/` | `04_augmented_de_comparison.Rmd` | option B: the six TC contrasts on each recount of `04_sequence-alignment` step 06 and their comparison with the previous record and the RefSeq control |
| `_superseded/` | retired scripts | the abandoned kallisto attempt (`14-kallisto-ng/`); README inside |
| `knit_html/` | the runner | HTML reports and logs (git-ignored) |

The isoform-to-gene map, the augmented annotations and the genome recounts (until 2026-10-03
`02_isoform-gene-map/`, `05_augmented-annotation/` and `06_genome-recount/` here) are in
`04_sequence-alignment/03_analyses/` (`isoform-gene-map/`, `augmented-annotation/`,
`genome-recount/`), with the retired CDS-based map (`_superseded/02_isoform-gene-map_cds/`).
