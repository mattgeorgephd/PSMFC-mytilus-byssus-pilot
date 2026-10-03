# 01_code/_superseded

Retired scripts of the Iso-Seq branch, kept as records. Do not run them: they write to the
current output folders.

| Script | What it did | Replaced by |
|---|---|---|
| `isoseq-as-reference/` | the 2023-2024 kallisto runs of the trimmed reads against the full Iso-Seq transcriptome, with a fixed fragment length and no isoform-to-gene table, and their DE visualisations (TC and the retired LC contrasts) | `04_sequence-alignment/01_code/04_isoform_gene_map.Rmd`, then `../02_salmon_quant.Rmd` and `../03_isoseq_de_comparison.Rmd` |

The retired CDS-based isoform map (`02_isoform_gene_map_cds.Rmd`) moved with the current map
to `04_sequence-alignment/01_code/_superseded/` on 2026-10-03.
