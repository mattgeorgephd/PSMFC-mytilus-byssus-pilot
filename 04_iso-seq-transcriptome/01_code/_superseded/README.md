# 01_code/_superseded

Retired scripts of the Iso-Seq branch, kept as records. Do not run them: they write to the
current output folders.

| Script | What it did | Replaced by |
|---|---|---|
| `isoseq-as-reference/` | the 2023-2024 kallisto runs of the trimmed reads against the full Iso-Seq transcriptome, with a fixed fragment length and no isoform-to-gene table, and their DE visualisations (TC and the retired LC contrasts) | `../02_isoform_gene_map.Rmd` to `../04_isoseq_de_comparison.Rmd` |
| `02_isoform_gene_map_cds.Rmd` | the first isoform-to-gene map (2026-10-01), built while the genome was not reachable: each isoform aligned to the annotation's CDS (minimap2 `-x asm20`), identity at least 0.90 over at least 100 bases. It reaches coding genes only and misses isoforms that end before the CDS. Output: `../../03_analyses/_superseded/02_isoform-gene-map_cds/` | `../02_isoform_gene_map.Rmd` (spliced alignment to the genome and exon overlap with every annotated gene), which compares its map with this one |
