# goslims_genome

`<code>_sigs_ID.tab` and `<code>_sigs_ID.BP_per_gene.tab` from the earlier GO slim script.
They are incomplete: each gene contributed only its first-listed GO ID (the others kept a
leading space that `GSEABase::GOCollection()` drops), and genes were looked up through only
the first GO ID of each slim term, so the tables hold 36-65% of the gene-to-slim links implied
by their own `GO.IDs` column. The corrected tables are in `../../goslims/`.
