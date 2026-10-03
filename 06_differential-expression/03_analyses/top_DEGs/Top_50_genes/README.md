# Top_50_genes

`<T><X>_topgenes.csv` for the six TC contrasts (FOA, FOW, FDO, GOA, GOW, GDO): the 25
most up-regulated and the 25 most down-regulated annotated DEGs (by apeglm log2 fold change; fewer
on a side that has fewer), with their BLAST / UniProt annotation. Written by `../../../01_code/07_top_degs.Rmd` (moved here from
`08_gene-annotation/03_analyses/Top_gene_summaries/`, so this folder writes only to itself).

One row per gene, carrying its best UniProt hit (highest bitscore; ties in table order), the
rule `07_enrichment` and `09` use. Until 2026-10-01 a gene with several hits (one per
transcript) took whichever came first in the table; that changed the shown hit for a handful
of genes, among them LOC134695253 in GDO, labelled `CO6A3` (collagen VI) rather than its best
hit `LRP6`.

`<code>_top50.png` draws each table; a bar's label is `gene_symbol`: the first gene name
UniProt lists that is not a locus tag, upper-cased, or a shortened protein name where every
name is a locus tag or there is none (`MCOR_42823`, a *M. coruscus* locus tag, becomes `ACDC`,
the byssal adhesive it is). Where two genes in one figure share a label, the LOC ID is added.
`geneID` is the first half of the UniProt entry name, a mnemonic that is often not the gene
symbol (`CO6A3` for COL6A3); it is kept in the table but no longer used as a label.
