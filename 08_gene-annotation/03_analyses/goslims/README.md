# goslims

Written by `../../01_code/01_go_slims.Rmd`.

| File | Contents |
|---|---|
| `<code>_sigs_ID.tab` | per TC contrast (FOA, FOW, FDO, GOA, GOW, GDO), one row per BP slim term: genes in it (all, up, down), the GO IDs that mapped, the genes' UniProt entry-name mnemonics (`Genes`, e.g. `CO6A3` for COL6A3) and LOC IDs |
| `<code>_sigs_ID.BP_per_gene.tab` | one row per DEG with its slim terms |
| `goslim_summary_TC.csv` | the six contrasts in one long table, the data behind the heatmap |
| `goslim_TC_heatmap.png` | DEGs per slim term and contrast, up (red) and down (blue); terms with fewer than five DEGs over all contrasts are left out of the figure |
| `goslim_provenance.txt` | the slim and `GO.db` versions behind the run |

A gene counts under every slim term any of its GO IDs falls under, so the counts do not add
up to the number of DEGs.
