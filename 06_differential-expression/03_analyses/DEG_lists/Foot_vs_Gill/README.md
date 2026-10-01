# Foot_vs_Gill

Gill against foot within the controls, written by `../../../01_code/04_shrinkage_filtration.Rmd`
from the fits of `03_deseq_contrasts.Rmd`. Positive log2 fold changes are higher in gill.

| Contrast | Samples |
|---|---|
| `FG_TC` | the 12 day-3 treatment controls, foot and gill |

The day-0 version (`FG_LC`) is retired with the other day-0 contrasts; its tables are in
`../../_superseded/LC_contrasts/Foot_vs_Gill/`.

Per contrast: `<code>_apeglm.csv` (every gene kept by the count filter), `<code>_siggene.csv`
(padj < 0.05), `<code>_filter_counts.csv` and `<code>_MA_plots.pdf`. Volcano plots:
`../../figures/volcano_FG.png`. GO enrichment of these contrasts: `07_enrichment`.
