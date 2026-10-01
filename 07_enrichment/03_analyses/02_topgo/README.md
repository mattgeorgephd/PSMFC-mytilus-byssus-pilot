# 02_topgo

topGO `weight01` with Fisher's exact test, the enrichment of record; written by
`../../01_code/02_topgo.Rmd`.

| File | Contents |
|---|---|
| `topgo_enriched.csv` | every term with `weight01` p < 0.01 in any run (7 contrasts x up/down x BP/MF/CC): annotated and significant genes, expected, `p_weight01`, `p_classic`, and the significant genes in the term (`genes`, matched to the term by GO ID) |
| `topgo_all_terms_TC_<BP,MF,CC>.csv` | every tested term of the 12 TC runs, p-values only, one file per ontology; read by steps 05 and 06 |
| `topgo_run_summary.csv` | per run: annotated genes, DEGs, terms tested, terms enriched |
| `topgo_<TC,FG>_<BP,MF,CC>_dotplot.png` | the five most significant terms of each run, one figure per family and ontology |
