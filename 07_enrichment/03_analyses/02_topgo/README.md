# 02_topgo

topGO `weight01` with Fisher's exact test, the enrichment of record; written by
`../../01_code/02_topgo.Rmd`.

| File | Contents |
|---|---|
| `topgo_enriched.csv` | every term with `weight01` p < 0.01 in any run (16 contrasts x up/down x BP/MF/CC): annotated and significant genes, expected, `p_weight01`, `p_classic`, `n_mt_encoded`, the significant genes |
| `topgo_all_terms_TC_BP.csv` | every tested BP term of the 12 TC runs, p-values only; read by steps 05 and 06 |
| `topgo_run_summary.csv` | per run: annotated genes, DEGs, terms tested, terms enriched |
| `topgo_TC_BP_dotplot.png`, `topgo_LC_BP_dotplot.png`, `topgo_FG_BP_dotplot.png` | the five most significant BP terms of each run, by family. Triangles: half or more of the term's DEGs are mitochondrially encoded protein LOCs |
