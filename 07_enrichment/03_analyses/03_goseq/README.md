# 03_goseq

goseq (Wallenius, weighted by median transcript length), the oyster template's method; written
by `../../01_code/03_goseq.Rmd`.

| File | Contents |
|---|---|
| `goseq_enriched.csv` | terms with BH-adjusted p < 0.05 in any run, with genes in term, DEGs in term, `p_over`, `padj`, `n_mt_encoded`, the DEGs |
| `goseq_all_terms_TC_BP.csv` | every tested BP term of the TC runs; read by step 06 |
| `goseq_run_summary.csv` | per run: genes, DEGs, terms tested and enriched |
| `goseq_pwf_TC_BP.png` | the probability weighting functions of the TC runs: no consistent length trend |
| `goseq_<TC,LC,FG>_BP_dotplot.png` | the five most significant BP terms of each run |
