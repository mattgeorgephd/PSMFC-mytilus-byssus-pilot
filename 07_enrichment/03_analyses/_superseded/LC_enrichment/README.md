# LC_enrichment

GO enrichment of the retired LC contrasts (each stressor, and the day-3 treatment control,
against the day-0 lab control), as written by `../../../01_code/` at commit 5ba5614, before
the LC contrasts were retired (2026-10-01). Nothing reads these files.

The day-0 animals are not a control: their feet were dissected differently from the day-3
feet (see `06_differential-expression/01_code/02_define_contrasts.Rmd`). These results also
predate the separation of the mitochondrial loci, so they include them.

| File | Contents |
|---|---|
| `topgo_LC_BP_dotplot.png`, `goseq_LC_BP_dotplot.png`, `clusterprofiler_LC_BP_dotplot.png`, `rrvgo_LC_BP_parents.png` | the LC dot plots (biological process) of each method |
| `LC_<table>.csv` | the LC rows of the result tables of that commit: `topgo_enriched`, `topgo_run_summary`, `goseq_enriched`, `goseq_run_summary`, `clusterprofiler_enriched`, `clusterprofiler_run_summary`, `rrvgo_parents`, `rrvgo_reduced_terms`, `gene_sets_summary` |

What replaced them: nothing; the TC contrasts are the only stressor contrasts, and
`FG_TC` the only foot-vs-gill contrast.
