# LC_contrasts

The lab-control (LC) contrasts, retired on 2026-10-01: each stressor against the day-0 lab
controls, the day-3 treatment control against the day-0 lab controls, and gill vs foot within
the day-0 lab controls (`FG_LC`). Written by `../../../01_code/04_shrinkage_filtration.Rmd` and
`11_deg_figures_all_contrasts.Rmd` at commit 5ba5614, the last run that defined them.

**Why they were retired.** The day-0 lab-control animals are not true lab controls: their feet
were dissected differently from the day-3 feet (two pieces, `F_PG` and `F` in the RNA isolation
log, against one piece, `F_PG`, at day 3), as the experimenter confirmed. A day-0 vs day-3
difference therefore mixes the part of the foot sampled into the effect, and the day-0
animals are used as a control nowhere in the analysis. The contrasts of record were always
the day-3 treatment-control (TC) contrasts, which do not change.

| Item | Contents |
|---|---|
| `LC_contrasts.csv`, `LC_contrast_samples.csv`, `LC_DEG_counts.csv` | the nine retired contrasts' definitions, samples and DEG counts (the LC rows of `../../DEG_lists/contrasts.csv`, `contrast_samples.csv` and `DEG_counts.csv` from that run) |
| `Foot/`, `Gill/`, `Foot_vs_Gill/` | per contrast: `<code>_apeglm.csv`, `<code>_siggene.csv`, `<code>_filter_counts.csv`, `<code>_MA_plots.pdf` |
| `DEG_overlap_TC_LC.csv`, `figures/TC_vs_LC_overlap.png` | DEGs found against the day-3 control, the day-0 control or both |
| `figures/volcano_LC.png` | volcano plots of the LC contrasts |

These tables still include the mitochondrial loci, which the current genome analysis leaves out.
