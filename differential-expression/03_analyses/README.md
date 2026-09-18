# 03_analyses

| Subfolder / file | Produced by | Contents |
|-----------|-------------|----------|
| `DEG_lists/Foot/`, `DEG_lists/Gill/` | `02_5_DESeq_*`, `03-*_Shrinkage_filtration` | Per-contrast count and treatment tables, full apeglm results (`<X>_TC_apeglm.csv`; `GOA_TC.csv`), DEG lists (`<X>_TC_siggene*.csv`, `padj < 0.05`), MA plots |
| `DEG_lists/DEG_provenance_check.csv` | `03_5-DEG_table_provenance_check` | Per contrast (six TC, eight LC; `family` column): does each result table on disk reproduce from the committed inputs |
| `DEG_lists/secretion_state_by_arm.csv`, `secretion_state_summary.csv`, `secretion_state_day0_vs_day3.png` | `01_7-secretion_state` | Byssal secretion state (`on`/`off`) counts by arm with Fisher tests; day-0 vs day-3 medians of the module score and each module gene; figure. The per-sample states are in `02_data/secretion_state.csv` |
| `DEG_lists/sensitivity_fourlevel/` | `02_6-DESeq_fourlevel_sensitivity` | One four-level DESeq2 model per tissue: `<X>_TC4_apeglm.csv`, `<X>_TC4_siggene.csv`, and `fourlevel_vs_pairwise_DEG_counts.csv` against the pairwise lists of record |
| `DEG_lists/sensitivity_secretion/` | `02_7-DESeq_foot_TC_secretion_sensitivity` (gated, off by default) | Foot stressor contrasts refitted with the secretion state as covariate; `secretion_sensitivity_DEG_counts.csv` |
| `DEG_lists/GOterms_genome/` | `04-File_joining` | DEG lists joined to the BLAST/UniProt/GO annotation: `_sigs_merged` (gene x hit), `_sigs_ID`, `_sigs_unID`; `clean_zenodo_files/` from `19-DEG_list_cleanup` |
| `DEG_lists/DEG_join_summary.csv` | `04-File_joining` | Per contrast: DEG genes, merged rows, annotated / unannotated genes, mitochondrial DEGs. **Report `n_DEG_genes`, not merged rows.** |
| `DEG_lists/DAVID_lists/`, `DEG_lists/REVIGO_lists/` | `enrichment/01_code/07-*`, `09-*` | Accession and GO-ID lists submitted to DAVID / REVIGO |
| `DEG_lists/goslims_genome/` | `gene-annotation/01_code/06-get_GOSlims` | GO-slim mappings per contrast |

`DEG_lists/` is an intermediate: the DESeq scripts write it and the joining, venn, volcano,
counts, and top-genes scripts read it back. Enrichment (in `../../enrichment/`) and the
gene-mechanics pipeline (`../../gene-mechanics-correlation/`, which builds its DEG union from
the `_TC_siggene*` files and its annotation map from `GOterms_genome/*_sigs_ID.csv`) also
read it.
