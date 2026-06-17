#!/usr/bin/env bash
# Phase 3c: differential-expression, enrichment, gene-annotation, gene-mechanics-correlation.
# Completes the absorb of Grace's repo and dissolves Byssus-expression-analysis/, analyses/,
# code/, and data/. Scripts move as a faithful archive (here() rewrites are phase 4).
# Run from repo ROOT on a branch. git mv preserves history; nothing is deleted.
set -euo pipefail
[ -d .git ] || { echo "Run from repo root."; exit 1; }
B=Byssus-expression-analysis
for t in differential-expression enrichment gene-annotation gene-mechanics-correlation; do
  [ -e "$t" ] && { echo "ERROR: '$t' exists. Reset: git reset --hard && git clean -fd"; exit 1; }
done
[ -d "$B/code" ] || { echo "ERROR: $B/code missing. Current clone of main?"; exit 1; }

rproj() { cat > "$1" <<'RPROJ'
Version: 1.0

RestoreWorkspace: No
SaveWorkspace: No
AlwaysSaveHistory: No

EnableCodeIndexing: Yes
UseSpacesForTab: Yes
NumSpacesForTab: 2
Encoding: UTF-8
RPROJ
}

# ===================================================== differential-expression
D=differential-expression
mkdir -p "$D/01_code" "$D/02_data" "$D/03_analyses"
for s in 01-diff-exp-analysis 01_5-gene_count_matrix \
         02_5_DESeq_Foot_TC_genome 02_5_DESeq_Gill_LC_genome 02_5_DESeq_Gill_TC_genome \
         02_5_DESeq_foot_LC_genome 02_5_DESeq_tissuevtissue_control \
         03-LC_Shrinkage_filtration 03-TC_shrinkage_filtration 04-File_joining \
         12-DEG_venn 12-Volcano-plots 15-number_DEGS 16-top_DEGs 19-DEG_list_cleanup; do
  git mv "$B/code/$s.Rmd" "$D/01_code/$s.Rmd"
done
git mv "$B/output/gene_count_matrix_clean.csv" "$D/02_data/gene_count_matrix_clean.csv"
git mv "$B/output/gene_count_matrix_clean"     "$D/02_data/gene_count_matrix_clean"
git mv "$B/output/transcript_count_matrix.csv" "$D/02_data/transcript_count_matrix.csv"
git mv "$B/output/treatmentinfo_clean.csv"     "$D/02_data/treatmentinfo_clean.csv"
git mv "$B/data/PSMFC-mytilus-byssus-pilot-RNA-tagseq_raw.csv" "$D/02_data/PSMFC-mytilus-byssus-pilot-RNA-tagseq_raw.csv"
git mv data/psmfc_mussel_rna_summary.csv "$D/02_data/psmfc_mussel_rna_summary.csv"
git mv "$B/output/DEG_lists" "$D/03_analyses/DEG_lists"
rproj "$D/differential-expression.Rproj"

# ================================================================= enrichment
E=enrichment
mkdir -p "$E/01_code" "$E/02_data" "$E/03_analyses"
for s in 07-GOenrichment_listcreationforDAVID 08-func_enrichment_DAVID_visual 09-GOenrichment_listcreation_REVIGO; do
  git mv "$B/code/$s.Rmd" "$E/01_code/$s.Rmd"
done
git mv "$B/output/Func_annot_DAVID"     "$E/03_analyses/Func_annot_DAVID"
git mv "$B/output/Revigo_results"       "$E/03_analyses/Revigo_results"
git mv "$B/output/uniprot_BG_DAVID.txt" "$E/03_analyses/uniprot_BG_DAVID.txt"
rproj "$E/enrichment.Rproj"

# ============================================================= gene-annotation
A=gene-annotation
mkdir -p "$A/01_code/_superseded" "$A/02_data" "$A/03_analyses"
for s in Annotation 06-get_GOSlims 17-uniprot_summaries 18-ortholog-lists; do
  git mv "$B/code/$s.Rmd" "$A/01_code/$s.Rmd"
done
git mv code/06-annotation.Rmd "$A/01_code/_superseded/06-annotation.Rmd"
git mv code/06-annotation.md  "$A/01_code/_superseded/06-annotation.md"
git mv "$B/data/Foot_proteins.txt"     "$A/02_data/Foot_proteins.txt"
git mv "$B/output/Top_gene_summaries"  "$A/03_analyses/Top_gene_summaries"
rproj "$A/gene-annotation.Rproj"

# ==================================================== gene-mechanics-correlation
G=gene-mechanics-correlation
mkdir -p "$G/01_code" "$G/02_data" "$G/03_analyses"
git mv "$B/code/11-byssal_thread_by_sample.Rmd" "$G/01_code/11-byssal_thread_by_sample.Rmd"
git mv code/20-gene_mechanics_correlation.Rmd   "$G/01_code/20-gene_mechanics_correlation.Rmd"
for f in HIF_GCM HSP_GCM perox_GCM foot_byss_GCM; do
  git mv "$B/output/$f.csv" "$G/02_data/$f.csv"
done
git mv "$B/output/foot_byss_gene_plot.pdf" "$G/03_analyses/foot_byss_gene_plot.pdf"
git mv "$B/output/gill_byss_gene_plot.pdf" "$G/03_analyses/gill_byss_gene_plot.pdf"
rproj "$G/gene-mechanics-correlation.Rproj"

# ===================================================== provenance + cleanup
mkdir -p _migration/byssus-provenance
git mv "$B/README.md"      _migration/byssus-provenance/README.md
git mv "$B/code/README.md" _migration/byssus-provenance/code-README.md
git mv "$B/data/README.md" _migration/byssus-provenance/data-README.md
mkdir -p _to_delete
git mv "$B/.Rhistory"               _to_delete/byssus.Rhistory.bak
git mv "$B/project-template.Rproj"  _to_delete/byssus-project-template.Rproj.bak
git mv "$B/code/.ipynb_checkpoints" _to_delete/byssus-ipynb-checkpoints
git mv "$B/data/.gitignore"         _to_delete/byssus-data.gitignore.bak
git mv data/README.md               _to_delete/root-data-README.md.bak
git mv data/.gitignore              _to_delete/root-data.gitignore.bak
git mv analyses/.gitignore          _to_delete/analyses.gitignore.bak

rmdir "$B/code" "$B/data" "$B/output" "$B" code data analyses 2>/dev/null || true
echo "Phase 3c staged."
