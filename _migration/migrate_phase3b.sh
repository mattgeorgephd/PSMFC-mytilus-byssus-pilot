#!/usr/bin/env bash
# Phase 3b: sequence-alignment, blast, iso-seq-transcriptome (absorb of Grace's repo, part 1).
# HPC scripts are moved unchanged (faithful archive); only 05-IsoSeq is rewritten separately.
# Run from repo ROOT on a branch. git mv preserves history; nothing is deleted.
set -euo pipefail
[ -d .git ] || { echo "Run from repo root."; exit 1; }
for t in sequence-alignment blast iso-seq-transcriptome; do
  [ -e "$t" ] && { echo "ERROR: '$t' exists. Reset: git reset --hard && git clean -fd"; exit 1; }
done
B=Byssus-expression-analysis
for s in "$B/code/13-Hisat.Rmd" "$B/code/07-HiSat_GL.Rmd" code/07-kallisto.Rmd code/05-IsoSeq-transcriptome-check.Rmd \
         analyses/07-kallisto analyses/blast "$B/output/13-Hisat" "$B/output/fastqc" \
         "$B/output/LOC_GO_list.txt" "$B/output/g.spid.txt" \
         "$B/code/Initial analysis with isoseq as reference (not final)"; do
  [ -e "$s" ] || { echo "ERROR: source missing: $s"; exit 1; }
done

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

# ============================================================ sequence-alignment
mkdir -p sequence-alignment/01_code/_superseded sequence-alignment/02_data/sample-submission sequence-alignment/03_analyses/_superseded
git mv "$B/code/13-Hisat.Rmd"        sequence-alignment/01_code/13-Hisat.Rmd
git mv "$B/code/07-HiSat_GL.Rmd"     sequence-alignment/01_code/_superseded/07-HiSat_GL.Rmd
git mv code/07-kallisto.Rmd          sequence-alignment/01_code/_superseded/07-kallisto.Rmd
git mv code/07-kallisto.md           sequence-alignment/01_code/_superseded/07-kallisto.md
git mv tagseq/* sequence-alignment/02_data/sample-submission/
rmdir tagseq 2>/dev/null || true
git mv "$B/output/13-Hisat" sequence-alignment/03_analyses/hisat
git mv "$B/output/fastqc"   sequence-alignment/03_analyses/fastqc
git mv analyses/07-kallisto sequence-alignment/03_analyses/_superseded/kallisto
rproj sequence-alignment/sequence-alignment.Rproj

# ======================================================================== blast
mkdir -p blast/01_code blast/02_data blast/03_analyses
git mv "$B/code/Mtros-genome-blast.Rmd"     blast/01_code/Mtros-genome-blast.Rmd
git mv "$B/code/Isoseq_vs_genome_blast.Rmd" blast/01_code/Isoseq_vs_genome_blast.Rmd
git mv code/uniprot-retrieval.py            blast/01_code/uniprot-retrieval.py
git mv analyses/blast blast/03_analyses/transcriptome-uniprot
mkdir -p blast/03_analyses/genome-foot
git mv "$B/output/LOC_GO_list.txt" blast/03_analyses/genome-foot/LOC_GO_list.txt
git mv "$B/output/g.spid.txt"      blast/03_analyses/genome-foot/g.spid.txt
rproj blast/blast.Rproj

# ======================================================== iso-seq-transcriptome
mkdir -p iso-seq-transcriptome/01_code/_superseded iso-seq-transcriptome/02_data iso-seq-transcriptome/03_analyses/_superseded
git mv code/05-IsoSeq-transcriptome-check.Rmd       iso-seq-transcriptome/01_code/05-IsoSeq-transcriptome-check.Rmd
git mv code/05-IsoSeq-transcriptome-check.md        iso-seq-transcriptome/01_code/05-IsoSeq-transcriptome-check.md
git mv code/05-IsoSeq-transcriptome-check_files     iso-seq-transcriptome/01_code/05-IsoSeq-transcriptome-check_files
git mv "$B/code/Initial analysis with isoseq as reference (not final)" iso-seq-transcriptome/01_code/_superseded/isoseq-as-reference
git mv "$B/output/14-kallisto-ng" iso-seq-transcriptome/03_analyses/_superseded/14-kallisto-ng
rproj iso-seq-transcriptome/iso-seq-transcriptome.Rproj

echo "Phase 3b staged."
