# prepDE

The count matrices of the original HPC alignment (HISAT2 + StringTie + prepDE): the count
matrices of record until 2026-10-02, when `../featurecounts/` replaced them. They are still
written: `04_iso-seq-transcriptome` step 07 refits the gene matrix with `06`'s rules as the
previous record it compares with, `04_iso-seq-transcriptome` step 06 checks its recount against
it, `06_differential-expression` step 13 compares the mitochondrial counts with it, and step 04
here compares gene names with it.

| File | Produced by | Contents |
|---|---|---|
| `transcript_count_matrix.csv` | `prepDE.py` on the HPC, from the StringTie output of `../../01_code/01_hisat_stringtie.Rmd` | transcripts x 131 libraries; read-length 75. Moved here from `06_differential-expression/02_data/` |
| `gene_count_matrix.csv` | `../../01_code/02_prepDE.Rmd` | genes x 131 libraries. With the StringTie GTFs present (`../hisat/sample_list.txt`) it counts them with `_prepde.R`, an R port of `prepDE.py3`; without them it sums the committed transcript matrix to genes through `../hisat/t_data.ctab` and `../../02_data/strg_gene_ids.csv` |

Counts are `ceiling(coverage x transcript length / read length)`, as `prepDE.py` computes them.
