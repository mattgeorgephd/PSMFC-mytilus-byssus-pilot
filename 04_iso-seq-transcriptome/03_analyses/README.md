# 03_analyses

Outputs of the Iso-Seq branch (`../01_code/`); each subfolder has its own README.

| Folder | Produced by | Contents |
|---|---|---|
| `02_isoform-gene-map/` | `02_isoform_gene_map.Rmd` | each Iso-Seq isoform's genome gene (or novel locus, mitochondrial feature, or none), with a summary and the agreement with the retired CDS-based map |
| `03_salmon/` | `03_salmon_quant.Rmd` | salmon mapping summary per library and the gene counts (isoforms summed per feature with tximport); the index and per-library salmon output are git-ignored |
| `04_isoseq-de/` | `04_isoseq_de_comparison.Rmd` | the six TC contrasts fitted on the Iso-Seq gene counts and their comparison with the genome results of `06` |
| `05_augmented-annotation/` | `05_augmented_annotation.Rmd` | option B: the RefSeq annotation augmented with the isoforms (3' extension; full models), with the rules' results; the GFF and SAF files are git-ignored |
| `06_genome-recount/` | `06_genome_recount.Rmd` | option B: the reads realigned and counted on each annotation by two counters (six gene count matrices), mapping summary, checks; the HISAT2 index and per-library files are git-ignored |
| `07_augmented-de/` | `07_augmented_de_comparison.Rmd` | option B: the six TC contrasts on each recount and their comparison with the record and the RefSeq control |
| `_superseded/` | retired scripts | the abandoned kallisto attempt (`14-kallisto-ng/`) and the CDS-based isoform map (`02_isoform-gene-map_cds/`); README inside |
| `knit_html/` | the runner | HTML reports and logs (git-ignored) |
