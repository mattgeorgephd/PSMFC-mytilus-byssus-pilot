# _superseded

Earlier scripts, kept as the method record; the runner does not call them.

| Script | Replaced by |
|---|---|
| `01_5-gene_count_matrix.Rmd` | `../01_clean_count_matrix.Rmd` (it read a raw matrix that was not in the repository) |
| `02_5_DESeq_Foot_TC_genome.Rmd`, `02_5_DESeq_Gill_TC_genome.Rmd` | `../02_define_contrasts.Rmd` + `../03_deseq_contrasts.Rmd` (TC) |
| `02_5_DESeq_foot_LC_genome.Rmd`, `02_5_DESeq_Gill_LC_genome.Rmd` | the same, LC contrasts |
| `02_5_DESeq_tissuevtissue_control.Rmd` | the same, FG contrasts |
| `03-LC_Shrinkage_filtration.Rmd` | `../04_shrinkage_filtration.Rmd`, which handles every family |
| `03_5-DEG_table_provenance_check.Rmd` | no longer needed: the tables are now rebuilt from the rule-defined contrasts on every run |

The interactive DESeq2 scripts could only run in one session one after another;
`03-TC_shrinkage_filtration.Rmd` (now `../04_shrinkage_filtration.Rmd`) depended on objects
they left in memory. The rewrite reproduces all 16 of their sample lists exactly and, for the
14 TC and LC contrasts they fitted, the same tested genes and DEG sets.
