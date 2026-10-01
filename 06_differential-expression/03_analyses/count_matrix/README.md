# count_matrix

Written by `../../01_code/01_clean_count_matrix.Rmd`; read by every later step and by `07`
and `09`.

| File | Contents |
|---|---|
| `gene_count_matrix_clean.csv` | genes x 129 libraries: the prepDE gene matrix of `05_sequence-alignment` with columns renamed to library IDs (`T015F`), T051F and T051G removed (QC) and genes sorted by ID. Counts unchanged |
| `treatmentinfo_clean.csv` | one row per library: `sample.1`, `tissue` (F foot, G gill), `treatment`, `day`, `region` ("foot (phenol gland to tip)", "foot (without phenol gland)" for the twelve day-0 FX libraries, "gill") |
| `library_crosswalk.csv` | every library in the sample sheet with its region, whether it is in the matrix (or removed at QC, or not sequenced) and its row in the RNA isolation log (`isolation_log_sample`, `isolation_log_tissue`, `isolation_date`), matched on mussel, concentration, volume and yield |
| `mitochondrial_loci.csv` | the 143 mitochondrial rows of the matrix, which the genome fits (steps 03, 05) leave out and step 13 analyses: the 12 protein genes and 5 RNA genes of the mitochondrial genome (`source = mitogenome`) and 126 LOCs on unplaced scaffolds whose best BLAST hit (highest bitscore) is a mitochondrially encoded protein (`nuclear copy`). Columns: `gene` (matrix row), `LOC_ID` (`gene_key()` of `tools/gene_ids.R`), `source`, `mt_protein` (ND1 ... ATP6, or tRNA/rRNA), `protein_name`, `pident`, `seqid` |
