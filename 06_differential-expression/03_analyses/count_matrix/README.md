# count_matrix

Written by `../../01_code/01_clean_count_matrix.Rmd`; read by every later step and by `07`
and `09`.

| File | Contents |
|---|---|
| `gene_count_matrix_clean.csv` | genes x 129 libraries: the count matrix of record of `05_sequence-alignment` (featureCounts on the Iso-Seq-extended annotation, `03_analyses/featurecounts/`), T051F and T051G removed (QC) and genes sorted by ID. Counts unchanged. (Until 2026-10-02, the prepDE gene matrix of the HPC alignment.) |
| `treatmentinfo_clean.csv` | one row per library: `sample.1`, `tissue` (F foot, G gill), `treatment`, `day`, `region` ("foot (phenol gland to tip)", "foot (without phenol gland)" for the twelve day-0 FX libraries, "gill") |
| `library_crosswalk.csv` | every library in the sample sheet with its region, whether it is in the matrix (or removed at QC, or not sequenced) and its row in the RNA isolation log (`isolation_log_sample`, `isolation_log_tissue`, `isolation_date`), matched on mussel, concentration, volume and yield |
| `mitochondrial_loci.csv` | the 331 mitochondrial rows of the matrix, which the genome fits (steps 03, 05) leave out: the 12 protein genes and the 26 tRNAs and rRNAs of the mitochondrial genome (`source = mitogenome`; the previous matrix named 21 of the RNAs by StringTie number, so only 5 were listed until 2026-10-02) and 293 copies on unplaced scaffolds (`nuclear copy`): 126 protein-coding LOCs whose best BLAST hit (highest bitscore) is a mitochondrially encoded protein and which NCBI names after it, and 167 pseudogenes NCBI names after one (`05_sequence-alignment/02_data/annotation_mt_like_loci.csv`). Columns: `gene` (matrix row), `LOC_ID` (`gene_key()` of `tools/gene_ids.R`), `source`, `mt_protein` (ND1 ... ATP6, or tRNA/rRNA), `protein_name` (the BLAST hit, or NCBI's description), `pident`, `seqid`, `evidence` (`mitogenome`, `BLAST and annotation`, `annotation`) |
