# byssus_individuals

Written by `../../01_code/07_byssal_individuals.Rmd`, once per tissue (`_F` foot, `_G` gill),
from script 04's byssal genes and variance-stabilised counts (`../byssus_genes/`), script 01's
paired-sample manifest (`../gene_mechanics/`), `05` step 01's count matrix, sample table,
mitochondrial list and crosswalk to the RNA isolation log (`05/02_data/psmfc_mussel_rna_summary.csv`,
for the RNA yield), `05` step 14's PCA scores, the protein names of the annotation of record
(`03_blast` step 05) and the thread counts in `01_mussel-measurements/mussel-thread-production.xlsx`.
The question: within each arm, are the byssal genes high or low together in the same animals,
which mussels are high for most mussel foot proteins (mfp-1 to -6), and does that state go with
library size, thread production or thread strength? Every statistic is within arms, so no
treatment effect enters; the day-0 lab control is shown only as a reference.

| File | Contents |
|---|---|
| `byssal_animals_<T>.csv` | one row per day-3 animal: `sample`, `mussel`, `treatment`, `lib_size`, `genes_detected`, `mfp_vst` (mean VST of the expressed mfp genes), `mfp_read_pct` (their share of the library's reads, %), `mfp_n_high` (how many of them put the animal in its arm's top third), `mfp_high` (more than half of them), `mfp_on_posterior` and `mfp_on` (the posterior of the upper component of the two-component mixture on `mfp_vst`, and whether it is above 0.5), `score_<group>` (mean within-arm z-score of each gene group with at least 3 expressed genes), `threads` (made in the 3 days), day-3 and baseline `mean_force` and `pad_area`, `pca_pc1` (PC1 of `05` step 14's "every arm" set), `rna_yield_ng` (RNA yield of the animal's isolation) |
| `byssal_coordination_<T>.csv` | one row per gene group: `genes`, `mean_r` (mean pairwise correlation of its genes across the animals of the same arm), the within-arm permutation null's median and 95th percentile, `p` |
| `byssal_group_correlations_<T>.csv` | Spearman correlations between the group scores (`group_a`, `group_b`, `rho`) |
| `byssal_mfp_tests_<T>.csv` | the number of mfp-high animals, and of animals high for every expressed mfp gene, against the permutation null; the mixture against a single normal (delta BIC), the mfp-on animals by arm and the Fisher exact test of arm against mfp-on; the mfp mean against library size, genes detected and PC1; the mfp score against threads made (`log(1 + threads) ~ score + arm`), day-3 mean peak force and plaque area (`log(day 3) ~ score + arm + log(baseline)`) and RNA yield (`log(yield) ~ score + arm`, a technical check), as partial r, slope, SE and p. In gill a single row says why it was not tested (fewer than 3 expressed mfp genes) |
| `byssal_mfp_genome_<T>.csv` | foot: every nuclear gene with at least 10 counts in a third of the animals (the expressed mfp genes left out) against the mfp score within arms (`VST ~ score + arm`): `partial_r`, `slope` (VST per within-arm SD of the score), `p`, `q` (BH), `protein_name` (best hit), `byssal_group` (script 04's group, if any) |
| `byssal_mfp_genome_<T>.png` | foot: A, every gene tested (partial r against -log10 p; script 04's byssal genes ringed); B, the genes at q < 0.05 by name, filled by kind (byssal gene, other annotated gene, no BLAST hit) |
| `byssal_heatmap_structural_by_id_<T>.png`, `byssal_heatmap_all_by_id_<T>.png` | the byssal structural genes (mfp-1 to -6, thread core, other byssal structural), and every gene of script 04, for every day-3 animal; columns by mussel ID within each arm, the same order in every panel and in both figures, so that one animal can be followed across gene groups |
| `byssal_heatmap_structural_by_mfp_<T>.png`, `byssal_heatmap_all_by_mfp_<T>.png` | the same, columns by mean mfp expression within each arm, highest first (in gill, by mean byssal structural expression) |
| `byssal_individuals_<T>.png` | foot: A, mean mfp expression by arm with the mixture's cut (dashed) and the mfp-on animals labelled, the lab control beside; B, coordination of each group against its permutations; C, the group scores' correlations; D, mfp score against threads made. Gill: B and C |
| `RUN_provenance_<T>.txt` | settings, coordination, the mfp-high and mfp-on animals by arm, the mfp share of reads, the mfp tests, the genome-wide counts, `checks` (the recomputed VST is script 04's, the permutation nulls are centred near 0, the mfp genes are script 04's, the mixture converged, the partial correlations' p equal `lm`'s), commit, input MD5s, R and package versions |

In the heatmaps a `+` after a row label marks an mfp gene below script 04's expression filter
(drawn, not tested) and `*` a hit below 40% identity; the top row counts the expressed mfp genes
for which the animal is in its arm's top third, and outlines the mfp-on animals.
