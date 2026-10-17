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
| `byssal_heatmap_original_v1_F.png`, `byssal_heatmap_original_v2_F.png` | foot: the first two heatmaps of these genes, drawn with the talk's code (October 2026), animals ordered within each arm by mean expression, no mussel IDs. v1: script 04's byssal structural genes that pass its filter (19 then, 22 now), genes by mean expression; v2: the talk's list (every mfp-1 to -6 gene with reads, `+` below the filter, `*` found by `03_blast` step 05; thread core and other byssal genes that pass the filter), animals by mean mfp expression. Each with the mean by arm (Kruskal-Wallis) and against threads made (Spearman), the talk's pooled tests. The port reproduces the talk's figures pixel for pixel (v1 on the byssal table of 0d8520c, which the talk used) |
| `byssal_individuals_<T>.png` | foot: A, mean mfp expression by arm with the mixture's cut (dashed) and the mfp-on animals labelled, the lab control beside; B, coordination of each group against its permutations; C, the group scores' correlations; D, mfp score against threads made. Gill: B and C |
| `RUN_provenance_<T>.txt` | settings, coordination, the mfp-high and mfp-on animals by arm, the mfp share of reads, the mfp tests, the genome-wide counts, `checks` (the recomputed VST is script 04's, the permutation nulls are centred near 0, the mfp genes are script 04's, the mixture converged, the partial correlations' p equal `lm`'s), commit, input MD5s, R and package versions |

In the heatmaps a `+` after a row label marks an mfp gene below script 04's expression filter
(drawn, not tested) and `*` a hit below 40% identity; the top row counts the expressed mfp genes
for which the animal is in its arm's top third, and outlines the mfp-on animals.

## Results (run of 610f1dc, 2026-10-17)

Foot: 46 day-3 animals (control 12, OA 12, OW 12, DO 10), 154 byssal genes of script 04, 140
expressed; 2,000 within-arm permutations (p 0.0005 is the smallest possible).

- **In every arm a few mussels express the mfp genes and the rest barely do.** Mean expression of
  the six expressed mfp genes (mfp-1, mfp-2, two mfp-3, mfp-4, mefp-5) is either near the floor of
  the transform or far above it: two normals fit far better than one (delta BIC 20.2), at VST 2.96
  and 7.92, about 30-fold apart. 11 of 46 animals express them: control T126 and T133, OA T015 and
  T118, OW T029, T031, T037 and T040, DO T049, T057 and T058 (VST 6.0 to 9.6; the others 5.6 or
  less). The share does not differ between arms (2, 2, 4 and 3; Fisher p 0.73). Four animals lie
  near the boundary (posterior for "on" 0.19 to 0.78: T037, T112, T127, T128). In the mfp-on
  animals the six genes take a median 0.57% of the library's reads (up to 0.94%), in the others
  0.00% (up to 0.13%).
- **The same animals are high for every mfp gene.** The six genes' mean within-arm correlation
  is 0.81 (permutations: 95th percentile 0.07). 10 animals are in their arm's top third for all
  six genes, where chance gives 0.1; 16 for at least four (5.4 by chance; the top third holds
  about four animals per arm for each gene, so this count is set largely by the rule).
- **The state reaches the other plaque and cuticle genes, less the thread core.** Mean within-arm
  correlation: other byssal structural genes 0.56 (foot proteins 9 to 17, YGH-rich proteins),
  thread core 0.26 (precollagens D, P and NG, thread matrix proteins), the other groups 0.02 to
  0.12. The animals' mfp score correlates with their score for the other byssal structural genes
  (Spearman 0.87), byssal accessory (0.68), collagen (0.56), thread core (0.50) and cross-linking
  and curing genes (0.46).
- **Genome-wide it is a byssal program, not a rise in all mRNA.** Of 7,791 nuclear genes tested
  against the mfp score within arms, 21 rise at q < 0.05 and none falls: 15 byssal genes (foot
  proteins 10, 11, 12 and 15, YGH-rich proteins 1 and 3, collagen-like proteins 5 to 7, byssal
  tyrosinase-like protein 1, TSP_1 domain-containing protein 1, precollagen D, nongradient
  byssal), 5 genes with no BLAST hit (partial r 0.56 to 0.95, candidates for byssal genes the
  annotation misses) and the ER co-chaperone ERdj4 (DnaJ B9, r 0.57). Tag-seq measures relative
  expression, so a rise of total mRNA cannot be seen; with at most about 1% of the reads in the
  program, normalisation shifts the other genes too little to show.
- **Not technical.** The mfp mean does not follow library size (Spearman 0.04), genes detected
  (0.02) or RNA yield (partial r -0.05, p 0.73). It is foot PC1 of `05` step 14's day-3 PCA
  (Spearman 0.92).
- **Not tied to threads or strength.** Within arms the mfp score goes with neither the threads
  made in the 3 days (partial r -0.12, p 0.44), day-3 mean peak force (0.15, p 0.35) nor plaque
  area (0.10, p 0.54), each adjusted for the animal's baseline where it has one.
- **Day 0.** Every lab-control foot expresses the mfp genes (median VST 10.8, 6.4 to 12.8), above
  most mfp-on animals at day 3. Their feet were dissected differently, so the comparison is not a
  measure of a change during the experiment.

What the state is cannot be settled here: an animal secreting new plaque material, or a foot
piece holding more of the plaque and cuticle glands (the pieces were cut by hand, and RNA yield
measures the size of the piece, not its make-up). The day-3 thread counts do not favour the
first: mfp-on animals made no more threads (median 9 against 7; within arms, on against off,
p 0.60).

Gill: 47 day-3 animals, 179 byssal genes, 157 expressed; one mfp gene passes the filter, so the mfp
tests are skipped. The thread core genes rise and fall together (mean within-arm r 0.50,
p 0.0005); the other groups 0.01 to 0.07.
