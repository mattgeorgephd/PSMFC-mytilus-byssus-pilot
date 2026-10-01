# 05_rrvgo

topGO's enriched BP terms (p < 0.01) of each run grouped into clusters of similar terms
(Wang similarity, cut at 0.7), each named after its most significant term; written by
`../../01_code/05_rrvgo.Rmd`. Replaces the REVIGO web submissions.

| File | Contents |
|---|---|
| `rrvgo_reduced_terms.csv` | every enriched BP term with its cluster, parent, score (-log10 p), annotated genes (`size`) and DEGs in term |
| `rrvgo_parents.csv` | one row per parent term and run: terms it stands for, smallest p, distinct DEGs, mitochondrially encoded protein DEGs |
| `rrvgo_<TC,LC,FG>_BP_parents.png` | the parent terms of every run side by side |
