# 05_rrvgo

topGO's enriched terms (p < 0.01) of each run, per ontology, grouped into clusters of similar
terms (Wang similarity computed with `GOSemSim::goSim()`, cut at 0.7), each named after its
most significant term; written by `../../01_code/05_rrvgo.Rmd`. Replaces the REVIGO web
submissions. A term the GO graph cannot place (all similarities missing) is kept as its own
cluster, with a message in the report, rather than dropped.

| File | Contents |
|---|---|
| `rrvgo_reduced_terms.csv` | every enriched term with its cluster, parent, score (-log10 p), annotated genes (`size`), DEGs in term and `p_weight01` |
| `rrvgo_parents.csv` | one row per parent term, run and ontology: terms it stands for, smallest p, distinct DEGs and their IDs |
| `rrvgo_<TC,FG>_<BP,MF,CC>_parents.png` | the parent terms of every run side by side, one figure per family and ontology |
