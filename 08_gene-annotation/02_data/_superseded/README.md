# 02_data/_superseded

Retired inputs, kept as records. No current script reads them.

| File | What it was | Replaced by |
|---|---|---|
| `goslim_generic_2023-07-27.obo` | the generic GO slim of GO release 2023-07-27 (141 terms; 72 biological process), from `geneontology/go-ontology` `src/ontology/go-edit.obo` at commit 17c29bb; `01_go_slims.Rmd` read it as `../goslim_generic.obo` with `GO.db` 3.18.0 until 2026-10-03 (moved here unchanged) | `../goslim_generic.obo`, the slim of GO release 2026-01-23, the release of `GO.db` 3.23.1. `../../01_code/_derive_goslim.R` rebuilds this file byte for byte from the 2023 commit |
