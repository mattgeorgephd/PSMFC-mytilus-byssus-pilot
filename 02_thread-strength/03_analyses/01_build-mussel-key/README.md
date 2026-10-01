# 01_build-mussel-key

Written by `../../01_code/01_build_mussel_key.Rmd` from
`../../../01_mussel-measurements/mussel-size-measurements.xlsx` (sheet `data`).

`mussel-treatment-key.csv`: one row per mussel tag (`T037`, zero-padded) with `species`,
`mussel_trt` (the arm the animal was assigned to: control, OA, OW, DO, lab_control,
desiccation), `group`, `days_in_trt` and `rna_sequenced`. Read by script 02, which joins it to
every trace. CSV rather than xlsx so changes show in `git diff`. Moved here from `02_data/`:
it is generated, so it belongs with the outputs.
