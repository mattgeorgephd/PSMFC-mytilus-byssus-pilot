# QC_plots

One JPEG per tensometer trace, in a subfolder per source folder of
`../../../02_data/tensometer_output/`, named after the trace file so retries never overwrite each
other. Each shows the raw force against time (blue) with the loess smoother used to find the
peak (red). Written by `../../../01_code/02_extract_tensometer_data.Rmd`; use them with
`note` and `failure` in the thread summary to spot bad pulls. Each subfolder's README names
the few plots left from an earlier run that have no trace file.

The plots are not committed (`.gitignore`, since 2026-10-05): the script writes all 383 on
every run, and they come out byte-different on every machine, so each run from another
computer added about 8 MB to the history. Run script 02 (or the folder's runner) to see them.
Committed here are the READMEs and the five plots left from the earlier run, which no script
writes and which could not be rebuilt.
