# QC_plots

One JPEG per tensometer trace, in a subfolder per source folder of
`../../../02_data/tensometer_output/`, named after the trace file so retries never overwrite each
other. Each shows the raw force against time (blue) with the loess smoother used to find the
peak (red). Written by `../../../01_code/02_extract_tensometer_data.Rmd`; use them with
`note` and `failure` in the thread summary to spot bad pulls. Each subfolder's README names
the few plots left from an earlier run that have no trace file.
