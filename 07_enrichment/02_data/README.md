# 02_data

No enrichment inputs are stored here. The enrichment scripts read their inputs cross-folder
through `../01_code/_paths.R`:

- DEG tables from `../../06_differential-expression/03_analyses/DEG_lists/GOterms_genome/`
- gene count matrices from `../../06_differential-expression/02_data/`
- gene-to-UniProt / GO mapping from `../../03_blast/03_analyses/genome-foot/`
