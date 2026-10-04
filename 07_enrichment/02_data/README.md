# 02_data

No inputs are stored here. The scripts read, through `../01_code/_paths.R`:

- the contrasts and their apeglm tables from
  `../../05_differential-expression/03_analyses/DEG_lists/`;
- the genome-wide BLAST / UniProt / GO table `LOC_GO_list.txt` (the genome blast of 2026,
  with UniProt release 2026_03 records) from `../../03_blast/03_analyses/genome-foot-sprot2026_03/`;
- the reference transcript lengths in `t_data.ctab` from
  `../../04_sequence-alignment/03_analyses/hisat/`.
