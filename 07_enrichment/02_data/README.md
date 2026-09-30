# 02_data

No enrichment inputs are stored here. The enrichment scripts read their inputs cross-folder:

- DEG lists from `../../differential-expression/03_analyses/DEG_lists/`
- gene-to-GO mapping from `../../blast/03_analyses/genome-foot/`

Phase 4 will wire these reads through a `repo_root` pointer.
