# 03_analyses

| Subfolder | Produced by | Contents |
|-----------|-------------|----------|
| `DEG_lists/` | `02_5_DESeq_*`, `03-*_Shrinkage_filtration`, `04-File_joining` | Significant-DEG and count/treatment tables per tissue x contrast, under `Foot/`, `Gill/`, and `GOterms_genome/` |

`DEG_lists/` is an intermediate: the DESeq scripts write it and the joining, venn, volcano,
counts, and top-genes scripts read it back. Enrichment (in `../../enrichment/`) also reads it.
