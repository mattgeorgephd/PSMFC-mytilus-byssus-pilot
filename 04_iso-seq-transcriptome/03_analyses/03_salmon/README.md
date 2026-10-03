# 03_salmon

Written by `../../01_code/03_salmon_quant.Rmd`: the 131 trimmed Tag-seq read files quantified
against the 411,251 Iso-Seq isoforms with salmon 1.10.3 (selective alignment, `-l A`,
`--noLengthCorrection`, no decoys; index k = 31), and the isoforms summed per feature of step
02's map with tximport (`countsFromAbundance = "no"`).

| File | Contents |
|---|---|
| `salmon_mapping_summary.csv` | per library: reads processed, reads mapped, percent mapped, the library type salmon detected |
| `gene_counts.csv.gz` | features (rows: gene keys, `novel:<locus>`, `mito:<isoform>`, `isoseq:<isoform>`) x libraries, summed read counts (EM-shared reads are fractional) |
| `read_classes_by_library.csv` | per library, the assigned reads by what their isoforms are: an annotated gene with the isoform ending within 100 bases of the gene's annotated 3' end (`gene_within_annotation`) or further out (`gene_3prime_extended`), a novel locus, an unassigned isoform, or mitochondrial; as reads and as a share (`pct_`) |
| `RUN_provenance.txt` | settings, salmon version, read source, input MD5s |
| `index/`, `index.log`, `quant/<library>/` | the salmon index and per-library output (git-ignored; the index takes about 15 minutes, each library 3 to 5 minutes on four cores) |

## Result (2026-10-02)

- Every library was detected as stranded, sense (`SF`), as 3' Tag-seq should be.
- salmon assigned a median 66.2% of each library's reads to the isoforms (53.6 to 77.8%;
  foot 70.3%, gill 62.3%), close to HISAT2's alignment rate to the genome for the same
  libraries (foot 70.1%, gill 63.7%; `../04_isoseq-de/FIG_mapping_rates.png`).
- The count matrix has 73,557 features, 65,041 of them with reads. salmon collapsed 27 isoforms
  that were exact duplicates of others.
- Where the assigned reads go (median share per library): 52.3% to isoforms that end within
  their gene's annotation, 18.4% to isoforms of annotated genes that run more than 100 bases
  past the gene's last annotated exon, 14.5% to novel loci, 7.0% to unassigned isoforms and
  7.4% to mitochondrial features. Tag-seq reads sit at the 3' end, so many reads of the
  extended isoforms, and all reads of novel loci, probably fall where the genome branch, which
  counts annotated exons, does not count them.
- The two libraries removed at QC stand out: T051G has 64% and T051F 49% of its assigned reads
  on mitochondrial features, against 4.9 to 10.4% in every other library.
