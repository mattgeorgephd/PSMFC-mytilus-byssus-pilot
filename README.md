# claude/blast-2026-parts (scratch; do not merge)

The finished parts of the genome blastx of pull request #3 (GO-db-2026), copied here as the
search runs so that a new session can resume it if the cloud container is lost. Nothing in
`main` or in #3 reads this branch; delete it once #3 carries the full blastx table
(`03_blast/03_analyses/genome-foot-sprot2026_03/`).

## The search

`03_blast/01_code/01_genome_blast.Rmd` (branch `claude/go-db-2026`), knit with:

    run: true, sprot_release: "2026_03",
    foot_fasta: "uniprotkb_mytilus_foot_2026_03_byssal.fasta",
    res_folder: "genome-foot-sprot2026_03", parts: 200, jobs: 1, threads: 4

- query: `cds_from_genomic.fasta` of RefSeq GCF_036588685.1 (53,269 CDSs), MD5 a4df841e16db68b41719d564986f706b
- database: UniProtKB/Swiss-Prot release 2026_03 (575,748 entries, UniProt `current_release` on
  2026-10-03) plus the 241 proteins of `03_blast/02_data/uniprotkb_mytilus_foot_2026_03_byssal.fasta`
  (MD5 7904e3f596a80d528f9e900fdc4f8b79): 575,989 sequences, 209,093,815 residues
- BLAST+ 2.15.0 blastx `-evalue 1E-20 -max_target_seqs 1 -outfmt 6`, through
  `03_blast/01_code/_blastx_parts.sh`

`parts/part_NNNN.tab` is the blastx output of part NNNN of 200 (the CDS FASTA split in order
into 200 files of about equal numbers of sequences, as `_blastx_parts.sh` splits it).

## Resuming

Copy `parts/part_*.tab` into `<work_dir>/output/genome_n_foot_blastx_r2026_03.tab.parts/` and
knit the script with the parameters above (or run `_blastx_parts.sh` with the same arguments).
The helper splits the CDS FASTA into the same 200 parts and searches only those without a
`.tab`; once all 200 are there it joins them in order. While UniProt serves release 2026_03
(until its next release, expected late October 2026) the database can be rebuilt from
`current_release`; after that, from `previous_releases/release-2026_03`.
