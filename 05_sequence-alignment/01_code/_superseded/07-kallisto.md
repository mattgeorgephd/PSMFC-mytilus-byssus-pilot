Kallisto Pseudo-alignment
================
Steven Roberts
14 December, 2023

- <a href="#1-download-tag-seq-data" id="toc-1-download-tag-seq-data">1
  Download tag-seq data</a>
- <a href="#2-get-transcriptome" id="toc-2-get-transcriptome">2 get
  transcriptome</a>
- <a href="#3-create-indices-for-transcriptome-with-kallisto"
  id="toc-3-create-indices-for-transcriptome-with-kallisto">3 Create
  indices for transcriptome with Kallisto</a>
- <a href="#4-kallisto-quantification"
  id="toc-4-kallisto-quantification">4 Kallisto quantification</a>
- <a href="#5-creating-count-matrix" id="toc-5-creating-count-matrix">5
  creating count matrix</a>
  - <a href="#501-note-sample-treatments-are-not-correct"
    id="toc-501-note-sample-treatments-are-not-correct">5.0.1 Note sample
    treatments are not correct.</a>

# 1 Download tag-seq data

``` bash
wget -r \
--no-directories --no-parent \
-P ../data \
-A .fastq.gz https://gannet.fish.washington.edu/panopea/PSMFC-mytilus-byssus-pilot/20220405-tagseq/ \
--no-check-certificate
```

# 2 get transcriptome

``` bash
cd ../data
curl -O https://owl.fish.washington.edu/halfshell/genomic-databank/Mtros-hq_transcripts.fasta
```

# 3 Create indices for transcriptome with Kallisto

``` bash
# Build index
/home/shared/kallisto_linux-v0.50.1/kallisto \
index \
-t 20 \
-i ../data/Mtros-hq_transcripts.index \
../data/Mtros-hq_transcripts.fasta
  
```

# 4 Kallisto quantification

``` bash
# Set the paths
DATA_DIRECTORY="../data"
KALLISTO_INDEX="../data/Mtros-hq_transcripts.index"
OUTPUT_DIRECTORY="../analyses/07-kallisto"


# Iterate over all .fq.gz files in the data directory
for FILE in "$DATA_DIRECTORY"/*.fastq.gz; do
    # Extract the base name of the file for naming the output folder
    BASENAME=$(basename "$FILE" _R1_001.fastq.gz)

    # Create output directory for this sample
    SAMPLE_OUTPUT="$OUTPUT_DIRECTORY/$BASENAME"
    mkdir -p "$SAMPLE_OUTPUT"

    # Run Kallisto quantification
    /home/shared/kallisto_linux-v0.50.1/kallisto quant -i "$KALLISTO_INDEX" -o "$SAMPLE_OUTPUT" \
        --single -t 20 -l 65 -s 2 "$FILE"
done

echo "Kallisto quantification complete."
```

# 5 creating count matrix

``` bash
perl /home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl \
--est_method kallisto \
    --gene_trans_map none \
    --out_prefix ../analyses/07-kallisto/test \
    --name_sample_by_basedir \
     ../analyses/07-kallisto/*/abundance.tsv
```

``` bash
head ../analyses/07-kallisto/test.isoform.counts.matrix
```

    ##  T001F_S100_L001 T001F_S100_L002 T001FX_S219_L001    T001FX_S219_L002    T001G_S112_L001 T001G_S112_L002 T002F_S101_L001 T002F_S101_L002 T002FX_S220_L001    T002FX_S220_L002    T002G_S113_L001 T002G_S113_L002 T003F_S102_L001 T003F_S102_L002 T003FX_S221_L001    T003FX_S221_L002    T003G_S114_L001 T003G_S114_L002 T004F_S103_L001 T004F_S103_L002 T004FX_S222_L001    T004FX_S222_L002    T004G_S115_L001 T004G_S115_L002 T005F_S104_L001 T005F_S104_L002 T005FX_S223_L001    T005FX_S223_L002    T005G_S116_L001 T005G_S116_L002 T006F_S105_L001 T006F_S105_L002 T006FX_S224_L001    T006FX_S224_L002    T006G_S117_L001 T006G_S117_L002 T007F_S106_L001 T007F_S106_L002 T007FX_S225_L001    T007FX_S225_L002    T007G_S118_L001 T007G_S118_L002 T008F_S107_L001 T008F_S107_L002 T008FX_S226_L001    T008FX_S226_L002    T008G_S119_L001 T008G_S119_L002 T009F_S108_L001 T009F_S108_L002 T009FX_S227_L001    T009FX_S227_L002    T009G_S120_L001 T009G_S120_L002 T010F_S109_L001 T010F_S109_L002 T010FX_S228_L001    T010FX_S228_L002    T010G_S121_L001 T010G_S121_L002 T011F_S110_L001 T011F_S110_L002 T011FX_S229_L001    T011FX_S229_L002    T011G_S122_L001 T011G_S122_L002 T012F_S111_L001 T012F_S111_L002
    ## MT_Pool_HQ_transcript/247050 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   14  18.039  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   1   0   0   0   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## MT_Pool_HQ_transcript/307636 0   0   0   0   1.0006  1   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   2.00661 0   0   0   0   0   0   1   0   1   0   0   1   0   0   0   0   0   1   1   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1.80716 1   0   0
    ## MT_Pool_HQ_transcript/1516   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## MT_Pool_HQ_transcript/200050 0   0   0   3.60284 4.28138 0   0   0   0   0   0   0   0   0   0   2.41717 1.07659 0   0   0   0   0   1.01944 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1.00726 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1.06984 0   0   0   0   0   0   0   0   0   0
    ## MT_Pool_HQ_transcript/42823  0   0   0   0   0   0   0   0   1.19923 0   0   0   0   0   0   0   0   1.46628 0   0   0   0   0   0   1.14306e-08 1   6.49272 2.45278e-08 11.3617 7.20403 0   0   0   0   0   0   0   0   0   0   0   0   0   2.15106e-08 0   0   0   0   0.00269595  0   0   0   2.178   3.07395 0   0   2   0   0.00268361  0   0   0   0   0   0   0   0   0
    ## MT_Pool_HQ_transcript/37520  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   4.07384 0   0   0   0   0   0   0
    ## MT_Pool_HQ_transcript/121901 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    ## MT_Pool_HQ_transcript/315546 0   0   0   3   2   1   4.79534 2   0   1   2   5.69027 0   0   1   3   0   5   0   1   1.86964 0.681234    4   3   0   0   0   0   6.94781 14.6353 0   0   1.31665 2.25    5.88393 7.95696 0   0   2   1   2   2.35035 0   1   0   0   10  7.78067 0   0   0   0   2.01093 4.93093 2   0   0   0   7.77766 10.5634 0   1   1   2   2   2.71926 1   1
    ## MT_Pool_HQ_transcript/162009 0   1.05177 0   0.563655    5.25029 0   2.42551 0   6.97861 3.90924 0   0   0   0   0   0   0   0   0   0   0   6.73435 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   4.3825  6.42377 4.39924 0   0

``` r
countmatrix <- read.delim("../analyses/07-kallisto/test.isoform.counts.matrix", header = TRUE, sep = '\t')
rownames(countmatrix) <- countmatrix$X
countmatrix <- countmatrix[,-1]
head(countmatrix)
```

    ##                              T001F_S100_L001 T001F_S100_L002 T001FX_S219_L001
    ## MT_Pool_HQ_transcript/247050               0               0                0
    ## MT_Pool_HQ_transcript/307636               0               0                0
    ## MT_Pool_HQ_transcript/1516                 0               0                0
    ## MT_Pool_HQ_transcript/200050               0               0                0
    ## MT_Pool_HQ_transcript/42823                0               0                0
    ## MT_Pool_HQ_transcript/37520                0               0                0
    ##                              T001FX_S219_L002 T001G_S112_L001 T001G_S112_L002
    ## MT_Pool_HQ_transcript/247050          0.00000         0.00000               0
    ## MT_Pool_HQ_transcript/307636          0.00000         1.00060               1
    ## MT_Pool_HQ_transcript/1516            0.00000         0.00000               0
    ## MT_Pool_HQ_transcript/200050          3.60284         4.28138               0
    ## MT_Pool_HQ_transcript/42823           0.00000         0.00000               0
    ## MT_Pool_HQ_transcript/37520           0.00000         0.00000               0
    ##                              T002F_S101_L001 T002F_S101_L002 T002FX_S220_L001
    ## MT_Pool_HQ_transcript/247050               0               0          0.00000
    ## MT_Pool_HQ_transcript/307636               0               0          0.00000
    ## MT_Pool_HQ_transcript/1516                 0               0          0.00000
    ## MT_Pool_HQ_transcript/200050               0               0          0.00000
    ## MT_Pool_HQ_transcript/42823                0               0          1.19923
    ## MT_Pool_HQ_transcript/37520                0               0          0.00000
    ##                              T002FX_S220_L002 T002G_S113_L001 T002G_S113_L002
    ## MT_Pool_HQ_transcript/247050                0               0               0
    ## MT_Pool_HQ_transcript/307636                0               0               0
    ## MT_Pool_HQ_transcript/1516                  0               0               0
    ## MT_Pool_HQ_transcript/200050                0               0               0
    ## MT_Pool_HQ_transcript/42823                 0               0               0
    ## MT_Pool_HQ_transcript/37520                 0               0               0
    ##                              T003F_S102_L001 T003F_S102_L002 T003FX_S221_L001
    ## MT_Pool_HQ_transcript/247050               0               0                0
    ## MT_Pool_HQ_transcript/307636               0               0                0
    ## MT_Pool_HQ_transcript/1516                 0               0                0
    ## MT_Pool_HQ_transcript/200050               0               0                0
    ## MT_Pool_HQ_transcript/42823                0               0                0
    ## MT_Pool_HQ_transcript/37520                0               0                0
    ##                              T003FX_S221_L002 T003G_S114_L001 T003G_S114_L002
    ## MT_Pool_HQ_transcript/247050          0.00000        14.00000        18.03900
    ## MT_Pool_HQ_transcript/307636          0.00000         0.00000         1.00000
    ## MT_Pool_HQ_transcript/1516            0.00000         0.00000         0.00000
    ## MT_Pool_HQ_transcript/200050          2.41717         1.07659         0.00000
    ## MT_Pool_HQ_transcript/42823           0.00000         0.00000         1.46628
    ## MT_Pool_HQ_transcript/37520           0.00000         0.00000         0.00000
    ##                              T004F_S103_L001 T004F_S103_L002 T004FX_S222_L001
    ## MT_Pool_HQ_transcript/247050               0               0                0
    ## MT_Pool_HQ_transcript/307636               0               0                0
    ## MT_Pool_HQ_transcript/1516                 0               0                0
    ## MT_Pool_HQ_transcript/200050               0               0                0
    ## MT_Pool_HQ_transcript/42823                0               0                0
    ## MT_Pool_HQ_transcript/37520                0               0                0
    ##                              T004FX_S222_L002 T004G_S115_L001 T004G_S115_L002
    ## MT_Pool_HQ_transcript/247050                0         0.00000         0.00000
    ## MT_Pool_HQ_transcript/307636                0         0.00000         2.00661
    ## MT_Pool_HQ_transcript/1516                  0         0.00000         0.00000
    ## MT_Pool_HQ_transcript/200050                0         1.01944         0.00000
    ## MT_Pool_HQ_transcript/42823                 0         0.00000         0.00000
    ## MT_Pool_HQ_transcript/37520                 0         0.00000         0.00000
    ##                              T005F_S104_L001 T005F_S104_L002 T005FX_S223_L001
    ## MT_Pool_HQ_transcript/247050     0.00000e+00               0          0.00000
    ## MT_Pool_HQ_transcript/307636     0.00000e+00               0          0.00000
    ## MT_Pool_HQ_transcript/1516       0.00000e+00               0          0.00000
    ## MT_Pool_HQ_transcript/200050     0.00000e+00               0          0.00000
    ## MT_Pool_HQ_transcript/42823      1.14306e-08               1          6.49272
    ## MT_Pool_HQ_transcript/37520      0.00000e+00               0          0.00000
    ##                              T005FX_S223_L002 T005G_S116_L001 T005G_S116_L002
    ## MT_Pool_HQ_transcript/247050      0.00000e+00          0.0000         0.00000
    ## MT_Pool_HQ_transcript/307636      0.00000e+00          0.0000         0.00000
    ## MT_Pool_HQ_transcript/1516        0.00000e+00          0.0000         0.00000
    ## MT_Pool_HQ_transcript/200050      0.00000e+00          0.0000         0.00000
    ## MT_Pool_HQ_transcript/42823       2.45278e-08         11.3617         7.20403
    ## MT_Pool_HQ_transcript/37520       0.00000e+00          0.0000         0.00000
    ##                              T006F_S105_L001 T006F_S105_L002 T006FX_S224_L001
    ## MT_Pool_HQ_transcript/247050               0               0                0
    ## MT_Pool_HQ_transcript/307636               1               0                1
    ## MT_Pool_HQ_transcript/1516                 0               0                0
    ## MT_Pool_HQ_transcript/200050               0               0                0
    ## MT_Pool_HQ_transcript/42823                0               0                0
    ## MT_Pool_HQ_transcript/37520                0               0                0
    ##                              T006FX_S224_L002 T006G_S117_L001 T006G_S117_L002
    ## MT_Pool_HQ_transcript/247050                0               0               0
    ## MT_Pool_HQ_transcript/307636                0               0               1
    ## MT_Pool_HQ_transcript/1516                  0               0               0
    ## MT_Pool_HQ_transcript/200050                0               0               0
    ## MT_Pool_HQ_transcript/42823                 0               0               0
    ## MT_Pool_HQ_transcript/37520                 0               0               0
    ##                              T007F_S106_L001 T007F_S106_L002 T007FX_S225_L001
    ## MT_Pool_HQ_transcript/247050               0               0                1
    ## MT_Pool_HQ_transcript/307636               0               0                0
    ## MT_Pool_HQ_transcript/1516                 0               0                0
    ## MT_Pool_HQ_transcript/200050               0               0                0
    ## MT_Pool_HQ_transcript/42823                0               0                0
    ## MT_Pool_HQ_transcript/37520                0               0                0
    ##                              T007FX_S225_L002 T007G_S118_L001 T007G_S118_L002
    ## MT_Pool_HQ_transcript/247050                0               0         0.00000
    ## MT_Pool_HQ_transcript/307636                0               0         1.00000
    ## MT_Pool_HQ_transcript/1516                  0               0         0.00000
    ## MT_Pool_HQ_transcript/200050                0               0         1.00726
    ## MT_Pool_HQ_transcript/42823                 0               0         0.00000
    ## MT_Pool_HQ_transcript/37520                 0               0         0.00000
    ##                              T008F_S107_L001 T008F_S107_L002 T008FX_S226_L001
    ## MT_Pool_HQ_transcript/247050               0     0.00000e+00                0
    ## MT_Pool_HQ_transcript/307636               1     0.00000e+00                0
    ## MT_Pool_HQ_transcript/1516                 0     0.00000e+00                0
    ## MT_Pool_HQ_transcript/200050               0     0.00000e+00                0
    ## MT_Pool_HQ_transcript/42823                0     2.15106e-08                0
    ## MT_Pool_HQ_transcript/37520                0     0.00000e+00                0
    ##                              T008FX_S226_L002 T008G_S119_L001 T008G_S119_L002
    ## MT_Pool_HQ_transcript/247050                0               0               0
    ## MT_Pool_HQ_transcript/307636                0               0               0
    ## MT_Pool_HQ_transcript/1516                  0               0               0
    ## MT_Pool_HQ_transcript/200050                0               0               0
    ## MT_Pool_HQ_transcript/42823                 0               0               0
    ## MT_Pool_HQ_transcript/37520                 0               0               0
    ##                              T009F_S108_L001 T009F_S108_L002 T009FX_S227_L001
    ## MT_Pool_HQ_transcript/247050      1.00000000               0                0
    ## MT_Pool_HQ_transcript/307636      1.00000000               1                0
    ## MT_Pool_HQ_transcript/1516        0.00000000               0                0
    ## MT_Pool_HQ_transcript/200050      0.00000000               0                0
    ## MT_Pool_HQ_transcript/42823       0.00269595               0                0
    ## MT_Pool_HQ_transcript/37520       0.00000000               0                0
    ##                              T009FX_S227_L002 T009G_S120_L001 T009G_S120_L002
    ## MT_Pool_HQ_transcript/247050                0           2.000         0.00000
    ## MT_Pool_HQ_transcript/307636                0           0.000         0.00000
    ## MT_Pool_HQ_transcript/1516                  0           0.000         0.00000
    ## MT_Pool_HQ_transcript/200050                0           0.000         0.00000
    ## MT_Pool_HQ_transcript/42823                 0           2.178         3.07395
    ## MT_Pool_HQ_transcript/37520                 0           0.000         0.00000
    ##                              T010F_S109_L001 T010F_S109_L002 T010FX_S228_L001
    ## MT_Pool_HQ_transcript/247050               0               0                0
    ## MT_Pool_HQ_transcript/307636               0               0                0
    ## MT_Pool_HQ_transcript/1516                 0               0                0
    ## MT_Pool_HQ_transcript/200050               0               0                0
    ## MT_Pool_HQ_transcript/42823                0               0                2
    ## MT_Pool_HQ_transcript/37520                0               0                0
    ##                              T010FX_S228_L002 T010G_S121_L001 T010G_S121_L002
    ## MT_Pool_HQ_transcript/247050          0.00000      0.00000000               0
    ## MT_Pool_HQ_transcript/307636          0.00000      0.00000000               0
    ## MT_Pool_HQ_transcript/1516            0.00000      0.00000000               0
    ## MT_Pool_HQ_transcript/200050          1.06984      0.00000000               0
    ## MT_Pool_HQ_transcript/42823           0.00000      0.00268361               0
    ## MT_Pool_HQ_transcript/37520           0.00000      0.00000000               0
    ##                              T011F_S110_L001 T011F_S110_L002 T011FX_S229_L001
    ## MT_Pool_HQ_transcript/247050         0.00000               0                0
    ## MT_Pool_HQ_transcript/307636         0.00000               0                0
    ## MT_Pool_HQ_transcript/1516           0.00000               0                0
    ## MT_Pool_HQ_transcript/200050         0.00000               0                0
    ## MT_Pool_HQ_transcript/42823          0.00000               0                0
    ## MT_Pool_HQ_transcript/37520          4.07384               0                0
    ##                              T011FX_S229_L002 T011G_S122_L001 T011G_S122_L002
    ## MT_Pool_HQ_transcript/247050                0         0.00000               0
    ## MT_Pool_HQ_transcript/307636                0         1.80716               1
    ## MT_Pool_HQ_transcript/1516                  0         0.00000               0
    ## MT_Pool_HQ_transcript/200050                0         0.00000               0
    ## MT_Pool_HQ_transcript/42823                 0         0.00000               0
    ## MT_Pool_HQ_transcript/37520                 0         0.00000               0
    ##                              T012F_S111_L001 T012F_S111_L002
    ## MT_Pool_HQ_transcript/247050               0               0
    ## MT_Pool_HQ_transcript/307636               0               0
    ## MT_Pool_HQ_transcript/1516                 0               0
    ## MT_Pool_HQ_transcript/200050               0               0
    ## MT_Pool_HQ_transcript/42823                0               0
    ## MT_Pool_HQ_transcript/37520                0               0

``` r
countmatrix <- round(countmatrix, 0)
str(countmatrix)
```

    ## 'data.frame':    411251 obs. of  68 variables:
    ##  $ T001F_S100_L001 : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T001F_S100_L002 : num  0 0 0 0 0 0 0 0 1 0 ...
    ##  $ T001FX_S219_L001: num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T001FX_S219_L002: num  0 0 0 4 0 0 0 3 1 0 ...
    ##  $ T001G_S112_L001 : num  0 1 0 4 0 0 0 2 5 0 ...
    ##  $ T001G_S112_L002 : num  0 1 0 0 0 0 0 1 0 0 ...
    ##  $ T002F_S101_L001 : num  0 0 0 0 0 0 0 5 2 0 ...
    ##  $ T002F_S101_L002 : num  0 0 0 0 0 0 0 2 0 0 ...
    ##  $ T002FX_S220_L001: num  0 0 0 0 1 0 0 0 7 2 ...
    ##  $ T002FX_S220_L002: num  0 0 0 0 0 0 0 1 4 3 ...
    ##  $ T002G_S113_L001 : num  0 0 0 0 0 0 0 2 0 0 ...
    ##  $ T002G_S113_L002 : num  0 0 0 0 0 0 0 6 0 0 ...
    ##  $ T003F_S102_L001 : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T003F_S102_L002 : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T003FX_S221_L001: num  0 0 0 0 0 0 0 1 0 0 ...
    ##  $ T003FX_S221_L002: num  0 0 0 2 0 0 0 3 0 0 ...
    ##  $ T003G_S114_L001 : num  14 0 0 1 0 0 0 0 0 0 ...
    ##  $ T003G_S114_L002 : num  18 1 0 0 1 0 0 5 0 2 ...
    ##  $ T004F_S103_L001 : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T004F_S103_L002 : num  0 0 0 0 0 0 0 1 0 0 ...
    ##  $ T004FX_S222_L001: num  0 0 0 0 0 0 0 2 0 0 ...
    ##  $ T004FX_S222_L002: num  0 0 0 0 0 0 0 1 7 0 ...
    ##  $ T004G_S115_L001 : num  0 0 0 1 0 0 0 4 0 0 ...
    ##  $ T004G_S115_L002 : num  0 2 0 0 0 0 0 3 0 1 ...
    ##  $ T005F_S104_L001 : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T005F_S104_L002 : num  0 0 0 0 1 0 0 0 0 0 ...
    ##  $ T005FX_S223_L001: num  0 0 0 0 6 0 0 0 0 0 ...
    ##  $ T005FX_S223_L002: num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T005G_S116_L001 : num  0 0 0 0 11 0 0 7 0 0 ...
    ##  $ T005G_S116_L002 : num  0 0 0 0 7 0 0 15 0 0 ...
    ##  $ T006F_S105_L001 : num  0 1 0 0 0 0 0 0 0 0 ...
    ##  $ T006F_S105_L002 : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T006FX_S224_L001: num  0 1 0 0 0 0 0 1 0 0 ...
    ##  $ T006FX_S224_L002: num  0 0 0 0 0 0 0 2 0 0 ...
    ##  $ T006G_S117_L001 : num  0 0 0 0 0 0 0 6 0 0 ...
    ##  $ T006G_S117_L002 : num  0 1 0 0 0 0 0 8 0 0 ...
    ##  $ T007F_S106_L001 : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T007F_S106_L002 : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T007FX_S225_L001: num  1 0 0 0 0 0 0 2 0 0 ...
    ##  $ T007FX_S225_L002: num  0 0 0 0 0 0 0 1 0 0 ...
    ##  $ T007G_S118_L001 : num  0 0 0 0 0 0 0 2 0 6 ...
    ##  $ T007G_S118_L002 : num  0 1 0 1 0 0 0 2 0 0 ...
    ##  $ T008F_S107_L001 : num  0 1 0 0 0 0 0 0 0 0 ...
    ##  $ T008F_S107_L002 : num  0 0 0 0 0 0 0 1 0 0 ...
    ##  $ T008FX_S226_L001: num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T008FX_S226_L002: num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T008G_S119_L001 : num  0 0 0 0 0 0 0 10 0 0 ...
    ##  $ T008G_S119_L002 : num  0 0 0 0 0 0 0 8 0 0 ...
    ##  $ T009F_S108_L001 : num  1 1 0 0 0 0 0 0 0 0 ...
    ##  $ T009F_S108_L002 : num  0 1 0 0 0 0 0 0 0 0 ...
    ##  $ T009FX_S227_L001: num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T009FX_S227_L002: num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ T009G_S120_L001 : num  2 0 0 0 2 0 0 2 0 2 ...
    ##  $ T009G_S120_L002 : num  0 0 0 0 3 0 0 5 0 0 ...
    ##  $ T010F_S109_L001 : num  0 0 0 0 0 0 0 2 0 0 ...
    ##  $ T010F_S109_L002 : num  0 0 0 0 0 0 0 0 0 1 ...
    ##  $ T010FX_S228_L001: num  0 0 0 0 2 0 0 0 0 1 ...
    ##  $ T010FX_S228_L002: num  0 0 0 1 0 0 0 0 0 0 ...
    ##  $ T010G_S121_L001 : num  0 0 0 0 0 0 0 8 0 0 ...
    ##  $ T010G_S121_L002 : num  0 0 0 0 0 0 0 11 0 0 ...
    ##  $ T011F_S110_L001 : num  0 0 0 0 0 4 0 0 0 0 ...
    ##  $ T011F_S110_L002 : num  0 0 0 0 0 0 0 1 0 0 ...
    ##  $ T011FX_S229_L001: num  0 0 0 0 0 0 0 1 0 0 ...
    ##  $ T011FX_S229_L002: num  0 0 0 0 0 0 0 2 4 0 ...
    ##  $ T011G_S122_L001 : num  0 2 0 0 0 0 0 2 6 5 ...
    ##  $ T011G_S122_L002 : num  0 1 0 0 0 0 0 3 4 0 ...
    ##  $ T012F_S111_L001 : num  0 0 0 0 0 0 0 1 0 5 ...
    ##  $ T012F_S111_L002 : num  0 0 0 0 0 0 0 1 0 0 ...

``` r
deseq2.colData <- data.frame(condition=factor(c(rep("control", 34), rep("desicated", 34))), 
                             type=factor(rep("single-read", 68)))
rownames(deseq2.colData) <- colnames(data)
deseq2.dds <- DESeqDataSetFromMatrix(countData = countmatrix,
                                     colData = deseq2.colData, 
                                     design = ~ condition)
```

``` r
deseq2.dds <- DESeq(deseq2.dds)
deseq2.res <- results(deseq2.dds)
deseq2.res <- deseq2.res[order(rownames(deseq2.res)), ]
```

### 5.0.1 Note sample treatments are not correct.

``` r
vsd <- vst(deseq2.dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")
```
