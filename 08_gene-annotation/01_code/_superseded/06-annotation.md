Annotation
================
Steven Roberts
07 August, 2023

In this following chunk where the fasta file is downloaded the
[release](https://www.uniprot.org/help/release-statistics) is noted and
the file name is modified accordingly.

``` bash
cd ../data

curl -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

mv uniprot_sprot.fasta.gz uniprot_sprot_r2023_02.fasta.gz
gunzip -k uniprot_sprot_r2023_02.fasta.gz
```

A protein blast database is then made.

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../data/uniprot_sprot_r2023_02.fasta \
-dbtype prot \
-out ../analyses/blast/uniprot_sprot_r2023_02
```

``` bash
head -4 ../data/Mtros-hq_transcripts.fasta


grep ">" ../data/Mtros-hq_transcripts.fasta | head
```

    ## >MT_Pool_HQ_transcript/0
    ## AGTGAGAGACAAGCCGGCTTGTTACAAAACATCTCCCTTTTCCTTGTCAACTCTTGGAGT
    ## ACGGACGGGACCTGACTTTTGACCACCGGAACTGACTTTCGGATTGATTATTATAATAAT
    ## AGTCAATCAGAAGTGATCAAGACCAGCATTTTGTGCTGAAGGACCGAACAGGAGCGTCCA
    ## >MT_Pool_HQ_transcript/0
    ## >MT_Pool_HQ_transcript/1
    ## >MT_Pool_HQ_transcript/2
    ## >MT_Pool_HQ_transcript/3
    ## >MT_Pool_HQ_transcript/4
    ## >MT_Pool_HQ_transcript/5
    ## >MT_Pool_HQ_transcript/6
    ## >MT_Pool_HQ_transcript/7
    ## >MT_Pool_HQ_transcript/8
    ## >MT_Pool_HQ_transcript/10

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastx \
-query ../data/Mtros-hq_transcripts.fasta \
-db ../analyses/blast/uniprot_sprot_r2023_02 \
-out ../analyses/blast/Mtros-hq-uniprot_blastx.tab \
-evalue 1E-20 \
-num_threads 40 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6
```

Here is what the output file looks like, and at this point we want to
get the UniProt Accession number for each gene

``` bash
head -2 ../analyses/blast/Mtros-hq-uniprot_blastx.tab
```

    ## MT_Pool_HQ_transcript/0  sp|D3ZHA0|FLNC_RAT  42.357  2486    1247    48  1977    8993    278 2724    0.0 1647
    ## MT_Pool_HQ_transcript/1  sp|P16157|ANK1_HUMAN    31.783  774 524 1   6347    8668    17  786 5.29e-119   425

``` r
blast <- read.csv("../analyses/blast/Mtros-hq-uniprot_blastx.tab", sep = '\t', header = FALSE)
```

``` r
blastsp <- blast %>%
  mutate(SPID = str_extract(V2, "(?<=\\|)[^\\|]*(?=\\|)"))
```

Now lets just write out SPIDs.

``` r
blast %>%
  mutate(SPID = str_extract(V2, "(?<=\\|)[^\\|]*(?=\\|)")) %>%
  distinct(V1, SPID, .keep_all = TRUE) %>%
  select(SPID) %>%
  write.table(file = "../analyses/blast/SPID.txt", sep = "\t", row.names = FALSE, quote = FALSE
) 
```

With a list of matching Swiss-Prot IDs, (technically UniProt Accession
number) we can go back to <https://www.uniprot.org> and grab
corresponding GO terms. This can be done via a web or using Python API.

``` bash
python3 uniprot-retrieval.py ../analyses/blast/SPID.txt
```

``` bash
mv uniprot-retrieval.tsv.gz ../analyses/blast/uniprot-retrieval.tsv.gz

gunzip ../analyses/blast/uniprot-retrieval.tsv.gz
```

Web

``` r
idmap <- read.csv("../analyses/blast/idmapping_2023_08_07.tsv", sep = '\t', header = TRUE, row.names=NULL)
```

Joining

``` r
left_join(blastsp, idmap, by = c("SPID" = "Entry")) %>%
  select(V1, V2, V11, Protein.names, Gene.Ontology..biological.process., Gene.Ontology.IDs) %>%
write.table(file = "../analyses/blast/Mtros_GO.txt", sep = #"\t", row.names = FALSE, quote = FALSE
  )
```

``` bash
head ../analyses/blast/Mtros_GO.txt
```

    ## "V1" "V2" "V11" "Protein.names" "Gene.Ontology..biological.process." "Gene.Ontology.IDs"
    ## "1" "MT_Pool_HQ_transcript/0" "sp|D3ZHA0|FLNC_RAT" 0 "Filamin-C (FLN-C) (ABP-280-like protein) (ABP-L) (Actin-binding-like protein) (Filamin-2) (Gamma-filamin)" "muscle cell development [GO:0055001]; sarcomere organization [GO:0045214]" "GO:0005737; GO:0005856; GO:0008092; GO:0016528; GO:0030018; GO:0030506; GO:0042383; GO:0043034; GO:0045214; GO:0051015; GO:0055001"
    ## "2" "MT_Pool_HQ_transcript/1" "sp|P16157|ANK1_HUMAN" 5.29e-119 "Ankyrin-1 (ANK-1) (Ankyrin-R) (Erythrocyte ankyrin)" "cytoskeleton organization [GO:0007010]; endoplasmic reticulum to Golgi vesicle-mediated transport [GO:0006888]; exocytosis [GO:0006887]; maintenance of epithelial cell apical/basal polarity [GO:0045199]; protein localization to plasma membrane [GO:0072659]; signal transduction [GO:0007165]" "GO:0005198; GO:0005200; GO:0005829; GO:0005856; GO:0005886; GO:0006887; GO:0006888; GO:0007010; GO:0007165; GO:0008093; GO:0009898; GO:0014731; GO:0016323; GO:0016529; GO:0019899; GO:0019903; GO:0030507; GO:0031430; GO:0043005; GO:0044325; GO:0045199; GO:0051117; GO:0072659"
    ## "3" "MT_Pool_HQ_transcript/2" "sp|P10079|FBP1_STRPU" 5.68e-88 "Fibropellin-1 (Epidermal growth factor-related protein 1) (Fibropellin-I) (SpEGF I) (UEGF-1)" "" "GO:0005509; GO:0005615; GO:0009374; GO:0031410; GO:0032579"
    ## "4" "MT_Pool_HQ_transcript/3" "sp|P10079|FBP1_STRPU" 2.05e-88 "Fibropellin-1 (Epidermal growth factor-related protein 1) (Fibropellin-I) (SpEGF I) (UEGF-1)" "" "GO:0005509; GO:0005615; GO:0009374; GO:0031410; GO:0032579"
    ## "5" "MT_Pool_HQ_transcript/4" "sp|P12276|FAS_CHICK" 0 "Fatty acid synthase (EC 2.3.1.85) [Includes: [Acyl-carrier-protein] S-acetyltransferase (EC 2.3.1.38); [Acyl-carrier-protein] S-malonyltransferase (EC 2.3.1.39); 3-oxoacyl-[acyl-carrier-protein] synthase (EC 2.3.1.41); 3-oxoacyl-[acyl-carrier-protein] reductase (EC 1.1.1.100); 3-hydroxyacyl-[acyl-carrier-protein] dehydratase (EC 4.2.1.59); Enoyl-[acyl-carrier-protein] reductase (EC 1.3.1.39); Acyl-[acyl-carrier-protein] hydrolase (EC 3.1.2.14)]" "fatty acid biosynthetic process [GO:0006633]; lactate metabolic process [GO:0006089]; positive regulation of appetite [GO:0032100]" "GO:0003697; GO:0004312; GO:0004313; GO:0004314; GO:0004315; GO:0004316; GO:0004317; GO:0004320; GO:0005737; GO:0006089; GO:0006633; GO:0008659; GO:0008693; GO:0016295; GO:0016296; GO:0031177; GO:0032100; GO:0047117; GO:0047451"
    ## "6" "MT_Pool_HQ_transcript/5" "sp|Q6V0I7|FAT4_HUMAN" 4.12e-81 "Protocadherin Fat 4 (hFat4) (Cadherin family member 14) (FAT tumor suppressor homolog 4) (Fat-like cadherin protein FAT-J)" "cerebral cortex development [GO:0021987]; heterophilic cell-cell adhesion via plasma membrane cell adhesion molecules [GO:0007157]; hippo signaling [GO:0035329]; homophilic cell adhesion via plasma membrane adhesion molecules [GO:0007156]; neurogenesis [GO:0022008]" "GO:0005509; GO:0005886; GO:0007156; GO:0007157; GO:0021987; GO:0022008; GO:0035329; GO:0070062"
    ## "7" "MT_Pool_HQ_transcript/6" "sp|Q24292|DS_DROME" 6.93e-80 "Protein dachsous (Adherin)" "establishment of ommatidial planar polarity [GO:0042067]; establishment of planar polarity [GO:0001736]; homophilic cell adhesion via plasma membrane adhesion molecules [GO:0007156]; imaginal disc morphogenesis [GO:0007560]; negative regulation of hippo signaling [GO:0035331]; positive regulation of hippo signaling [GO:0035332]; protein localization involved in establishment of planar polarity [GO:0090251]" "GO:0001736; GO:0005509; GO:0005886; GO:0007156; GO:0007560; GO:0016327; GO:0035331; GO:0035332; GO:0042067; GO:0070161; GO:0090251"
    ## "8" "MT_Pool_HQ_transcript/7" "sp|E9Q555|RN213_MOUSE" 0 "E3 ubiquitin-protein ligase RNF213 (EC 2.3.2.27) (EC 3.6.4.-) (E3 ubiquitin-lipopolysaccharide ligase RNF213) (EC 2.3.2.-) (Mysterin) (RING finger protein 213)" "angiogenesis [GO:0001525]; defense response to bacterium [GO:0042742]; immune system process [GO:0002376]; lipid droplet formation [GO:0140042]; lipid ubiquitination [GO:0120323]; negative regulation of non-canonical Wnt signaling pathway [GO:2000051]; protein autoubiquitination [GO:0051865]; protein K63-linked ubiquitination [GO:0070534]; protein ubiquitination [GO:0016567]; regulation of lipid metabolic process [GO:0019216]; sprouting angiogenesis [GO:0002040]; ubiquitin-dependent protein catabolic process [GO:0006511]; xenophagy [GO:0098792]" "GO:0001525; GO:0002040; GO:0002376; GO:0004842; GO:0005524; GO:0005730; GO:0005737; GO:0005811; GO:0005829; GO:0006511; GO:0016567; GO:0016887; GO:0019216; GO:0042742; GO:0046872; GO:0051865; GO:0061630; GO:0070534; GO:0098792; GO:0120323; GO:0140042; GO:2000051"
    ## "9" "MT_Pool_HQ_transcript/8" "sp|Q8TD57|DYH3_HUMAN" 0 "Dynein axonemal heavy chain 3 (Axonemal beta dynein heavy chain 3) (HsADHC3) (Ciliary dynein heavy chain 3) (Dnahc3-b)" "cilium-dependent cell motility [GO:0060285]; microtubule-based movement [GO:0007018]" "GO:0003777; GO:0005524; GO:0005858; GO:0005874; GO:0007018; GO:0008569; GO:0030286; GO:0045505; GO:0051959; GO:0060285"
