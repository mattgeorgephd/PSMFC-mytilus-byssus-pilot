## Molecular chaperones, their co-chaperones and the regulators of the two folding-stress
## responses, by family: the stress-protein genes of 09 (script 01's candidates and script 02's
## HSP_proteostasis module) since 2026-10-18 (Matt's choice, option d), in place of keywords
## on any part of the name ("heat shock", "Hsp", "BiP", "chaperone", ...), which had also caught
## histone, tubulin and ribosome-assembly chaperones and aliases of unrelated proteins.
##
## A gene's family comes from its best hit's recommended UniProt name (the text before the
## first parenthesis), so an alternative name cannot place it; the one exception is the J-domain
## proteins, recognised also by their HGNC alias "DnaJ homolog subfamily A/B/C member N", which
## UniProt keeps for J-proteins with other names (sacsin, GAK, SEC63, TIM14, HSCB, auxilin).
## Families and sources:
##   Hsp70/Hsp110, Hsp90, small HSPs, chaperonins (HSPD/HSPE, CCT): the human HSP families of
##     Kampinga et al. 2009, Cell Stress Chaperones 14:105 (HSPA, HSPH, HSPC, HSPB, HSPD/E, CCT);
##   J-domain proteins (DNAJ): Kampinga and Craig 2010, Nat Rev Mol Cell Biol 11:579;
##   Hsp100/Clp unfoldases: Schirmer et al. 1996, Trends Biochem Sci 21:289;
##   prefoldin: Vainberg et al. 1998, Cell 93:863;
##   Hsp70 co-chaperones (nucleotide-exchange factors, HIP, HOP, CHIP, SGT, the mtHsp70 escort
##     protein): Bracher and Verghese 2015, Front Mol Biosci 2:10, and Hartl et al. 2011,
##     Nature 475:324;
##   Hsp90 co-chaperones: Rohl et al. 2013, Trends Biochem Sci 38:253;
##   ER folding and quality control (calnexin, calreticulin, the PDI family, UGGT, the ER
##     PPIases, HSP47): Ellgaard and Helenius 2003, Nat Rev Mol Cell Biol 4:181;
##   regulators: heat-shock factor and its binding protein (Richter et al. 2010, Mol Cell
##     40:253) and the unfolded-protein-response sensors IRE1, PERK and ATF6 with XBP1 (Walter
##     and Ron 2011, Science 334:1081).
## Left out on purpose: chaperones of a single client or complex (tubulin-folding cofactors,
## MESD, NACHO, phosducin-like 3, histone chaperones, ribosomal-protein, proteasome and
## respiratory-complex assembly chaperones), proteases (Lon, ClpP, YME1L, HtrA2), cytosolic and
## spliceosomal PPIases, the PDI oxidases ERO1 and QSOX1, torsins, and proteins whose
## chaperone partnership is not their main role (FKBP8, TOM70, TSC1, FNIP1, R2TP).
CHAPERONE_FAMILIES <- c(
  "Hsp70/Hsp110"        = "^(heat shock (cognate )?70 kDa protein|heat shock cognate 71 kDa protein|heat shock protein 70\\b|heat shock protein (105|110) kDa|97 kDa heat shock protein|endoplasmic reticulum chaperone BiP|78 kDa glucose-regulated protein|stress-70 protein|hypoxia up-regulated protein 1|heat shock protein SS[A-Z][0-9]|chaperone protein DnaK)",
  "Hsp90"               = "^(heat shock protein HSP 90|heat shock protein (83|90)\\b|endoplasmin|heat shock protein 75 kDa, mitochondrial|ATP-dependent molecular chaperone HS[PC]82)",
  "small HSP"           = "^(heat shock protein beta-([1-9]|10)\\b|alpha-crystallin [AB] chain|small heat shock protein|heat shock protein (20|21|22|23|24|25|26|27|28|30)\\b|heat shock 27 kDa protein|heat shock protein hsp-16|protein lethal\\(2\\)essential for life)",
  "chaperonin"          = "^(60 kDa heat shock protein|10 kDa heat shock protein|T-complex protein 1 subunit|60 kDa chaperonin|10 kDa chaperonin)",
  "J-domain protein"    = "^(DnaJ homolog|chaperone protein DnaJ|DnaJ protein|DnaJ-like)",
  "Hsp100/Clp"          = "^(ATP-dependent clpX-like chaperone|caseinolytic peptidase B protein homolog|chaperone protein Clp[AB]\\b|heat shock protein 104)",
  "prefoldin"           = "^prefoldin subunit",
  "Hsp70 co-chaperone"  = "^(BAG family molecular chaperone regulator|Hsp70-binding protein|GrpE protein homolog|Hsc70-interacting protein|stress-induced-phosphoprotein|Hsp70-Hsp90 organi[sz]ing protein|E3 ubiquitin-protein ligase CHIP|DNL-type zinc finger protein|small glutamine-rich tetratricopeptide repeat-containing protein)",
  "Hsp90 co-chaperone"  = "^(Hsp90 co-chaperone Cdc37|prostaglandin E synthase 3|activator of 90 kDa heat shock protein ATPase|cysteine and histidine-rich domain-containing protein 1|protein SGT1 homolog|AH receptor-interacting protein|peptidyl-prolyl cis-trans isomerase (FKBP4|FKBP5|D)\\b|protein unc-45 homolog|Hsp70/Hsp90 co-chaperone CNS1|tetratricopeptide repeat protein 4\\b|nuclear migration protein nudC|serine/threonine-protein phosphatase 5\\b)",
  "ER folding"          = "^(calnexin|calreticulin|calmegin|protein disulfide-isomerase|endoplasmic reticulum resident protein (29|44)\\b|thioredoxin domain-containing protein 5\\b|UDP-glucose:glycoprotein glucosyltransferase|serpin H1\\b|peptidyl-prolyl cis-trans isomerase B\\b|peptidyl-prolyl cis-trans isomerase FKBP(2|7|9|10|11|14)\\b)",
  "stress regulator"    = "^(heat shock factor protein|heat shock transcription factor|heat shock factor-binding protein|serine/threonine-protein kinase/endoribonuclease IRE|eukaryotic translation initiation factor 2-alpha kinase 3|cyclic AMP-dependent transcription factor ATF-6|X-box-binding protein 1)")
J_PROTEIN_ALIAS <- "DnaJ homolog subfamily [ABC] member [0-9]+"

## the recommended name: the UniProt names string up to its first parenthesis or bracket
recommended_name <- function(protein_name) sub("\\s+[\\(\\[].*$", "", protein_name)

## each name's chaperone family, or NA; stops if a name falls in two families
chaperone_family <- function(protein_name) {
  rec <- recommended_name(protein_name)
  hit <- vapply(CHAPERONE_FAMILIES, function(rx) grepl(rx, rec, ignore.case = TRUE, perl = TRUE),
                logical(length(rec)))
  hit <- matrix(hit, nrow = length(rec))
  colnames(hit) <- names(CHAPERONE_FAMILIES)
  hit[, "J-domain protein"] <- hit[, "J-domain protein"] |
    grepl(J_PROTEIN_ALIAS, protein_name, ignore.case = TRUE, perl = TRUE)
  hit[is.na(protein_name), ] <- FALSE
  if (any(rowSums(hit) > 1))
    stop("names in two chaperone families: ", paste(protein_name[rowSums(hit) > 1], collapse = "; "))
  out <- rep(NA_character_, length(rec))
  out[rowSums(hit) == 1] <- colnames(hit)[max.col(hit[rowSums(hit) == 1, , drop = FALSE], ties.method = "first")]
  out
}
