## Helpers shared by the GO enrichment scripts (02-06). Sourced after _paths.R.
##
## The unit is the gene as it appears in the DESeq2 tables (e.g.
## "gene-LOC134682960|LOC134682960"); its annotation comes from the LOC key after the "|",
## through 01_go-inputs/gene_annotation.tsv written by script 01.

PADJ   <- 0.05                 # a DEG, as in 06_differential-expression
ONTS   <- c("BP", "MF", "CC")
MIN_GS <- 10                   # smallest GO term tested (annotated genes in the universe)
MAX_GS <- 500                  # largest, for clusterProfiler

contrasts <- read.csv(file.path(deg, "contrasts.csv"), colClasses = "character")

## Genes of one contrast: the universe is every gene DESeq2 gave an adjusted p (the genes it
## could have called; independent filtering leaves the rest at NA), split into up and down DEGs.
gene_sets <- function(code) {
  r <- contrasts[contrasts$code == code, ]
  stopifnot(nrow(r) == 1)
  d <- read.csv(file.path(deg, r$dir, paste0(code, "_apeglm.csv")))
  d <- d[!is.na(d$padj), ]
  list(universe = d$gene,
       up   = d$gene[d$padj < PADJ & d$log2FoldChange > 0],
       down = d$gene[d$padj < PADJ & d$log2FoldChange < 0])
}

## Gene annotation from script 01: one row per gene, GO IDs split by ontology.
read_annotation <- function() {
  a <- read.delim(file.path(out_inputs, "gene_annotation.tsv"), stringsAsFactors = FALSE,
                  colClasses = c(length = "numeric"))
  for (o in ONTS) a[[paste0("GO_", o)]][is.na(a[[paste0("GO_", o)]])] <- ""
  a
}

## gene -> direct GO IDs of one ontology, for the genes that have any (topGO's gene2GO)
gene2go_direct <- function(annot, ont, genes = annot$gene) {
  a <- annot[annot$gene %in% genes & nzchar(annot[[paste0("GO_", ont)]]), ]
  setNames(strsplit(a[[paste0("GO_", ont)]], ";", fixed = TRUE), a$gene)
}

## Long (GO, gene) table with every ancestor added, for goseq and clusterProfiler, which test
## the annotation as given (topGO walks the graph itself and takes the direct annotation).
gene2go_propagated <- function(annot, ont, genes = annot$gene) {
  direct <- gene2go_direct(annot, ont, genes)
  anc_map <- switch(ont, BP = GO.db::GOBPANCESTOR, MF = GO.db::GOMFANCESTOR, CC = GO.db::GOCCANCESTOR)
  ids <- unique(unlist(direct, use.names = FALSE))
  anc <- AnnotationDbi::mget(ids, anc_map, ifnotfound = NA)
  closure <- setNames(lapply(ids, function(i) unique(c(i, setdiff(anc[[i]], c("all", NA))))), ids)
  per_gene <- lapply(direct, function(g) unique(unlist(closure[g], use.names = FALSE)))
  data.frame(GO = unlist(per_gene, use.names = FALSE), gene = rep(names(per_gene), lengths(per_gene)))
}

go_term <- function(ids) {
  t <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, unique(ids), "TERM", "GOID"))
  unname(setNames(t$TERM, t$GOID)[ids])
}

## runs: one per contrast x direction, in contrasts.csv order
runs_for <- function(families) {
  cs <- contrasts[contrasts$family %in% families, ]
  data.frame(code = rep(cs$code, each = 2), family = rep(cs$family, each = 2),
             direction = rep(c("up", "down"), nrow(cs)))
}

ONT_NAME <- c(BP = "biological process", MF = "molecular function", CC = "cellular component")

## short labels for figures; the family (and so the reference group) goes in the title
contrast_label <- function(code) {
  tissue <- ifelse(substr(code, 1, 1) == "G", "Gill", "Foot")
  x <- substr(code, 2, 3)
  ifelse(code == "FG_TC", "Day-3 controls", paste(tissue, x))
}
FAMILY_TITLE <- c(TC = "stressor vs day-3 treatment control (TC, of record)",
                  FG = "gill vs foot in the day-3 controls (FG)")

## Dot plot of GO terms x contrasts, one panel per direction. `d` has one row per enriched term
## and run: code, direction, Term, p (the method's p for colour) and n (DEGs in term, or terms in
## a cluster). Colour is capped at P_CAP so one extreme term does not wash out the rest.
P_CAP <- 10
go_dotplot <- function(d, codes, title, caption, file, p_label = "-log10 p", n_label = "DEGs in term") {
  d <- d %>%
    mutate(direction = factor(ifelse(direction == "up", "Up-regulated", "Down-regulated"),
                              levels = c("Up-regulated", "Down-regulated")),
           contrast = factor(contrast_label(code), levels = unique(contrast_label(codes))),
           Term = factor(Term, levels = rev(unique(Term[order(match(code, codes), p)]))),
           logp = pmin(-log10(p), P_CAP))
  p <- ggplot(d, aes(contrast, Term)) +
    geom_point(aes(size = n, colour = logp)) +
    facet_wrap(~ direction, drop = FALSE) +
    scale_colour_gradient(low = "#cfcfcf", high = "#1a1a1a", limits = c(0, P_CAP),   # neutral: red and blue mean up and down
                          breaks = seq(0, P_CAP, 2.5), labels = c(seq(0, P_CAP - 2.5, 2.5), paste0(P_CAP, "+")),
                          name = p_label) +
    scale_size_area(max_size = 6, name = n_label) +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = title,
         caption = paste(strwrap(caption, width = 110), collapse = "\n")) +
    theme_psmfc(10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.caption = element_text(hjust = 0),
          plot.caption.position = "plot", plot.title.position = "plot")
  n_codes <- length(unique(codes))
  ggsave(file, p, width = 4.6 + 0.55 * 2 * n_codes, height = 2.4 + 0.2 * nlevels(d$Term),
         dpi = 300, limitsize = FALSE)
  p
}
