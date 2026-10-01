## Shared colours and theme for every figure in the repository (02, 06, 07, 08, 09).
## Base R + ggplot2. Scripts source this file instead of defining their own colours.
##
## Treatments keep the hues of the thread-strength figures: control, the reference group, in a
## neutral grey; OA green, OW orange, DO purple. The original values failed a colour-vision
## check (DO against control, and baseline against treatment control, were hard to tell apart
## even with full colour vision), so the same hues sit at steps that pass: worst all-pairs
## deltaE 9.2 (deutan) and 16.7 (normal vision), thread-condition extras included. OA's green
## is below 3:1 contrast on white, so every figure keeps a legend and, where points overlap,
## a redundant shape per treatment.
TREATMENT_COLORS <- c(control = "#6e6e6e", OA = "#1baf7a", OW = "#eb6834", DO = "#7b3294")
TREATMENT_SHAPES <- c(control = 16, OA = 15, OW = 17, DO = 18)
TREATMENT_LABELS <- c(control = "Control", OA = "Ocean acidification (OA)",
                      OW = "Ocean warming (OW)", DO = "Hypoxia (DO)")

## Thread condition (02_thread-strength): what a thread was built in. The arms above, plus
## the pre-exposure baseline (blue) and the day-0 lab reference, outside the experiment (light grey).
THREAD_TRT_COLORS <- c(lab_control = "#bdbdbd", baseline = "#2a78d6",
                       treatment_control = TREATMENT_COLORS[["control"]],
                       OW = TREATMENT_COLORS[["OW"]], OA = TREATMENT_COLORS[["OA"]],
                       DO = TREATMENT_COLORS[["DO"]])

## Tissue: hues kept apart from the treatment hues. Foot is the phenol gland to the tip of the
## foot, the region sampled in every animal. The twelve day-0 animals also have a library of
## the rest of the foot (IDs ending FX); REGION_* colours the three sampled regions in figures
## that include those libraries (all-pairs deltaE 17.6 CVD, 33.9 normal vision; the pink is
## below 3:1 contrast on white, so those figures keep a legend).
TISSUE_COLORS <- c(F = "#4a3aa7", G = "#008300")
TISSUE_LABELS <- c(F = "Foot", G = "Gill")
REGION_COLORS <- c("foot (phenol gland to tip)" = "#4a3aa7", "foot (without phenol gland)" = "#e87ba4",
                   gill = "#008300")
REGION_LABELS <- c("foot (phenol gland to tip)" = "Foot: phenol gland to tip",
                   "foot (without phenol gland)" = "Foot without phenol gland (day 0 only)",
                   gill = "Gill")

## Direction of change, in every figure: red up, blue down, a neutral grey for "not
## significant". Tissue and region colours above are never combined with treatment colours
## in one figure.
DIRECTION_COLORS <- c(up = "#e34948", down = "#2a78d6", ns = "#bdbdbd")

theme_psmfc <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70"),
                   legend.key = ggplot2::element_blank())
}
