# Title: generate plots PSMFC-mytilus-byssus-pilot
# Author: Matthew George; mattgeorgephd@gmail.com
# Date: 11/2021

## clear
rm(list=ls())

## Grab the WD from the file location
library(rstudioapi)
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path )); getwd()

## Load R packages
library(readxl)
library(ggplot2)
library(stringr)
library(tidyverse)
library(Johnson)
library(agricolae)
library(nlme)
library(multcomp)

## Set ggplot theme
my_theme <- theme(line              = element_line(size=1.5),
                  rect              = element_rect(size=1.5),
                  text              = element_text(size=14,color="black"),
                  panel.background  = element_blank(),
                  panel.grid.major  = element_blank(), 
                  panel.grid.minor  = element_blank(),
                  axis.text.x       = element_text(size=16,color="black"),
                  axis.text.y       = element_text(size=16,color="black"),
                  axis.title.x      = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                  axis.title.y      = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.ticks.x      = element_line(color="black"),
                  axis.ticks.y      = element_line(color="black"),
                  # axis.line         = element_line(color = "black", size = 0.1),
                  panel.border      = element_rect(color = "black", fill=NA, size=1.5),
                  legend.key        = element_blank()) # removes background of legend bullets



###########################################################################################################################
### [1] Boxplot - SMR vs. time by treatment - OA, OW, DO

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path )); setwd('respirometry/output'); getwd()

trt_list <- read_excel("processed_summary.xlsx", sheet = "trt_list", col_names = TRUE)

trt_plot            <- read_excel("processed_summary.xlsx", sheet = "SMR", col_names = TRUE)
trt_plot$timepoint  <- factor(trt_plot$timepoint, levels=c("baseline","final"),ordered=TRUE)
trt_plot$trt_list   <- factor(trt_plot$trt_list,levels=trt_list$trt_list,ordered=TRUE)
trt_plot$trt2_list   <- factor(trt_plot$trt2_list,levels=c("gallo-baseline","tross-baseline","gallo-final","tross-final"),ordered=TRUE)

trt_plot_OA <- trt_plot %>% filter(trt == "OA")
trt_plot_OW <- trt_plot %>% filter(trt == "OW")
trt_plot_DO <- trt_plot %>% filter(trt == "DO")

p1 <- ggplot(trt_plot_OA, aes(x=timepoint, y=SMR, group=as.factor(ID), color=species)) +
      # geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
      #              outlier.size=1, notch=FALSE) +
      # scale_fill_manual(values=c("orangered1","royalblue1")) +
      geom_point() +
      geom_line(size=1) +
      scale_color_manual(values=c("orangered1","royalblue1")) +
      # facet_wrap(~ID) +
      scale_y_continuous(breaks = seq(0, 0.3, 0.1), limits = c(0, 0.32)) +
      # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
      # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
      my_theme

p1

bp1 <- ggplot(trt_plot_OA, aes(x=timepoint, y=SMR, group=as.factor(trt2_list),fill=species)) +
      geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
                   outlier.size=1, notch=FALSE) +
      scale_fill_manual(values=c("orangered1","royalblue1")) +
      # geom_point() +
      # geom_line(size=1) +
      # scale_color_manual(values=c("orangered1","royalblue1")) +
      # facet_wrap(~ID) +
      scale_y_continuous(breaks = seq(0, 0.3, 0.1), limits = c(0, 0.32)) +
      # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
      # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
      my_theme

bp1

p2 <- ggplot(trt_plot_OW, aes(x=timepoint, y=SMR, group=as.factor(ID), color=species)) +
      # geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
      #              outlier.size=1, notch=FALSE) +
      # scale_fill_manual(values=c("orangered1","royalblue1")) +
      geom_point() +
      geom_line(size=1) +
      scale_color_manual(values=c("orangered1","royalblue1")) +
      # facet_wrap(~ID) +
      scale_y_continuous(breaks = seq(0, 0.3, 0.1), limits = c(0, 0.32)) +
      # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
      # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
      my_theme

p2

bp2 <- ggplot(trt_plot_OW, aes(x=timepoint, y=SMR, group=as.factor(trt2_list),fill=species)) +
        geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
                     outlier.size=1, notch=FALSE) +
        scale_fill_manual(values=c("orangered1","royalblue1")) +
        # geom_point() +
        # geom_line(size=1) +
        # scale_color_manual(values=c("orangered1","royalblue1")) +
        # facet_wrap(~ID) +
        scale_y_continuous(breaks = seq(0, 0.3, 0.1), limits = c(0, 0.32)) +
        # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
        # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
        my_theme

bp2

p3 <- ggplot(trt_plot_DO, aes(x=timepoint, y=SMR, group=as.factor(ID), color=species)) +
      # geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
      #              outlier.size=1, notch=FALSE) +
      # scale_fill_manual(values=c("orangered1","royalblue1")) +
      geom_point() +
      geom_line(size=1) +
      scale_color_manual(values=c("orangered1","royalblue1")) +
      # facet_wrap(~ID) +
      scale_y_continuous(breaks = seq(0, 0.6, 0.1), limits = c(0, 0.62)) +
      # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
      # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
      my_theme

p3

bp3 <- ggplot(trt_plot_DO, aes(x=timepoint, y=SMR, group=as.factor(trt2_list),fill=species)) +
      geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
                   outlier.size=1, notch=FALSE) +
      scale_fill_manual(values=c("orangered1","royalblue1")) +
      # geom_point() +
      # geom_line(size=1) +
      # scale_color_manual(values=c("orangered1","royalblue1")) +
      # facet_wrap(~ID) +
      scale_y_continuous(breaks = seq(0, 0.6, 0.1), limits = c(0, 0.62)) +
      # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
      # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
      my_theme

bp3

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path )); setwd('plots'); getwd()

ggsave("BOXPLOT_SMR_OA.png",
       plot   = bp1,
       dpi    = 600,
       device = "png",
       width  = 4,
       height = 4,
       units  = "in")

ggsave("BOXPLOT_SMR_OW.png",
       plot   = bp2,
       dpi    = 600,
       device = "png",
       width  = 4,
       height = 4,
       units  = "in")

ggsave("BOXPLOT_SMR_DO.png",
       plot   = bp3,
       dpi    = 600,
       device = "png",
       width  = 4,
       height = 4,
       units  = "in")

ggsave("LINEGRAPH_SMR_OA.png",
       plot   = p1,
       dpi    = 600,
       device = "png",
       width  = 4,
       height = 4,
       units  = "in")

ggsave("LINEGRAPH_SMR_OW.png",
       plot   = p2,
       dpi    = 600,
       device = "png",
       width  = 4,
       height = 4,
       units  = "in")

ggsave("LINEGRAPH_SMR_DO.png",
       plot   = p3,
       dpi    = 600,
       device = "png",
       width  = 4,
       height = 4,
       units  = "in")

 
###########################################################################################################################
### [2] Boxplot - SMR vs. time by treatment - OA, OW, DO

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path )); setwd('morphometrics'); getwd()

trt_list <- read_excel("morphometrics.xlsx", sheet = "trt_list", col_names = TRUE)

trt_plot              <- read_excel("morphometrics.xlsx", sheet = "thread_count", col_names = TRUE)
trt_plot$species      <- factor(trt_plot$species)
trt_plot$trt          <- factor(trt_plot$trt)
trt_plot$time         <- factor(trt_plot$time,levels=c("before","after"),ordered=TRUE)
trt_plot$combo         <- factor(trt_plot$combo)
# trt_plot$trt2_list   <- factor(trt_plot$trt2_list,levels=c("gallo-baseline","tross-baseline","gallo-final","tross-final"),ordered=TRUE)

trt_plot_CC <- trt_plot %>% filter(trt == "control")
trt_plot_OA <- trt_plot %>% filter(trt == "OA")
trt_plot_OW <- trt_plot %>% filter(trt == "OW")
trt_plot_DO <- trt_plot %>% filter(trt == "DO")

p1 <- ggplot(trt_plot_OA, aes(x=time, y=thread, group=as.factor(ID), color=species)) +
  # geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
  #              outlier.size=1, notch=FALSE) +
  # scale_fill_manual(values=c("orangered1","royalblue1")) +
  geom_point() +
  geom_line(size=1) +
  scale_color_manual(values=c("orangered1","royalblue1")) +
  # facet_wrap(~ID) +
  scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
  # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
  # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
  my_theme

p1

bp1 <- ggplot(trt_plot_OA, aes(x=time, y=thread,fill=species)) +
  geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
               outlier.size=1, notch=FALSE) +
  scale_fill_manual(values=c("orangered1","royalblue1")) +
  # geom_point() +
  # geom_line(size=1) +
  # scale_color_manual(values=c("orangered1","royalblue1")) +
  # facet_wrap(~ID) +
  scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
  # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
  # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
  my_theme

bp1

p2 <- ggplot(trt_plot_OW, aes(x=time, y=thread, group=as.factor(ID), color=species)) +
  # geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
  #              outlier.size=1, notch=FALSE) +
  # scale_fill_manual(values=c("orangered1","royalblue1")) +
  geom_point() +
  geom_line(size=1) +
  scale_color_manual(values=c("orangered1","royalblue1")) +
  # facet_wrap(~ID) +
  scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
  # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
  # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
  my_theme

p2

bp2 <- ggplot(trt_plot_OW, aes(x=time, y=thread,fill=species)) +
  geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
               outlier.size=1, notch=FALSE) +
  scale_fill_manual(values=c("orangered1","royalblue1")) +
  # geom_point() +
  # geom_line(size=1) +
  # scale_color_manual(values=c("orangered1","royalblue1")) +
  # facet_wrap(~ID) +
  scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
  # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
  # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
  my_theme

bp2

p3 <- ggplot(trt_plot_DO, aes(x=time, y=thread, group=as.factor(ID), color=species)) +
  # geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
  #              outlier.size=1, notch=FALSE) +
  # scale_fill_manual(values=c("orangered1","royalblue1")) +
  geom_point() +
  geom_line(size=1) +
  scale_color_manual(values=c("royalblue1")) +
  # facet_wrap(~ID) +
  scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
  # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
  # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
  my_theme

p3

bp3 <- ggplot(trt_plot_DO, aes(x=time, y=thread,fill=species)) +
  geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
               outlier.size=1, notch=FALSE) +
  scale_fill_manual(values=c("royalblue1")) +
  # geom_point() +
  # geom_line(size=1) +
  # scale_color_manual(values=c("orangered1","royalblue1")) +
  # facet_wrap(~ID) +
  scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
  # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
  # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
  my_theme

bp3

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path )); setwd('plots'); getwd()

ggsave("BOXPLOT_thread_OA.png",
       plot   = bp1,
       dpi    = 600,
       device = "png",
       width  = 4,
       height = 4,
       units  = "in")

ggsave("BOXPLOT_thread_OW.png",
       plot   = bp2,
       dpi    = 600,
       device = "png",
       width  = 4,
       height = 4,
       units  = "in")

ggsave("BOXPLOT_thread_DO.png",
       plot   = bp3,
       dpi    = 600,
       device = "png",
       width  = 4,
       height = 4,
       units  = "in")

ggsave("LINEGRAPH_thread_OA.png",
       plot   = p1,
       dpi    = 600,
       device = "png",
       width  = 4,
       height = 4,
       units  = "in")

ggsave("LINEGRAPH_thread_OW.png",
       plot   = p2,
       dpi    = 600,
       device = "png",
       width  = 4,
       height = 4,
       units  = "in")

ggsave("LINEGRAPH_thread_DO.png",
       plot   = p3,
       dpi    = 600,
       device = "png",
       width  = 4,
       height = 4,
       units  = "in")


###########################################################################################################################
### [1] Boxplot - Condition vs. time by treatment - OA, OW, DO

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path )); setwd('morphometrics'); getwd()

# trt_list <- read_excel("processed_summary.xlsx", sheet = "trt_list", col_names = TRUE)

trt_plot            <- read_excel("morphometrics.xlsx", sheet = "condition", col_names = TRUE)
# trt_plot$timepoint  <- factor(trt_plot$timepoint, levels=c("baseline","final"),ordered=TRUE)
trt_plot$trt         <- factor(trt_plot$trt,levels=c('baseline','final','OA','OW','DO'),ordered=TRUE)
trt_plot$species     <- factor(trt_plot$species,levels=c("gallo","tross"),ordered=TRUE)

# trt_plot_OA <- trt_plot %>% filter(trt == "OA")
# trt_plot_OW <- trt_plot %>% filter(trt == "OW")
# trt_plot_DO <- trt_plot %>% filter(trt == "DO")

bp4 <- ggplot(trt_plot, aes(x=trt, y=CI,group=combo,fill=species)) +
        geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
                     outlier.size=1, notch=FALSE) +
        scale_fill_manual(values=c("orangered1","royalblue1")) +
        # facet_wrap(~ID) +
        # scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
        # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
        # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
        my_theme

bp4

bp5 <- ggplot(trt_plot, aes(x=trt, y=GI,group=combo,fill=species)) +
        geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
                     outlier.size=1, notch=FALSE) +
        scale_fill_manual(values=c("orangered1","royalblue1")) +
        # facet_wrap(~ID) +
        # scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
        # scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0, 22)) +
        # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) +
        my_theme

bp5

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path )); setwd('plots'); getwd()

ggsave("BOXPLOT_CI.png",
       plot   = bp4,
       dpi    = 600,
       device = "png",
       width  = 6,
       height = 4,
       units  = "in")

ggsave("BOXPLOT_GI.png",
       plot   = bp5,
       dpi    = 600,
       device = "png",
       width  = 6,
       height = 4,
       units  = "in")

# ### STATS
# 
# library(tidyverse)
# library(ggpubr)
# library(rstatix)
# 
# 
# 
# x = MR_plot$SMR
# qqnorm(x) # check for normality
# qqline(x) # Draw the line
# result <- shapiro.test(x) # p-value fail = good, don't need transformation
# print(result$p.value)
# if(result$p.value<0.05)     {
#   x_johnson <- RE.Johnson(x) # transform
#   x_transformed = x_johnson$transformed
#   qqnorm(x_transformed) # check linearity of tranformed data
#   qqline(x_transformed)
#   print(shapiro.test(x_transformed))
#   x <- x_transformed  
#   print("transformed!",quote=FALSE)}
# shapiro.test(x)
# 
# MR_plot$SMR <- x
# 
# 
# summary(my_anova <- aov(SMR ~ ploidy*timepoint*|, data = MR_plot))
# TukeyHSD(my_anova, ordered = TRUE)
# 
# 
# res.aov <- anova_test(
#   data = MR_plot, dv = SMR, wid = as.factor(MR_plot$ID),
#   within = c(timepoint, ploidy)
# )
# get_anova_table(res.aov)

######################################################################################################
## STATS - heated treatment\


library(tidyverse)
library(ggpubr)
library(rstatix)

trt_plot_OA <- trt_plot_OA %>% filter(species == "gallo")             

res.aov <- anova_test(data = trt_plot_DO, dv = SMR, wid = ID, within = timepoint)
get_anova_table(res.aov)

## Transform data
x = trt_plot_OA$thread
qqnorm(x) # check for normality
qqline(x) # Draw the line
result <- shapiro.test(x) # p-value fail = good, don't need transformation
print(result$p.value)
if(result$p.value<0.05)     {
  x_johnson <- RE.Johnson(x) # transform
  x_transformed = x_johnson$transformed
  qqnorm(x_transformed) # check linearity of tranformed data
  qqline(x_transformed)
  print(shapiro.test(x_transformed))
  x <- x_transformed
  print("transformed!",quote=FALSE)}
shapiro.test(x)

trt_plot_OA$SMR <- x

# trt_plot$timepoint <- as.factor(trt_plot$timepoint)
trt_plot$trt <- as.factor(trt_plot$trt)
trt_plot$ID <- as.factor(trt_plot$ID)

# Statistical testings
summary(aov(thread ~ time, data = trt_plot_DO))

tx <- with(trt_plot, interaction(ploidy,timepoint))
amod <- aov(SMR ~ tx, data=MR_plot)
HSD.test(amod, "tx", group=TRUE, console=TRUE)

# FR
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path )); setwd('feeding_rate'); getwd()

trt_list <- read_excel("FR.xlsx", sheet = "trt_list", col_names = TRUE)

FR_plot           <- read_excel("FR.xlsx", sheet = "heat", col_names = TRUE)
# FR_plot           <- filter(FR_plot, FR_plot$death == "yes")
FR_plot$ploidy    <- factor(FR_plot$ploidy, levels=c("D","T"),ordered=TRUE)
# MR_plot$temp      <- factor(MR_plot$temp, levels=c("10","30"),ordered=TRUE)
FR_plot$timepoint <- factor(FR_plot$timepoint, levels=c("-10","1","2","6","10"),ordered=TRUE)
FR_plot$trt       <- factor(FR_plot$trt,levels=trt_list$trt_list,ordered=TRUE)

## Transform data
x = FR_plot$size_corrected_FR
qqnorm(x) # check for normality
qqline(x) # Draw the line
result <- shapiro.test(x) # p-value fail = good, don't need transformation
print(result$p.value)
if(result$p.value<0.05)     {
  x_johnson <- RE.Johnson(x) # transform
  x_transformed = x_johnson$transformed
  qqnorm(x_transformed) # check linearity of tranformed data
  qqline(x_transformed)
  print(shapiro.test(x_transformed))
  x <- x_transformed
  print("transformed!",quote=FALSE)}
shapiro.test(x)

FR_plot$size_corrected_FR <- x

# Statistical testings
summary(aov(size_corrected_FR ~ ploidy*timepoint, data = FR_plot))

tx <- with(FR_plot, interaction(ploidy,timepoint))
amod <- aov(size_corrected_FR ~ tx, data= FR_plot)
HSD.test(amod, "tx", group=TRUE, console=TRUE)


# summary(my_anova <- aov(mortality ~ temp*p, data = MR_plot))
# TukeyHSD(my_anova, ordered = TRUE)
# 
# summary(my_anova <- aov(umol_L_hr ~ timepoint*temp*ploidy, data = MR_plot))
# TukeyHSD(my_anova, ordered = TRUE)
