# Title: Material Testing Data Extractor
# Author: Monica Klopp
# Date: 05/2022

# Set working directory
home <- 'C:/Users/Carrington_Lab/Desktop/PSMFC-mytilus-byssus-pilot'
setwd(home); getwd()

# Install packages (run first time only)
# install.packages("tidyverse")

# Load packages
require(ggplot2)

# Imports the data
data <- read.delim('C:/Users/Carrington_Lab/Desktop/PSMFC-mytilus-byssus-pilot/thread_strength/tensometer_output/T01_01.txt', # location of file  
                   header = FALSE, # is there a header? - NaNs broke it, turn it off
                   sep = "\t", # how is the data separated?
                   dec = ".")

# Flip data from row to column format
data_transpose = t(data) # transpose

# Fix headers
data_snipped <- data_transpose[,2:4] # keep all rows, grab only column 2 through 4
colnames(data_snipped) <- c('Time', 'Displacement', 'Force')

# Convert to dataframe
data_snipped <- as.data.frame(data_snipped) # change character vector to dataframe
data_snipped <- lapply(data_snipped, function(x) as.numeric(data_snipped)) # change characters to numeric [this doesnt work yet - maybe we should avoid this whole issue by not importing text with delim]

# Plot data using ggplot 
p1  <- ggplot(data_snipped, # what data set to plot?
              aes(x=Time, y=Force)) + # aes stands for aesthetic mapping
       geom_line() +
       geom_point()

p1
              
              
              
#               +
#        geom_boxplot(colour = "grey30", size = 0.8,outlier.colour="grey30", outlier.shape = 16,
#                      outlier.size=1, notch=FALSE) +
#        scale_fill_manual(values=trt_list$trt_colors[5:24]) +
#        # geom_point() +
#        scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.52)) +
#        # scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 32)) + 
#        my_theme
# 
# 
# p1


# Find maximum force

# Find maximum slope (extra credit)
