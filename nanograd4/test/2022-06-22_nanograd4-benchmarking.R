#!/usr/bin/env Rscript

# Plot benchmarking data for nanograd4 runtime 
# written by AS on 2022-06-21

# data on onedrive; https://anu365-my.sharepoint.com/:t:/g/personal/u6081208_anu_edu_au/EX0rgbsN03lLs2XAraSljd4BnI7d_RwAg43IbVgS2imiWg?e=aVV08C

#################################################

library(tidyverse)
library(scales)
library(ggrepel)

run_time <- read_tsv("/Users/AJlocal/localGadiData/2022-06-22_nanograd4-benchmarking/2022-06-22_nanograd4-runtime-benchmarking.txt", col_names = T, skip = 1) 
run_time

cor.test(run_time$Reads, run_time$Run_time)

scientific_10 <- function(x) {
  gsub("e", " x 10^", scientific_format()(x))
}

dev.new(noRStudioGD = TRUE)

ggplot(run_time, aes(x = Reads, y = Run_time)) + 
  geom_point() +   
  theme_bw() +
  theme(text = element_text(size=16)) + 
  # ggtitle(str_wrap("Nanograd runtime on a laptop for human transcriptomes", 60)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  xlab("Number of reads in library") + 
  ylab("Runtime (seconds)") + 
  geom_smooth(method = "lm") + 
  scale_x_log10() + 
  coord_cartesian(ylim=c(0, 800))

