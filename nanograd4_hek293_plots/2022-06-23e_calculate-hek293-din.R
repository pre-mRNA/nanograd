#!/usr/bin/env Rscript

# written by AJ Sethi on 2022-06-23
# calculate DIN scores for all the HEK293 samples
# based on nanograd4 final output at /g/data/lf10/as7425/nanograd/analysis/2022-06-23_hek293-nanograd4-analysis/nanograd4_out

suppressPackageStartupMessages({
library(tidyverse)
})

######################################
######################################

# set the working directoryw
setwd("/g/data/lf10/as7425/nanograd/analysis/2022-06-23_hek293-nanograd4-analysis/nanograd4_out")

# list all files
files <- dir(pattern = "*.txt")

# omit the temporary files
files <- files[seq(1, 28, 2)]

# read files in
data <- files %>%
  map(read_tsv() %>%    # read in all the files individually, using
  mutate(df_ratio = df / tx_len) %>%
  summarise(sum = sum(df_ratio), n = n()) %>%
  mutate(DIN = 30 * sum/n) %>%
  select(DIN))
data
