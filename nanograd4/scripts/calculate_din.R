#!/usr/bin/env Rscript

# written by AJ Sethi on 2022-06-23
# calculate DIN scores for all the HEK293 samples
# based on nanograd4 final output at /g/data/lf10/as7425/nanograd/analysis/2022-06-23_hek293-nanograd4-analysis/nanograd4_out


args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("\nUsage: Rscript calculate_din.R /path/to/nanograd4/summarised_transcripts.txt")
}

suppressPackageStartupMessages({
library(tidyverse)
})

######################################
######################################

# read files in
data <- read_tsv(args[1], col_names = T)

DIN <- data %>%
  mutate(df_ratio = df / tx_len) %>%
  summarise(sum = sum(df_ratio), n = n()) %>%
  mutate(DIN = 30 * sum/n) %>%
  select(DIN)

print(DIN)
