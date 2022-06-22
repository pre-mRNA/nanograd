#!/usr/bin/env Rscript

# analyse merged sequencing summary and degradation data 
# written by AS on 2022-06-22
# test the updated version of nanograd4 which relies on transcriptome rather than genome annotations 

#################################################

library(tidyverse)

data_path <- "~/localGadiData/2022-06-22_nanograd4-benchmarking/all_degradation_combined.txt.gz"

# import data 
input <- read_tsv(data_path, col_names = T, col_types = "ffddfddfffddf") %>% 
  mutate(condition = factor(condition, levels = c("wt_rep1", "wt_rep2", "deg_rep1", "deg_rep2")))

#################################################

# start by plotting basic properties of the data 
wt1 <- input %>% filter(condition == "wt_rep1")

# remove input 
rm(input)

# make a data summary that we can check on to look at indiviudal transcripts 
ctrl_summary <- wt1 %>% 
  select(-condition) %>% 
  group_by(gene_name, transcript_id, transcript_biotype, tx_len) %>% 
  summarise(reads = n(), df_var = var(df), df_mean = mean(df), df_median = median(df))
ctrl_summary

# plot the mean-variance relationship 
ggplot(ctrl_summary %>% filter(reads > 100), aes(x = df_mean, y = df_var)) + geom_point() + scale_x_log10() + scale_y_log10()

# plot the tendency of the end_reasons over time 
ggplot(wt1 %>% filter(), aes(x = start_time, y = map_len)) + 
  geom_hex()

# plot the mapped_length vs duration 
ggplot(wt1, aes(x = duration, y = map_len, color = end_reason)) + 
    geom_point() + 
    facet_wrap(~end_reason, nrow = 2)

# plot the sequence_length vs map_length for the 4 end_reasons 
ggplot(wt1, aes(x = start_time, y = map_len)) + 
  geom_hex() + 
  facet_wrap(~end_reason, nrow = 2)

# histogram
ggplot(wt1, aes(x = map_len)) + 
  geom_histogram() + 
  scale_y_log10() + 
  facet_wrap(~end_reason, nrow = 2)

#################################################

# pick a gene, GAPDH, and plot data for that gene 
ctrl_gap <- wt1 %>% filter(transcript_id == "ENST00000229239")

# plot the sequence_length vs map_length for the 4 end_reasons 
ggplot(ctrl_gap, aes(x = seq_len, y = map_len, color = end_reason)) + geom_point()

# plot the sequence_length vs map_length for the 4 end_reasons 
ggplot(ctrl_gap, aes(x = duration, y = map_len, color = end_reason)) + geom_point()




#################################################

# pick a gene, GAPDH, and plot data for that gene 
ctrl_r2 <- wt1 %>% filter(transcript_id == "ENST00000343262")

# plot the sequence_length vs map_length for the 4 end_reasons 
ggplot(ctrl_r2, aes(x = seq_len, y = map_len, color = end_reason)) + geom_point()

# plot the sequence_length vs map_length for the 4 end_reasons 
ggplot(ctrl_r2, aes(x = duration, y = map_len, color = end_reason)) + geom_point()

# plot the sequence_length vs map_length for the 4 end_reasons 
ggplot(ctrl_r2, aes(x = start_time, y = map_len, color = end_reason)) + 
  geom_point()

#################################################

# END 

# scratch area: 
a <- head(input)
