#!/usr/bin/env Rscript

# analyse merged sequencing summary and degradation data 
# written by AS on 2022-06-21

#################################################

library(tidyverse)

data_path <- "~/localGadiData/2022-06-21_degradation-all-data/2022-06-21_all_degradation_combined.txt.gz"

# import data 
input <- read_tsv(data_path, col_names = T, col_types = "ffddfddfffddf") %>% 
  mutate(condition = factor(condition, levels = c("wt_rep1", "wt_rep2", "deg_rep1", "deg_rep2")))
  
# save data as unzipped text file 
# write_tsv(input, "/Users/AJlocal/localGadiData/2022-06-21_degradation-all-data/2022-06-21_all_degradation_combined.txt", col_names = T)

#################################################

# make a histogram of read duration 
ggplot(input %>% slice_sample(n = 10000), aes(x = duration, color = condition)) + 
  geom_histogram() + 
  scale_y_log10() + 
  xlim(5,20) + 
  facet_wrap(~condition, nrow = 2)

# check the median duration using summarise 
cond_end <- input %>% 
  slice_sample(n = 10000) %>% 
  group_by(condition, end_reason) %>% 
  summarise(median_duration = median(duration), median_seq_length = median(seq_len), median_align_length = median(map_len)) %>% 
  mutate(sequence_rate = median_seq_length / median_duration, align_rate = median_align_length / median_duration) %>% 
  arrange(desc(align_rate))
  unite(group, c(condition, end_reason), sep="_") %>% 
  arrange(desc(median_duration))

# plot the median duration of each condition_end-reason
ggplot(cond_end, aes(x = group, y = median_duration)) + geom_point()

# compare start time to sequencing duration between the 4 conditions 
dev.new(noRStudioGD = TRUE)
ggplot(input %>% slice_sample(n = 100000), aes(x = start_time, y = duration, color = condition)) + 
  geom_point() + 
  facet_wrap(~condition, nrow = 2) 
#################################################

# use group_by and summarise() to look at transcript and condition-specific enrichment of read_end reason 

# summarise end_reason data 
end_reason_data <- input %>% 
  group_by(condition, end_reason) %>% 
  summarise(n = n())
end_reason_data

# count the read end reasons for every transcript in every condition
tx_cond_reason <- input %>%
  # slice_sample(n = 25000) %>% 
  group_by(condition, transcript_id) %>% 
  # mutate(n_reads = n()) %>% # when mutating grouped data, one mutate is performed per_group and no data is lost 
  # filter(n_reads > 20) %>%  # when desigining a heuristic like this, it's worth looking at the data distribution and testing different cutoffs \
  # select(-n_reads) %>% 
  group_by(condition, transcript_id, end_reason) %>% 
  summarise(count = n()) %>% 
  pivot_wider(names_from = end_reason, values_from = count) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  mutate(tx_total_count = sum(signal_positive + unblock_mux_change + signal_negative))

# get the ratio of unblock to signal_positive for every transcript in tx_cond_reason
unblock_tibble <- tx_cond_reason %>% 
  mutate(unblock_to_positive = unblock_mux_change / signal_positive)

# select wt_rep1 
wt_rep1_unblock <- unblock_tibble %>% filter(condition == "wt_rep1")

# 
ggplot(wt_rep1_unblock, aes(x = tx_total_count, y = unblock_to_positive)) + 
  geom_point()

##################################################
##################################################

# export a plot of unblock to signal positive termination events 

# distribution of unblock to positive ratio 
ggplot(unblock_tibble %>% filter(tx_total_count > 20), aes(x = unblock_to_positive, color = condition)) + 
  scale_x_log10() +
  stat_ecdf() +   
  theme_bw() +
  theme(text = element_text(size=16)) + 
  ggtitle(str_wrap("Ratio of unblock termination to signal positive termination per transcript", 60)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  xlab("Ratio of unblock termination to signal positive termination") + 
  ylab("Transcript density")

# save the plot as both a vector and as raster
ggsave("/Users/AJlocal/Documents/nanograd/plots/2022-06-21_unblock-ratio-distribution_degradation.png", plot = last_plot(), scale = 1, width = 8.5, height = 4.5, unit = "in", dpi = 4500, bg = NULL)

##################################################
##################################################

ggplot(aes(x = value, fill = stat)) + 
  geom_histogram(bins = 100, alpha = 1, position = "identity") + 
  facet_wrap(~stat, nrow = 2) +




# plot transcript count vs unblock to positive ratio 
ggplot(unblock_tibble, aes(x = tx_total_count, y = unblock_to_positive)) + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_hex()



# quickly, plot the distribution of counts per transcript within a condition using our temp data 
ggplot(tx_cond_reason, aes(x = tx_total_count, fill = condition)) + 
  geom_histogram() + 
  facet_wrap(~condition, nrow = 2) + 
  scale_x_log10() + 
  scale_y_log10()

tx_cond_reason

#################################################

# explore HMGN2
hmg_data <- input %>% filter(condition == "wt_rep1") %>% filter(transcript_id == "ENST00000361427") %>% 
  mutate(translocation_rate = duration / map_len)

ggplot(hmg_data, aes(x = seq_len, y = map_len, color = translocation_rate)) + 
  geom_point()

ggplot(input %>% filter(transcript_id == "ENST00000361427"), aes(x = map_len, color = condition)) + stat_ecdf()

# calculate transcript integrity factors 
input %>% filter(transcript_id == "ENST00000361427") %>% mutate(df = map_len/(1940-670)) %>% group_by(condition) %>% summarise(df = median(df))


ggplot(input %>% filter(transcript_id == "ENST00000361427") %>% mutate(df = map_len/(1940-670)), aes(x = df, color = condition)) + 
  stat_ecdf() + 
  theme_bw() +
  theme(text = element_text(size=16)) + 
  ggtitle(str_wrap("HMGN2 per-read integrity numbers", 60)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  xlab("Read-specific integrity numbers") + 
  ylab("Fraction of data") + 
  coord_cartesian(ylim = c(0,1), xlim=c(0,1))

# make a plot of end_reason 
dev.new(noRStudioGD = TRUE)

# plot to a device 
ggplot(input %>% filter(end_reason != "mux_change"), aes(x = df, color = end_reason)) + 
  stat_ecdf() + 
  theme_bw() +
  theme(text = element_text(size=16)) + 
  ggtitle(str_wrap("Overall read integrity numbers", 60)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  xlab("Read-specific integrity numbers") + 
  ylab("Fraction of data") + 
  coord_cartesian(ylim = c(0,1), xlim = c(0,1)) + 
  facet_wrap(~condition, nrow = 2)

#################################################

# plot transcript integrity factor for all libraries
input %>% group_by(condition, transcript_id) %>% 
  summarise(df = median(df), count = n()) %>% 
  ggplot(., aes(x = df, color = condition)) +
  stat_ecdf() + 
  theme_bw() +
  theme(text = element_text(size=16)) + 
  ggtitle(str_wrap("Transcriptome-wide integrity factors", 60)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  xlab("Transcript-specific integrity factors") + 
  ylab("Fraction of transcripts") + 
  coord_cartesian(ylim = c(0,1), xlim = c(0,1))




