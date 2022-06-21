#!/usr/bin/env Rscript

# Plot unannotated and annotated signal data for the SBRS talk 
# written by AS on 2022-06-21

#################################################

library(tidyverse)

# read in nanopolish eventalign file (path for JCSMR iMac)

a <- read_tsv("~/localGadiData/2022-06-17_signal-plots-SXL/26017490-ad90-47a2-a66c-f4251c3880c5_events.txt", col_names = T) %>%
  mutate(current = samples) %>%
  rowwise() %>%
  mutate(base = substr(reference_kmer, start = 3, stop = 3))

# split data so that each signal (col 14) has a new line
b <- separate_rows(a, current, convert = TRUE)
b <- b %>% mutate(id = seq(1:nrow(b))) %>%  # add time data
  mutate(time = id / 3200) %>% # convert each observation to a time value by dividing by 3200
  mutate(event_status = case_when(
    model_kmer == "NNNNN" ~ "Missed event", # mutate to add 'missed_event' label to certain kmers
    model_kmer != "NNNNN" ~ "Aligned event")) # otherwise, label 'aligned event'

# plot the current from 0 to 2 seconds 
# one type of plot, where color dictates base and vertical lines dictate resquiggle failure
dev.new(noRStudioGD = TRUE)
ggplot(b %>% arrange(event_status) %>% mutate(time = time - 3), aes(x = time, y = current)) +
  geom_point(size = 0.4, alpha = 0.6) +
  theme(text=element_text(size=24)) +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) +
  geom_line(aes(x = time, y = current), alpha = 0.1) +
  xlab("Time (seconds)") +
  ylab("Current (mA)") +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) + 
  xlim(0.55,0.85)

# iterate the above so color is only for missed events
# one type of plot, where color dictates base and vertical lines dictate resquiggle failure
dev.new(noRStudioGD = TRUE)
ggplot(b %>% arrange(event_status) %>% mutate(time = time - 3) %>% mutate(Base = base), aes(x = time, y = current)) +
  geom_point(aes(color = Base), size = 0.4, alpha = 0.6) +
  theme(text=element_text(size=24)) +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) +
  geom_line(aes(x = time, y = current), alpha = 0.1) +
  xlab("Time (seconds)") +
  ylab("Current (mA)") +
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) + 
  xlim(0.5,1)
