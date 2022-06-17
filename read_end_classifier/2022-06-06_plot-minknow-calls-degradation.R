#!/bin/user/env/R

# Aim: Plot the minknow read end classification calls we previously extracted 
# Written by AJ Sethi on 2022-06-06

##################################################
##################################################

library(tidyverse)

# import data 
a <- read_delim("/Users/AJlocal/Documents/Gadi/toLocal/ss1.txt", col_names = F, delim = " ", col_types = "nf") %>% 
  rename(state = 2, count = 1) %>% 
  mutate(sample = "a") %>% 
  mutate(percent = 100 * count / sum(count))

b <- read_delim("/Users/AJlocal/Documents/Gadi/toLocal/ss2.txt", col_names = F, delim = " ", col_types = "nf") %>% 
  rename(state = 2, count = 1) %>% 
  mutate(sample = "b") %>% 
  mutate(percent = 100 * count / sum(count))

merged <- bind_rows(a, b)

# plot 
ggplot(merged, aes(x = state, y = percent, fill = sample)) + 
  geom_bar(position="dodge", stat='identity') + 
  theme_light() + 
  ggtitle(str_wrap("Read end classifications in sample A and B", 40)) +
  xlab("Minknow read end classifications") + 
  ylab("Percentage of reads") 
