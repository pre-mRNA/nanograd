library(tidyverse)

import <- read_delim("/Users/asethi/localGadiData/2022-03-15_develop-nanograd4/test_out.txt.temp", delim = " ", col_names = F) %>% 
  rename(read = 1, transcript = 2, start = 3, cigar = 4, seq_len = 5) %>% 
  select(-cigar)

# first, find the most abundant transcript 
txabund <- import %>% group_by(transcript) %>% summarise(n = n()) %>% arrange(desc(n))

# plot individual reads 

# rps21
ggplot(import %>% filter(transcript == "ENSMUST00000148629.2"), aes(x = start)) + geom_histogram() + scale_y_log10()

# mt rRNA
ggplot(import %>% filter(transcript == "ENSMUST00000082390.1"), aes(x = start)) + geom_histogram() + scale_y_log10()

