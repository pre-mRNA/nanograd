

# preprocess using 
# for i in *; do cat $i | tr " " "\t" |  awk '{$1=$1;print}' | tr " " "\t" | tr " " "\t" |  awk '{$1=$1;print}' | tr " " "\t" > "${i}_fixed.txt"


# read a file 
setwd("/Users/asethi/localGadiData/2022-07-28_align-transcriptome-nanocount")
library(tidyverse)

deg_1 <- read_tsv("counts_filtAlignments_sorted_5mM_MgCl_degrdation_pass1.bam.txt_fixed.txt", col_names = F) %>% 
  rename(count = 1, transcript = 2) %>% 
  mutate(sample = "degraded_rep1", condition = "degraded")

deg_2 <- read_tsv("counts_filtAlignments_sorted_5mM_MgCl_degrdation_pass2.bam.txt_fixed.txt", col_names = F) %>% 
  rename(count = 1, transcript = 2) %>% 
  mutate(sample = "degraded_rep1", condition = "degraded")

wt_1 <- read_tsv("counts_filtAlignments_sorted_undegraded_hek293_pass1.bam.txt_fixed.txt", col_names = F) %>% 
  rename(count = 1, transcript = 2) %>% 
  mutate(sample = "undegraded_rep1", condition = "undegraded")

wt_2 <- read_tsv("counts_filtAlignments_sorted_undegraded_hek293_pass2.bam.txt_fixed.txt", col_names = F) %>% 
  rename(count = 1, transcript = 2) %>% 
  mutate(sample = "undegraded_rep1", condition = "undegraded")

# bind rows
import <- full_join(wt_1, wt_2, by = "transcript") %>% 
  full_join(., deg_1, by = "transcript") %>% 
  full_join(., deg_2, by = "transcript") %>% 
  select(transcript, count.x, count.y, count.x.x, count.y.y) %>% 
  rename(wt_1 = count.x, wt_2 = count.y, deg_1 = count.x.x, deg_2 = count.y.y) %>% 
  replace(is.na(.), 0)
