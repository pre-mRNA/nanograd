#!/bin/user/env/R

# Aim: Compare WT TIN with published RNA half lives to see if TIN is affected by RNA turnover 
# Written by AS on 2022-09-08
# Last modified by AS on 2022-09-08

# Published half-live data at https://codeocean.com/capsule/7351682/tree/v1

###########################

# load libraries 

library(tidyverse)
library(rtracklayer)

###########################


# write a function to import data 
import_tin <- function(dataPath, mySample, myCondition){
  
  # read the tsv 
  read_tsv(dataPath, col_names = T, col_types = "fddfffddddd") %>% 
    dplyr::rename(tin = df_ratio) %>% 
    # dplyr::select(transcript_id, tin) %>% 
    mutate(condition = myCondition, sample = mySample)
}

# read in the data 
merged_tin <- bind_rows(import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.undegraded_hek293_pass1.fastq.gz.sorted_nanograd5_out.txt.tmp", "control_rep1", "undegraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.undegraded_hek293_pass2.fastq.gz.sorted_nanograd5_out.txt.tmp", "control_rep2", "undegraded")) %>% 
  filter(cov > 10) %>% 
  dplyr::select(transcript_id, tin, sample, condition)

merged_tin

###########################

# map transcripts to genes 
gene_transcript_map <- rtracklayer::import("~/localGadiData/2022-05-05_KW-Salmon-DTU/Homo_sapiens.GRCh38.104.chr.gtf") %>%
  as_tibble() %>%
  dplyr::select(gene_id, transcript_id) %>%
  na.omit() %>%
  dplyr::distinct()

# make a wide matrix of TIN 
tin_matrix <- merged_tin %>% 
  select(transcript_id, tin, sample) %>% 
  pivot_wider(names_from = sample, values_from = tin) 

genewise_tin <- inner_join(tin_matrix, gene_transcript_map, by = "transcript_id") %>% 
  group_by(gene_id) %>% 
  mutate(n = n()) %>% 
  filter(n == 1) %>% 
  ungroup() %>% 
  select(-n) %>% 
  select(gene_id, control_rep1, control_rep2)

###########################

# import our published half lives 
public_half_lives <- read_csv("./data/published_half_lives.csv", col_names = T) %>% 
  dplyr::rename("gene_id" = 1)

library(corrplot)

# make correlation plots for common transcripts 
dev.new(width = 9, heigth = 9, noRStudioGD = T, unit = "cm")
public_half_lives %>% 
  select(-gene_id) %>% 
  stats::cor(., use = "pairwise.complete.obs", method = "spearman") %>% 
  corrplot(., method = "pie", type = "full", order = "original", tl.col = "black", cl.lim=c(0,1), diag = F, addrect = T)

# check the correlation between libraries 
cor_data <- tin_matrix_complete %>% 
  select(-transcript_id) %>% 
  cor(., method = "spearman") 

###########################

# filter for published hek293 gene half lives 
hek293_known_half_lives <- public_half_lives %>% 
  select(gene_id, contains("293T")) %>% 
  na.omit()

hist(hek293_known_half_lives$Wu_et_al_293T)

tmp <- left_join(genewise_tin, hek293_known_half_lives, by = "gene_id")

dev.new(width = 9, heigth = 9, noRStudioGD = T, unit = "cm")

ggplot(tmp, aes(x = control_rep1, y = Wu_et_al_293T)) + 
  geom_point(size = 2, alpha = 0.5) + 
  geom_smooth(method = "lm") + 
  xlab("HEK293 TIN") + 
  ylab("Published gene decay rate")

cor.test(tmp$control_rep1, tmp$Wu_et_al_293T, method = "spearman")
