#!/bin/user/env/R

# Aim: Merge nanograd5 transcript integrity values from first 6 degradation samples 
# Share data with Bhavika and Kat for heatmaps and RNA degradation rates
# Written by AS on 2022-08-23
# Last modified by AS on 2022-08-23

# download data from Gadi
# at /g/data/lf10/as7425/nanograd/analysis/2022-08-23_nanograd5-degradation-first6

###########################

# load libraries 

library(tidyverse)
library(ComplexHeatmap)

###########################


# write a function to import data 
import_tin <- function(dataPath, mySample, myCondition){
  
  # read the tsv 
  read_tsv(dataPath, col_names = T, col_types = "fddfffddddd") %>% 
    dplyr::rename(tin = df_ratio) %>% 
    dplyr::select(transcript_id, tin) %>% 
    mutate(condition = myCondition, sample = mySample)
}

# read in the data 
merged_tin <- bind_rows(import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.undegraded_hek293_pass1.fastq.gz.sorted_nanograd5_out.txt.tmp", "control_rep1", "undegraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.undegraded_hek293_pass2.fastq.gz.sorted_nanograd5_out.txt.tmp", "control_rep2", "undegraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/mild_degradataion_rep1.fastq_nanograd5_out.txt.tmp", "mild_rep1", "mildly_degraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/mild_degradataion_rep2.fastq_nanograd5_out.txt.tmp", "mild_rep2", "mildly_degraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.5mM_MgCl_degrdation_pass1.fastq.gz.sorted_nanograd5_out.txt.tmp", "heavy_rep1", "heavy_degradation"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.5mM_MgCl_degrdation_pass2.fastq.gz.sorted_nanograd5_out.txt.tmp", "heavy_rep2", "heavy_degradation"))

merged_tin

###########################

# make a wide matrix of TIN 
tin_matrix <- merged_tin %>% 
  select(transcript_id, tin, sample) %>% 
  pivot_wider(names_from = sample, values_from = tin)

# note that we have 55450 transcripts with NAs being present 
# let's filter out NAs 

tin_matrix_complete <- tin_matrix %>% 
  na.omit() 
nrow(tin_matrix_complete)
# 290 transcripts after filtering 

###########################

# transcript level correlation plot 
library(corrplot)

# make correlation plots for common transcripts 
dev.new(width = 9, heigth = 9, noRStudioGD = T, unit = "cm")
tin_matrix_complete %>% 
  select(-transcript_id) %>% 
  cor(., method = "spearman") %>% 
  corrplot(., method = "pie", type = "full", order = "original", tl.col = "black", cl.lim=c(0,1), diag = F, addrect = T)

# check the correlation between libraries 
cor_data <- tin_matrix_complete %>% 
  select(-transcript_id) %>% 
  cor(., method = "spearman") 

###########################

# pool the replicates by condition  

# write a function to import data 
import_tin <- function(dataPath, mySample, myCondition){
  
  # read the tsv 
  read_tsv(dataPath, col_names = T, col_types = "fddfffddddd") %>% 
    dplyr::rename(tin = df_ratio) %>% 
    dplyr::select(transcript_id, tin, cov) %>% 
    mutate(condition = myCondition, sample = mySample)
}

# read in the data 
merged_tin <- bind_rows(import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.undegraded_hek293_pass1.fastq.gz.sorted_nanograd5_out.txt.tmp", "control_rep1", "undegraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.undegraded_hek293_pass2.fastq.gz.sorted_nanograd5_out.txt.tmp", "control_rep2", "undegraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/mild_degradataion_rep1.fastq_nanograd5_out.txt.tmp", "mild_rep1", "mildly_degraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/mild_degradataion_rep2.fastq_nanograd5_out.txt.tmp", "mild_rep2", "mildly_degraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.5mM_MgCl_degrdation_pass1.fastq.gz.sorted_nanograd5_out.txt.tmp", "heavy_rep1", "heavy_degradation"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.5mM_MgCl_degrdation_pass2.fastq.gz.sorted_nanograd5_out.txt.tmp", "heavy_rep2", "heavy_degradation"))

# take the average per condition te
merged_per_condition_tin <- merged_tin %>% 
  group_by(transcript_id, condition) %>% 
  summarise(tin = weighted.mean(tin, cov), sum_cov = sum(cov)) %>%
  filter(sum_cov > 10) %>% # impose an arbitrary cutoff for sum_cov and this needs to be tuned 
  select(-sum_cov) %>% 
  pivot_wider(names_from = condition, values_from = tin) %>% 
  na.omit() %>% 
  dplyr::select(transcript_id, undegraded, mildly_degraded, heavy_degradation)

merged_per_condition_tin # 3516 transcripts 

###########################

# make a complex heatmap 

# we first want to convert to marix format
merged_matrix <- merged_per_condition_tin %>% ungroup() %>% select(-transcript_id) %>% as.matrix()
rownames(merged_matrix) <- merged_per_condition_tin$transcript_id


# perform clustering and make a heatmap 
dev.new(width = 12, height = 12, noRStudioGD = T, unit = "cm")
Heatmap(merged_matrix)

# save the data 
write_tsv(merged_per_condition_tin, "./data/2022-08-23_degradation-first6-perConditionTin.txt.gz", col_names = T)








