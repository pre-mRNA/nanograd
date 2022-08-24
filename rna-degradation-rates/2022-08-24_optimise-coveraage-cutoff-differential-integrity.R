#!/bin/user/env/R

# Aim: Optimize the correlation of technical replicates against the number of transcripts in the analysis when pooling techical replicates for the degradation analysis
# Share data with Bhavika and Kat for heatmaps and RNA degradation rates
# Written by AS on 2022-08-24
# Last modified by AS on 2022-08-24

# download data from Gadi
# at /g/data/lf10/as7425/nanograd/analysis/2022-08-23_nanograd5-degradation-first6

###########################

# load libraries 

library(tidyverse)

###########################

# import raw data and prepare for optimization 

# write a function to import data 
import_tin <- function(dataPath, mySample, myCondition, myRep){
  
  # read the tsv 
  read_tsv(dataPath, col_names = T, col_types = "fddfffddddd") %>% 
    dplyr::rename(tin = df_ratio) %>% 
    dplyr::select(transcript_id, tin, cov) %>% 
    mutate(condition = myCondition, sample = mySample, rep = myRep)
}

# read in the data 
merged_tin <- bind_rows(import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.undegraded_hek293_pass1.fastq.gz.sorted_nanograd5_out.txt.tmp", "control_rep1", "undegraded", "rep1"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.undegraded_hek293_pass2.fastq.gz.sorted_nanograd5_out.txt.tmp", "control_rep2", "undegraded", "rep2"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/mild_degradataion_rep1.fastq_nanograd5_out.txt.tmp", "mild_rep1", "mildly_degraded", "rep1"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/mild_degradataion_rep2.fastq_nanograd5_out.txt.tmp", "mild_rep2", "mildly_degraded", "rep2"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.5mM_MgCl_degrdation_pass1.fastq.gz.sorted_nanograd5_out.txt.tmp", "heavy_rep1", "heavy_degradation", "rep1"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.5mM_MgCl_degrdation_pass2.fastq.gz.sorted_nanograd5_out.txt.tmp", "heavy_rep2", "heavy_degradation", "rep2"))

# exploratory data plot 
ggpplot(merged_din)

###########################

# using a loop for different coverage cutoffs, test the number of transcripts remaining against the correlation of technical replicates 
# Investigate the results using different plotting approaches 

# the first thing is to have a single output for correlation 

iterate_cutoff <- function(cutoff){

  # filtering the data
  filt <- merged_tin %>% 
    select(transcript_id, tin, rep, condition, cov) %>% 
    group_by(transcript_id, condition) %>% 
    mutate(cov = sum(cov)) %>% 
    filter(cov > cutoff)
  
  # calculating correlation 
  cor_data <- filt %>% 
  pivot_wider(names_from = rep, values_from = tin) %>% 
  na.omit() %>% 
  select(rep1, rep2)
  cor_out <- cor.test(cor_data$rep1, cor_data$rep2)[[4]] %>% as_tibble()
  
  # calculating number of transcripts that are observed in all conditions 
  transcript_count <- filt %>% 
    select(transcript_id, rep, condition, tin) %>% 
    pivot_wider(names_from = c(condition), values_from = tin) %>% 
    group_by(transcript_id) %>% 
    summarise(undegraded = sum(undegraded), mildly_degraded = sum(mildly_degraded), heavy_degradation = sum(heavy_degradation)) %>% 
    na.omit() %>% 
    nrow()
  
  return(c(cutoff, cor_out, transcript_count) %>% unlist() %>% as_tibble())
}

# let's make a range of cutoffs to test 
cutoff_range <- seq(0,20,1) %>% as_tibble()

# apply the loop to our range of cutoffs 
cutoff_output <- apply(cutoff_range, 1, iterate_cutoff) %>% 
  bind_cols() %>% 
  add_rownames() %>% 
  gather(var, value, -rowname) %>% 
  spread(rowname, value) %>% 
  select(-var) %>% 
  rename(cutoff = 1, cor = 2, tx_count = 3)


# plot the outputs 
dev.new(width = 9, height = 9, noRStudioGD = TRUE, unit = "cm", scale = 2)

ggplot(cutoff_output, aes(x = tx_count, y = cor, color = cutoff)) + 
  geom_point() + 
  ggtitle("Pooled coverage cutoff vs number of transcripts") + 
  xlab("Number of transcripts in all of the 3 conditions") + 
  ylab("Spearman correlation of transcript integrity between technical replicates") + 
  theme_minimal()

ggsave("./plots/2022-08-24_differential-integrity-pool-replicates-cutoff.png", plot = last_plot(), height = 15, width = 15, unit = "cm")


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