#!/bin/user/env/R

# Aim: Make metaplots of transcript coverage for Agin's degradation data 
# Written by AS on 2022-08-29
# Last modified by AS on 2022-08-29

# download data from Gadi
# at /g/data/lf10/as7425/nanograd/analysis/2022-08-29_degradation-first6-metaplots/trimmed

###########################

# change data path depending on where doing analysis 

# data_path <- "~/localGadiData/2022-08-29_degradation-coverage-metaplots"
# data_path <- "/g/data/lf10/as7425/nanograd/analysis/2022-08-29_degradation-first6-metaplots/trimmed"

###########################

# load libraries 

library(tidyverse)

###########################

# import raw data and aggregate

# write a function to import data 
import_coverage <- function(dataPath, mySample, myCondition, myRep){
  
  # read the tsv 
  read_tsv(dataPath, col_names = T, col_types = "dfdd", num_threads = 4) %>% 
    rename(coverage = score) %>% 
    mutate(condition = myCondition, sample = mySample)
}

# read in the data 
merged_coverage <- bind_rows(import_coverage(paste(data_path, "/trimmed_annotated_all.undegraded_hek293_pass1.fastq.gz.sorted.bam_coverage.txt.gz", sep = ""), "control_rep1", "undegraded", "rep1"),
                        import_coverage(paste(data_path, "/trimmed_annotated_all.undegraded_hek293_pass2.fastq.gz.sorted.bam_coverage.txt.gz", sep = ""), "control_rep2", "undegraded", "rep2"),
                        import_coverage(paste(data_path, "/trimmed_annotated_primary_mild_degradataion_rep1.fastq.bam_coverage.txt.gz", sep = ""), "mild_rep1", "mild_degradation", "rep1"),
                        import_coverage(paste(data_path, "/trimmed_annotated_primary_mild_degradataion_rep2.fastq.bam_coverage.txt.gz", sep = ""), "mild_rep2", "mild_degradation", "rep2"),
                        import_coverage(paste(data_path, "/trimmed_annotated_all.5mM_MgCl_degrdation_pass1.fastq.gz.sorted.bam_coverage.txt.gz", sep = ""), "heavy_rep1", "heavy_degradation", "rep1"),
                        import_coverage(paste(data_path, "/trimmed_annotated_all.5mM_MgCl_degrdation_pass2.fastq.gz.sorted.bam_coverage.txt.gz", sep = ""), "heavy_rep2", "heavy_degradation", "rep2"))


###########################

# summarise coverage across windows 

# define break width 
breaks <- seq(0,3,0.025) 

# get cumulative coverage 
coverage_breaks <- merged_coverage %>% 
  filter(transcript_biotype == "protein_coding") %>% 
  mutate(interval = cut(rel_pos, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)) %>%
  group_by(interval, condition) %>% 
  summarise(total_coverage = sum(coverage)) 
Sys.time()

# save coverage_breaks for protein coding genes 
write_tsv(coverage_breaks, paste(data_path, "/protein_coding_coverage_breaks.txt", sep = ""), col_names = T)

###########################

# plot coverage in metatranscript bins 
coverage_breaks <- read_tsv(paste(data_path, "/protein_coding_coverage_breaks.txt", sep = ""), col_names = T) %>% 
  group_by(condition) %>% 
  mutate(cum_cov = sum(total_coverage)) %>% 
  ungroup() %>% 
  mutate(norm_cov = total_coverage / cum_cov) %>% 
  mutate(condition = factor(condition, levels = c("undegraded", "mildly_degraded", "heavy_degradation"))) %>% 
  mutate(metagene = case_when(interval < 40 ~ "5UTR",
                              (interval > 40 & interval < 80) ~ "CDS",
                              interval > 80 ~ "3UTR")) %>% 
  mutate(metagene = factor(metagene, levels = c("5UTR", "CDS", "3UTR"))) 

# scatterplot
ggplot(coverage_breaks %>% na.omit(metagene) %>% mutate(interval = interval %% 40), aes(x = interval, y = norm_cov, color = condition)) + 
  geom_point(alpha = 1) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Protein coding transcripts") + 
  xlab("Relative metatranscriptomic location") + 
  ylab("Proportion of transcriptomic coverage") + 
  facet_wrap(~metagene, nrow = 3, ncol = 1, scales="free_y")

# barchart 
ggplot(coverage_breaks, aes(x = interval, y = norm_cov, color = condition, fill = condition)) + 
  geom_bar(alpha = 1, stat = 'identity', position = 'dodge') + 
  geom_vline(xintercept = c(80,40), col = "black", alpha = 0.5) + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Protein coding transcripts") + 
  xlab("Relative metatranscriptomic location") + 
  xlim(40,80)
  ylab("Proportion of transcriptomic coverage") # + 
  geom_smooth(method = "lm", alpha = 0.2) # iterate over loess span 
  geom_smooth(span = 0.35, alpha = 0.2) + # iterate over loess span 










# merge by condition 
ggplot(merged_coverage %>% slice_sample(n = 200000000), aes(x = rel_pos, y = coverage)) +
  geom_bin2d() + 
  geom_vline(xintercept = 1:2, col = "grey") + 
  theme_minimal() + 
  theme(text = element_text(size=24)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Metatranscript coverage density") + 
  xlab("Metatranscript location") + ylab("Coverage") + 
  facet_wrap(~condition, nrow = 3) + 
  xlim(0,3) + 
  scale_y_log10()

ggplot(merged_coverage %>% slice_sample(n = 20000000), aes(x = rel_pos, fill = condition)) +
  geom_histogram() + 
  geom_vline(xintercept = 1:2, col = "grey") + 
  theme_bw() + 
  theme(text = element_text(size=24)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Metatranscript coverage density") + 
  xlab("Metatranscript location") + ylab("Coverage") + 
  facet_wrap(~condition, nrow = 3) 


