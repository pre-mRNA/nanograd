#!/bin/user/env/R

# Aim: Make histograms of read length distribution for Agin's first 6 degradation samples 
# Written by AS on 2022-08-29
# Last modified by AS on 2022-08-29

# download data from Gadi
# at "/g/data/lf10/as7425/nanograd/analysis/2022-08-29_read-length-distribution"


###########################

# get data path for gadi or local use 

data_path <- "~/localGadiData/2022-08-29_nanograd-first6-read-length-distribution"
# data_path <- "/g/data/lf10/as7425/nanograd/analysis/2022-08-29_degradation-first6-metaplots/trimmed"

###########################

# process BAM alignments to get a distribution of observed read lengths 

# write in a function to import observed read lengths
import_length <- function(path, mySample, myCondition, myRep){
  read_tsv(path, col_names = F) %>% 
    dplyr::rename(read_len = 1, bam = 2) %>% 
    dplyr::mutate(condition = myCondition) %>% 
    dplyr::mutate(sample = mySample)
}

# import first 6 samples 
merged_length <- bind_rows(import_length(paste(data_path, "/all.undegraded_hek293_pass1.fastq.gz.sorted.bam_alignment_length.txt", sep = ""), "control_rep1", "undegraded", "rep1"),
                             import_length(paste(data_path, "/all.undegraded_hek293_pass2.fastq.gz.sorted.bam_alignment_length.txt", sep = ""), "control_rep2", "undegraded", "rep2"),
                             import_length(paste(data_path, "/primary_mild_degradataion_rep1.fastq.bam_alignment_length.txt", sep = ""), "mild_rep1", "mild_degraded", "rep1"),
                             import_length(paste(data_path, "/primary_mild_degradataion_rep2.fastq.bam_alignment_length.txt", sep = ""), "mild_rep2", "mild_degraded", "rep2"),
                             import_length(paste(data_path, "/all.5mM_MgCl_degrdation_pass1.fastq.gz.sorted.bam_alignment_length.txt", sep = ""), "heavy_rep1", "heavy_degraded", "rep1"),
                             import_length(paste(data_path, "/all.5mM_MgCl_degrdation_pass2.fastq.gz.sorted.bam_alignment_length.txt", sep = ""), "heavy_rep2", "heavy_degraded", "rep2"))
                             

# plot for observed lengths 
ggplot(merged_length, aes(x = read_len, color = sample)) + 
  geom_density(adjust = 1) + 
  xlim(5000,10000)

