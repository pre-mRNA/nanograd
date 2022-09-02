#!/bin/user/env/R

# Aim: Use linear regression to estimate transcript-specific degradation rates based on pooled TIN scores from Agin's 3 experimental conditions
# Written by AS, KW, BK, AR on 2022-08-24
# Last modified by AS on 2022-08-24

# undegraded sample heated for 0 mintues 
# mild-degraded sample heated for 2.5 minutes 
# harsh-degraded sample heated for 7.5 minutes

###########################

# load packages 
library(MASS)
library(tidyverse)

# load in the data 
tin_scores <- read_tsv("./data/2022-08-24_degradationFirst6-poolTin-filt12.txt.gz", col_names = T, col_types = "fddd")

# do simple linear modelling for each transcript 
tin_models <- tin_scores %>% 
  gather('month', 'value', -transcript_id) %>% 
  group_by(transcript_id) %>% 
  mutate(n = row_number()) %>% 
  arrange(transcript_id, n) %>% 
  nest() %>% 
  mutate(model = map(data, fit_model)) %>% 
  mutate(slope = map_dbl(model, get_slope)) %>% 
  dplyr::select(transcript_id, slope)


###########################

# use old code to get transcript metadata

# write a function to import data 
import_tin <- function(dataPath, mySample, myCondition){
  
  # read the tsv 
  read_tsv(dataPath, col_names = T, col_types = "fddfffddddd")
}

# read in the data 
merged_tin <- bind_rows(import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.undegraded_hek293_pass1.fastq.gz.sorted_nanograd5_out.txt.tmp", "control_rep1", "undegraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.undegraded_hek293_pass2.fastq.gz.sorted_nanograd5_out.txt.tmp", "control_rep2", "undegraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/mild_degradataion_rep1.fastq_nanograd5_out.txt.tmp", "mild_rep1", "mildly_degraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/mild_degradataion_rep2.fastq_nanograd5_out.txt.tmp", "mild_rep2", "mildly_degraded"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.5mM_MgCl_degrdation_pass1.fastq.gz.sorted_nanograd5_out.txt.tmp", "heavy_rep1", "heavy_degradation"),
                        import_tin("~/localGadiData/2022-08-23_nanograd5-degradation-first6/all.5mM_MgCl_degrdation_pass2.fastq.gz.sorted_nanograd5_out.txt.tmp", "heavy_rep2", "heavy_degradation"))

annotation <- merged_tin %>% dplyr::select(transcript_id, gene_name, gene_id, transcript_biotype, tx_len, cds_len, utr5_len, utr3_len) %>% unique()

###########################

# add the annotation data to the linear models 
annotated_slopes <- left_join(tin_models, annotation, by = "transcript_id")

###########################

# is degradation associated with transcript_length 
cor.test(annotated_slopes$slope, annotated_slopes$tx_len)

# correlation of 0.2, p < 2.2e-16 
# correlation isn't massive 
# option A: other factors specify the rate of deradation 
# option B: something is going on is that we haven't figued out 


# plotting section 



# plot 1. distribution of slopes 
dev.new(width = 8, height = 8, noRStudioGD = TRUE, unit = "cm")

ggplot(annotated_slopes, aes(x = slope, color = transcript_biotype, fill = transcript_biotype)) + geom_histogram() + facet_grid(~transcript_biotype) + scale_y_log10()

# plot 3. slope against transcript length 
ggplot(annotated_slopes, aes(x = tx_len, y = slope, color = transcript_biotype)) + geom_point() + 
  geom_smooth(method = "lm") + 
  geom_smooth(method = "loess")

filt_slopes <- annotated_slopes %>% filter(tx_len < 2200 & tx_len > 1800)
ggplot(filt_slopes, aes(x = tx_len, y = slope, color = transcript_biotype)) + 
  geom_point() 

# plot 2. 
ggplot(annotated_slopes, aes(x = transcript_biotype, y = slope, color = transcript_biotype, fill = transcript_biotype)) + geom_violin()

###########################


# tutorial section 
# from https://community.rstudio.com/t/linear-regression-by-rows/21452/3

library(dplyr)
library(tidyr)
library(purrr)
library(broom)

data <- data.frame(stringsAsFactors=FALSE,
                   ITEM = c("item A", "item B"),
                   JAN = c(10, 30),
                   FEB = c(20, 35),
                   MAR = c(15, 37),
                   APR = c(17, 38)
)

fit_model <- function(df) lm(value ~ n, data = df)
get_slope <- function(mod) tidy(mod)$estimate[2]

data %>% 
  gather('month', 'value', -ITEM) %>% 
  group_by(ITEM) %>% 
  mutate(n = row_number()) %>% 
  arrange(ITEM, n) %>% 
  nest() %>% 
  mutate(model = map(data, fit_model)) %>% 
  mutate(slope = map_dbl(model, get_slope))