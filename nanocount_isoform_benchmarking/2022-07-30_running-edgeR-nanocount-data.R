#!/bin/bash

# written by AJ and BK on 2022-07-30

###########################################################

# preprcess data 
# for i in *; do cat $i|tr " " "\t"| awk '{$1=$1;print}' | tr " " "\t" | tr " " "\t" |  awk '{$1=$1;print}' | tr " " "\t" > "${i}_fixed.txt"

# load library
library(tidyverse)
library(edgeR)

# read data
# AJ's Computer
setwd("/Users/asethi/localGadiData/2022-07-28_align-transcriptome-nanocount")

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


# Bhavika's Computer
setwd("//wsl.localhost/Ubuntu/home/bhavika_kumar/localGadiData/2022-07-28_Nanocount-data")

deg_1 <- read_tsv("counts_filtAlignments_sorted_5mM_MgCl_degrdation_pass1.bam.txt_fixed.txt", col_names=F) %>%
  rename(count = 1, transcript = 2)%>%
  mutate(sample = "degraded_rep1", condition = "degraded")

deg_2 <- read_tsv("counts_filtAlignments_sorted_5mM_MgCl_degrdation_pass2.bam.txt_fixed.txt", col_names=F) %>%
  rename(count = 1, transcript = 2)%>%
  mutate(sample = "degraded_rep1", condition = "degraded")

wt_1 <- read_tsv("counts_filtAlignments_sorted_undegraded_hek293_pass1.bam.txt_fixed.txt", col_names=F) %>%
  rename(count = 1, transcript = 2)%>%
  mutate(sample = "undegraded_rep1", condition = "undegraded")

wt_2 <- read_tsv("counts_filtAlignments_sorted_undegraded_hek293_pass2.bam.txt_fixed.txt", col_names=F) %>%
  rename(count = 1, transcript = 2)%>%
  mutate(sample = "undegraded_rep1", condition = "undegraded")

# bind rows
import <- full_join(wt_1, wt_2, by="transcript") %>%
  full_join(., deg_1, by="transcript") %>%
  full_join(., deg_2, by="transcript") %>%
  select(transcript, count.x, count.y, count.x.x, count.y.y) %>%
  rename(wt_1 = count.x, wt_2 = count.y, deg_1 = count.x.x, deg_2 = count.y.y) %>%
  replace(is.na(.), 0)

# convert import to matrix

testCutoff <- function(cutoff){
  
counts_import_matrix <- as.data.frame(import)
data_clean <- counts_import_matrix[,-1]
row.names(data_clean) <- counts_import_matrix[,1]

# log transform the counts
cpm_log <- cpm(data_clean, log =TRUE)

#get median count
median_log2_cpm <- apply(cpm_log, 1, median)

####################################################################################################


  
  # take the expression cutoff from the input 
  expr_cutoff <- cutoff 
  

# filter for transcripts where median log2 expression is greater than cutoff
data_clean <- data_clean[median_log2_cpm > expr_cutoff, ]

# print number of transcripts in the analysis
rows <- nrow(data_clean)
print(paste("you have", rows, "transcripts remaining in the analysis"))

# define the groups for comparison 
group <- c("control","control","degraded", "degraded")

# design the experiment 
design <- model.matrix(~group)
design

# make DGE object 
y <- DGEList(counts=data_clean,group=factor(group))

# normalize counts using edgeR 
y <- calcNormFactors(y)
y$samples

# estimate dispersion 
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

# extract the differential expression data and convert it to a tibble and
DEtib <- lrt[[14]] %>% 
  as_tibble(rownames = "transcript_id") %>% 
  dplyr::rename(p.raw = PValue)
  
# correct the p-value
DEtib$p.adj <- p.adjust(DEtib$p.raw, method = "fdr")

total_transcripts <- nrow(DEtib)
sig_transcripts <- nrow(DEtib %>% filter(p.adj < 0.05 & abs(logFC) > 1))
specificity <- total_transcripts/(total_transcripts + sig_transcripts)
print(paste(total_transcripts, " total transcripts; ", sig_transcripts, " significant transcripts;", specificity, " specificity"))


return(c(cutoff, total_transcripts, sig_transcripts, specificity) %>% as_tibble())
}

##############################################################################################################

# initialize parameters for the loop 

# make a vector of values to run the loop over 
ints <- seq(1, 10, 0.1) %>% as_tibble()

# run the loop 
iteration_loop_out <- apply(ints, 1, testCutoff) %>% 
  bind_cols() %>%  
  add_rownames() %>% 
  gather(var, value, -rowname) %>% 
  spread(rowname, value) %>% 
  dplyr::select(-var) %>% 
  dplyr::rename(cutoff = 1, total_genes = 2, sig_genes = 3, specificity = 4) %>% 
  arrange(cutoff)

# plot the results 
ggplot(iteration_loop_out, aes(x = sig_genes, y = specificity, color = cutoff)) + geom_point()
ggplot(iteration_loop_out, aes(x = total_genes, y = sig_genes, color = cutoff)) + geom_point() + scale_x_log10()

# assign to a variable, uncorrected 
uncorrected_DEtranscripts <- DEtib

##############################################################################################################

# now repeat, but incorporate DIN 

# take the expression cutoff from the input 
expr_cutoff <- 5

# filter for genes where the median log2 expression is greater than the cutoff 
data_clean <- data_clean[median_log2_cpm > expr_cutoff, ]

# print how many transcripts we still have in the analysis 
rows <- nrow(data_clean)
print(paste("you have ", rows, " transcripts remaining in the analysis"))

# define the groups for comparison 
group <- c("control","control","degraded", "degraded")

# add DIN 
din <- c(11.8, 11.8, 8.24, 8.11)

# design the experiment 
design <- model.matrix(~group+din)
design

# make DGE object 
y <- DGEList(counts=data_clean,group=factor(group))

# normalize counts using edgeR 
y <- calcNormFactors(y)
y$samples

# estimate dispersion 
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

# extract the differential expression data and convert it to a tibble 
DEtib <- lrt[[14]] %>% 
  as_tibble(rownames = "transcript_id") %>% 
  dplyr::rename(p.raw = PValue)
  
# correct the p-value
DEtib$p.adj <- p.adjust(DEtib$p.raw, method = "fdr")


total_transcripts <- nrow(DEtib)
sig_transcripts <- nrow(DEtib %>% filter(p.adj < 0.05 & abs(logFC) > 1))
specificity <- total_transcripts/(total_transcripts + sig_transcripts)
print(paste(total_transcripts, " total transcripts; ", sig_transcripts, " significant transcripts;", specificity, " specificity"))



# assign to a variable, corrected 
corrected_DEtranscripts <- DEtib

##############################################################################################################

# join the corrected and uncorrected analyses 
joined_DEtranscripts <- full_join(uncorrected_DEtranscripts, corrected_DEtranscripts, by = "transcript_id", suffix=c("_uncorrected", "_corrected")) %>% 
  mutate(deltaCPM = logCPM_corrected - logCPM_uncorrected, 
         deltapadj = p.adj_corrected -  p.adj_uncorrected,
         deltalogFC = logFC_corrected - logFC_uncorrected)

ggplot(joined_DEtranscripts, aes(x = deltaCPM)) + geom_histogram() + scale_y_log10()
ggplot(joined_DEtranscripts, aes(x = deltaCPM)) + stat_ecdf()

ggplot(joined_DEtranscripts, aes(x = deltalogFC)) + geom_histogram()
ggplot(joined_DEtranscripts, aes(x = deltalogFC)) + stat_ecdf()

ggplot(joined_DEtranscripts, aes(x = deltapadj)) + geom_histogram()
ggplot(joined_DEtranscripts, aes(x = deltapadj)) + stat_ecdf()

ggplot(joined_DEtranscripts, aes(x = deltalogFC, y = deltapadj)) + geom_point() 

ggplot(joined_DEtranscripts, aes(x = deltaCPM, y = deltapadj)) + geom_point()
