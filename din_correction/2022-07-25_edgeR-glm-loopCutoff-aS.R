#!/bin/bash

# written by AJ abd BK on 2022-07-25

###########################################################

# load libraries 
library(tidyverse)
library(edgeR)
library(EnhancedVolcano)

# JCSMR iMac: 
counts_file <- "/Users/AJlocal/localGadiData/2022-06-22_HEK293-degradation-first4-AR_liqa-genome-alignments_BK/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

# Bhavika's laptop: 
# counts_file <- "d:/Users/Sujata Kumar/Desktop/Project/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

########## liqa transcript abundance

# JCSMR iMac
liqa_transcript_counts <- "/Users/AJlocal/localGadiData/2022-07-17_LIQA_isoform_counts/undegraded_hek293_pass1_primary.txt"

# liqa_transcript_counts <- "//wsl.localhost/Ubuntu/home/bhavika_kumar/localGadiData/isoform_expression/undegraded_hek293_pass1_primary.txt"

########## biomart transcript lengths

# JCSMR iMac
bm_tx_lengths <- "/Users/AJlocal/localGadiData/2022-07-17_LIQA_isoform_counts/2022-07-17_biomart-human-transcript-lengths.txt"

# bm_tx_lengths <- "d:/Users/Sujata Kumar/Desktop/Project/2022-07-17_biomart-human-transcript-lengths.txt"

# import and preprocess the counts data from featureCounts + uLTRA 

# read in the counts file 
raw_counts <- read_tsv(counts_file, col_names = T, skip = 1, col_types = "fccccdfdddd") %>%
  dplyr::rename(gene_id = Geneid, deg_rep1 = 8, deg_rep2 = 9, wt_rep1 = 10, wt_rep2 = 11) %>%
  mutate(across(c(Start, End, Chr, Strand), gsub, pattern = ";.*", replacement = "")) %>% 
  mutate(chr = Chr) %>% 
  select(gene_id, gene_name, chr, wt_rep1, wt_rep2, deg_rep1, deg_rep2)


# plot gene coverage histograms 
# raw_counts %>% pivot_longer(cols = contains("_rep"), names_to = "library", values_to = "count") %>% ggplot(., aes(x = count, fill = library)) + geom_histogram() + facet_wrap(~library) + scale_x_log10() + scale_y_log10() + xlab("Count per gene") + ylab("Density")

####################################################################################################

# use LIQA output to calculate the most likely isoform and transcript length for each gene 

liqa <- read_tsv(liqa_transcript_counts, col_names = T, col_types = "ffddd") %>% 
  rename(transcript_id = IsoformName) %>% 
  rename(transcript_abundance = ReadPerGene_corrected) %>% 
  select(transcript_id, transcript_abundance)

# from bioMart, for humans, download a table of gene id, transcript_id, and transcirpt length (including UTRs and CDS)
# uploaded to OneDrive as 2022-07-17_biomart-transcript-lengths.txt

bm <- read_tsv(bm_tx_lengths, col_names = T) %>% 
  rename(gene_id = 1, transcript_id = 3, transcript_length = 5) %>% 
  select(gene_id, transcript_id, transcript_length)

# merge the data 
merged <- left_join(bm, liqa, by = "transcript_id") %>% 
  group_by(gene_id) %>% 
  slice(which.max(transcript_abundance))

# figure out which genes are not in merged
`%ni%` <- Negate(`%in%`)
outside_genes <- bm %>% filter(gene_id %ni% merged$gene_id) %>% 
  group_by(gene_id) %>% 
  summarise(transcript_length = mean(transcript_length))

# merge merged with outside_genes to have a length for each gene 
final_length_key <- bind_rows(merged %>% select(gene_id, transcript_length), outside_genes)

####################################################################################################
# attach gene length to raw_counts 

counts <- inner_join(raw_counts, final_length_key, by = "gene_id") %>% 
  select(gene_id, wt_rep1, wt_rep2, deg_rep1, deg_rep2)

name_key <- inner_join(raw_counts, final_length_key, by = "gene_id") %>% select(gene_id, gene_name)

####################################################################################################

# convert raw counts to matrix format 
testCutoff <- function(cutoff){
  
counts_import_matrix <- as.data.frame(counts)
data_clean <- counts_import_matrix[,-1]
row.names(data_clean) <- counts_import_matrix[,1]

# log transform the counts 
cpm_log <- cpm(data_clean, log = TRUE)

# get median count
median_log2_cpm <- apply(cpm_log, 1, median)

# plot median log cpm 
#plot hist(median_log2_cpm)

####################################################################################################

# create loop for expression cutoff 



  # take the expression cutoff from the input 
  expr_cutoff <- cutoff

  # filter for genes where the median log2 expression is greater than the cutoff 
  data_clean <- data_clean[median_log2_cpm > expr_cutoff, ]

  # print how many genes we still have in the analysis 
  rows <- nrow(data_clean)
  print(paste("you have ", rows, " genes remaining in the analysis"))

  # define the groups for comparison 
  group <- c("control","control","degraded", "degraded")
  
   din <- c(11.8, 11.8, 8.24, 8.11)
  #din <- c(10, 10, 8, 8)
  
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

  # extract the differential expression data and convert it to a tibble and add gene names 
  DEtib <- lrt[[14]] %>% 
    as_tibble(rownames = "gene_id") %>% 
    dplyr::rename(p.raw = PValue) %>% 
    right_join(name_key, ., by = "gene_id")

  # correct the p-value
  DEtib$p.adj <- p.adjust(DEtib$p.raw, method = "fdr")


  total_genes <- nrow(DEtib)
  sig_genes <- nrow(DEtib %>% filter(p.adj < 0.05 & abs(logFC) > 1))
  specificity <- total_genes/(total_genes + sig_genes)
  print(paste(total_genes, " total genes; ", sig_genes, " significant genes;", specificity, " specificity"))
  
  return(c(cutoff, total_genes, sig_genes, specificity) %>% as_tibble())
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

##############################################################################################################

# pick the cutoff 
# 5.5cpm logCPM cutoff, 4721 genes,	61 DEG,	0.987244 specificity

# new cutoff
# 5.0 cpm prefilter,	6463 genes,	0 DEGs,	1.0000000 specificity 











