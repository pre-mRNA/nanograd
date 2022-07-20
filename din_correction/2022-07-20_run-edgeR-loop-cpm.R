#!/bin/bash

# written by AJ Sethi on 2022-07-17
# aim: 

####################################################################################################
####################################################################################################

# load libraries 
library(tidyverse)
library(edgeR)
library(EnhancedVolcano)

####################################################################################################

# set paths for all external data 

########## featurecounts gene expression data 

# JCSMR iMac: 
counts_file <- "/Users/AJlocal/localGadiData/2022-06-22_HEK293-degradation-first4-AR_liqa-genome-alignments_BK/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

# Bhavika's laptop: 
# counts_file <- "d:/Users/Sujata Kumar/Desktop/Project/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

########## liqa transcript abundance

# JCSMR iMac
liqa_transcript_counts <- "/Users/AJlocal/localGadiData/2022-07-17_LIQA_isoform_counts/undegraded_hek293_pass1_primary.txt"

########## biomart transcript lengths

# JCSMR iMac
bm_tx_lengths <- "/Users/AJlocal/localGadiData/2022-07-17_LIQA_isoform_counts/2022-07-17_biomart-human-transcript-lengths.txt"

####################################################################################################

# import and preprocess the counts data from featureCounts + uLTRA 

# read in the counts file 
raw_counts <- read_tsv(counts_file, col_names = T, skip = 1, col_types = "fccccdfdddd") %>%
  dplyr::rename(gene_id = Geneid, deg_rep1 = 8, deg_rep2 = 9, wt_rep1 = 10, wt_rep2 = 11) %>%
  mutate(across(c(Start, End, Chr, Strand), gsub, pattern = ";.*", replacement = "")) %>% 
  mutate(chr = Chr) %>% 
  select(gene_id, gene_name, chr, wt_rep1, wt_rep2, deg_rep1, deg_rep2)

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

counts_import_matrix <- as.data.frame(counts)
counts_matrix <- counts_import_matrix[,-1]
row.names(counts_matrix) <- counts_import_matrix[,1]

# define the groups for comparison 
group <- c("control","control","degraded", "degraded")

# make DGE object 
d <- DGEList(counts=counts_matrix,group=factor(group))
d

####################################################################################################

# normalize counts using edgeR 

# Normalization
dt <- calcNormFactors(d,method="TMM") # further testing possible 

####################################################################################################

# make a PCA plot 
#plot plotMDS(dt, gene.selection="common")

####################################################################################################

# estimate dispersion and run edgeR 
d1 <- estimateCommonDisp(dt, verbose=T)
d1 <- estimateTagwiseDisp(d1)

# run edgeR
design.mat <- model.matrix(~ 0 + dt$samples$group)
colnames(design.mat) <- levels(dt$samples$group)

# estimate dispersion 
d2 <- estimateGLMCommonDisp(dt,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="auto")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)

# do edgeR using GLM
fit <- glmFit(d2,design.mat)

# plot BCV 
#plot plotBCV(d2)

##############################################################################################################

# get results from edgeR 

# creating the contrast 
raw_result <- glmLRT(fit,contrast=c(1,-1))

# extract the differential expression data and convert it to a tibble and add gene names 
DEtib <- raw_result[[14]] %>% 
  as_tibble(rownames = "gene_id") %>% 
  dplyr::rename(p.raw = PValue) %>% 
  right_join(name_key, ., by = "gene_id")

##############################################################################################################

# create a loop to test different logCPM cutoffs and return the number of significant genes and the specificity 

# first, write and test the loop 
iterate_cpm <- function(cutoff){
  
  a <- DE_highexpression <- DEtib %>% filter(logCPM > cutoff)
  
  a$p.adj <- p.adjust(a$p.raw, method = "fdr")
  
  
  total_genes <- nrow(a)
  sig_genes <- nrow(a %>% filter(p.adj < 0.05 & abs(logFC) > 1))
  specificity <- total_genes/(total_genes + sig_genes)
  # print(paste(total_genes, " total genes; ", sig_genes, " significant genes;", specificity, " specificity"))
  
  return(c(cutoff, total_genes, sig_genes, specificity) %>% as_tibble())
}

# second, make a vector of values to run the loop over 
ints <- seq(1, 10, 0.1) %>% as_tibble()

# run the loop over our values 
iteration_loop_out <- apply(ints, 1, iterate_cpm) %>% 
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

