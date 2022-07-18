#!/bin/bash

# written by AJ Sethi on 2022-07-17
# aim: 

####################################################################################################
####################################################################################################

# modules
library(tidyverse)
library(edgeR)

####################################################################################################

# import and preprocess the counts data from featureCounts + uLTRA 

# path to counts file 

# JCSMR iMac: 
# counts_file <- "/Users/AJlocal/localGadiData/2022-06-22_HEK293-degradation-first4-AR_liqa-genome-alignments_BK/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

# Bhavika's laptop: 
# counts_file <- "d:/Users/Sujata Kumar/Desktop/Project/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

# read in the counts file 
raw_counts <- read_tsv(counts_file, col_names = T, skip = 1, col_types = "fccccdfdddd") %>%
  dplyr::rename(gene_id = Geneid, deg_rep1 = 8, deg_rep2 = 9, wt_rep1 = 10, wt_rep2 = 11) %>%
  mutate(across(c(Start, End, Chr, Strand), gsub, pattern = ";.*", replacement = "")) %>% 
  mutate(chr = Chr) %>% 
  select(gene_id, gene_name, chr, wt_rep1, wt_rep2, deg_rep1, deg_rep2)

####################################################################################################

# use LIQA output to calculate the most likely isoform and transcript length for each gene 

liqa <- read_tsv("/Users/AJlocal/localGadiData/2022-07-17_LIQA_isoform_counts/undegraded_hek293_pass1_primary.txt", col_names = T, col_types = "ffddd") 

# from bioMart, for humans, download a table of gene id, transcript_id, and transcirpt length (including UTRs and CDS)

####################################################################################################




# creating a list with samples and  gene_id and length
y <- DGEList(counts=raw_counts[,8:11], genes=raw_counts("gene_id","Length"))

# Normalizing the data
y <- calcNormfactors(y)

