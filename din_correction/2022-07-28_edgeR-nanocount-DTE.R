#!/bin/user/env/R

# analyse liqa output from degradation experiment
# written by AS on 28-7-2022


# download data from gadi using a shell command
# scp -r as7425@gadi.nci.org.au:/g/data/lf10/as7425/nanograd/analysis/2022-07-28_align-transcriptome-nanocount/*txt $(pwd)

###########################

# load libraries 

library(tidyverse)
library(edgeR)
library(EnhancedVolcano)

###########################

# read in RNA-Seq libraries 

readNanoCount <- function(path, sample_id, treatment){
  return(read_tsv(path, col_names = T, col_types = "cddd") %>% 
           mutate(sample = sample_id, condition = treatment))
}

nanocount <- bind_rows(readNanoCount("/Users/AJlocal/localGadiData/2022-07-28_align-transcriptome-nanocount/isoformCounts_undegraded_hek293_pass1.txt", "WT_rep1", "WT"),
readNanoCount("/Users/AJlocal/localGadiData/2022-07-28_align-transcriptome-nanocount/isoformCounts_undegraded_hek293_pass2.txt", "WT_rep2", "WT"),
readNanoCount("/Users/AJlocal/localGadiData/2022-07-28_align-transcriptome-nanocount/isoformCounts_5mM_MgCl_degrdation_pass1.txt", "deg_rep1", "deg"),
readNanoCount("/Users/AJlocal/localGadiData/2022-07-28_align-transcriptome-nanocount/isoformCounts_5mM_MgCl_degrdation_pass2.txt", "deg_rep2", "deg"))
