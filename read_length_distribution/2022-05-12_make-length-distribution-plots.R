#!/bin/bash

# written by AJ Sethi on 2022-05-05
# Plot observed read lengths (BAM) vs expected read lenths (transcript annotation) per-library


####################################################################################################
####################################################################################################

# modules
library(tidyverse)

####################################################################################################

# process LIQA data and annotation to get a distribution of expected read lengths 

# write a function to import liqa output 
importE <- function(path, myCondition, mySample){
  read_tsv(path, col_names = T) %>% 
    dplyr::mutate(sample = mySample) %>% 
    dplyr::mutate(condition = myCondition)
}

# import first 4 samples expected read length data from liqa
e_wt_1 <- importE("~/localGadiData/2022-05-12_liqa/undegraded_hek293_pass1_primary.txt", "wt", "wt_1")
e_wt_2 <- importE("~/localGadiData/2022-05-12_liqa/undegraded_hek293_pass2_primary.txt", "wt", "wt_2")
e_deg_1 <- importE("~/localGadiData/2022-05-12_liqa/5mM_MgCl_degrdation_pass1_primary.txt", "deg", "deg_1")
e_deg_2 <- importE("~/localGadiData/2022-05-12_liqa/5mM_MgCl_degrdation_pass2_primary.txt", "deg", "deg_2")

# bind rows 
e_comp <- bind_rows(e_wt_1, e_wt_2, e_deg_1, e_deg_2) %>% 
  dplyr::rename(transcript = IsoformName, 
         count = ReadPerGene_corrected) %>% 
  dplyr::select(transcript, count, condition) %>%
  mutate_at(vars(count), ~ as.integer(round(.x)))

# link human annotation metadata
human_metadata <- read_tsv("~/localGadiData/2022-05-11_nanograd-metaplots-R/GRCh38.104.metadata.txt", col_names = T) %>% 
  select(transcript, tx_len) %>% 
  separate(transcript, into = "transcript", remove = T)

# join observed counts with metadata
e_count_final <- inner_join(e_comp, human_metadata, by = "transcript") %>% 
  uncount(count)

# plot for expected lengths
ggplot(e_count_final, aes(x = tx_len, color = condition)) + geom_density(adjust = 3) + xlim(0,7500)

####################################################################################################

# process BAM alignments to get a distribution of observed read lengths 

# write in a function to import observed read lengths
importO <- function(path, myCondition, mySample){
  read_tsv(path, col_names = F) %>% 
    dplyr::rename(read_len = 1, bam = 2) %>% 
    dplyr::mutate(condition = myCondition) %>% 
    dplyr::mutate(sample = mySample)
}

# import first 4 samples
o_wt_1 <- importO("~/localGadiData/2022-05-12_read-length-distribution/undegraded_hek293_pass1_primary.bam_alignment_length.txt", "wt", "wt_1")
o_wt_2 <- importO("~/localGadiData/2022-05-12_read-length-distribution/undegraded_hek293_pass2_primary.bam_alignment_length.txt", "wt", "wt_2")
o_deg_1 <- importO("~/localGadiData/2022-05-12_read-length-distribution/5mM_MgCl_degrdation_pass1_primary.bam_alignment_length.txt", "deg", "deg_1")
o_deg_2 <- importO("~/localGadiData/2022-05-12_read-length-distribution/5mM_MgCl_degrdation_pass2_primary.bam_alignment_length.txt", "deg", "deg_2")

# bind rows 
o_comp <- bind_rows(o_wt_1, o_wt_2, o_deg_1, o_deg_2)

# plot for observed lengths 
ggplot(o_comp, aes(x = read_len, color = condition)) + geom_density(adjust = 3) + xlim(0,7500)


####################################################################################################

# Merge expected and observed distributions 
eando <- bind_rows(o_comp %>% 
                     mutate(type = "observed") %>% 
                     mutate(type = ifelse(condition == "wt", "Observed read lengths in WT", "Observed read lengths in degraded")),
                   e_count_final %>% 
                     mutate(type = "expected", read_len = tx_len) %>% 
                     filter(condition == "wt") %>% 
                     mutate(type = "Expected read lengths in WT"))


# set factor levels 
eando$type <- factor(eando$type, levels = c("Observed read lengths in WT", "Observed read lengths in degraded", "Expected read lengths in WT"))

# prepare to plot 
dev.new(width=8.5, height=4.5, unit="in")

# plot expected and observed read lengths
g <- ggplot(eando %>% filter(), aes(x = read_len)) + 
  geom_density(adjust = 3, aes(color = type), size = 0.7) + 
  xlim(0,10000) + 
  theme_light() + 
  ggtitle(str_wrap("Expected and observed alignment lengths in HEK293 poly(A)+ RNA", 40)) +
  xlab("Alignment length (nt)") + 
  ylab("Density") + 
  theme(text = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        plot.title= element_text(size = 16, hjust = 0.5), 
        legend.title = element_text(colour="black", size=0, face="bold"),
        legend.spacing.x = unit(0.2, 'cm'), 
        legend.spacing.y = unit(0.15, 'cm'), 
        legend.position=c(0.62,0.6), 
        legend.background=element_blank()) + 
  scale_colour_manual(values = c("#7474E8", "#9E3954", "#ADA81D"))

# save
ggsave("~/Desktop/2022-05-12_read-length-distribution.png", plot = g, scale = 1, width = 4.5, height = 4.5, unit = "in", dpi = 300, bg = NULL)

eando2 <- eando %>% filter() %>% mutate(xlab = str_wrap(eando$type, width = 12))

g1 <- ggplot(eando2, aes(x = xlab, y = read_len)) + 
  geom_boxplot(aes(fill = type)) + 
  theme_light() + 
  ylim(0,5000) + 
  ggtitle(str_wrap("Alignment lengths in HEK293 poly(A)+ RNA", 40)) +
  ylab("Length (nt)") + 
  xlab("sample") + 
  theme(text = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        plot.title= element_text(size = 16, hjust = 0.5), 
        legend.title = element_text(colour="black", size=0, face="bold"),
        legend.spacing.x = unit(0.2, 'cm'), 
        legend.spacing.y = unit(0.15, 'cm'), 
        legend.background=element_blank()) + 
  scale_colour_manual(values = c("#7474E8", "#9E3954", "#ADA81D"))

g1

# save
ggsave("~/Desktop/2022-05-12_histogram.png", plot = g1, scale = 1, width = 4.5, height = 4.5, unit = "in", dpi = 300, bg = NULL)


