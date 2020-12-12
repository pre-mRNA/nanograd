# written by A.J. Sethi on 2020-11-21
# calculates statistics about polyA clusters

library(tidyverse)

# parse the arguments
args <- commandArgs()
outdir <- args[7]
infile <- args[8]

# import the data
cImport <- read_tsv(infile, comment = "#", col_names=F, col_types = "fddcdf?") %>%
  rename(chromosome=1, start=2, end=3, name=4, score=5, strand=6, cluster=7)

# count the support level for each cluster
highConfidenceClusters <- cImport %>%
  group_by(cluster) %>%
  summarise(count = n()) %>%
  filter(count > 20) %>%
  select(cluster)

# identify putative polyA sites and put them in BED format -- START
by_cluster <- cImport %>%
  group_by(cluster) %>%
  summarise(position = round(mean(start), 0), chromosome=first(chromosome), strand=first(strand))

# create a bed-styled list of polyA sites
polyA_total_bed <- by_cluster %>%
  mutate(start = position -1 ) %>%
  mutate(end = position + 1) %>%
  add_column(score = "1") %>%
  select(chromosome, start, end, cluster, score, strand)
# identify putative polyA sites and put them in BED format -- END

# filter for bed intervals in the high confidence targets
highConfidenceBed <- subset(polyA_total_bed, cluster %in% highConfidenceClusters$cluster)

# write out the "high confidence clusters" to the output directory
write_tsv(highConfidenceBed, outdir, col_names = F)
