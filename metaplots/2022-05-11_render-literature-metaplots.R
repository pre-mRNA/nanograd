#!/usr/bin/env Rscript

# written by AJ Sethi on 2022-05-11

# Aim: To make metaplots of transcript coverage in published studies, using samtools depth output as the metaplot input
# recycle some of the CHEUI annotate code to calculate metaplot positions

# provide metaplot data
sdepth_out <- "/Users/asethi/localGadiData/2022-05-11_nanograd-metaplots-R/2022-11-05_samtools-depth.txt.gz"
# gadi
# sdepth_out <- "/g/data/lf10/as7425/nanograd/analysis/2022-05-11_literature_metaplots/metaplot_input/2022-11-05_samtools-depth.txt.gz"

# also provide the human annotation in GTF format
human_annotation <- "/Users/asethi/localGadiData/2022-05-11_nanograd-metaplots-R/Homo_sapiens.GRCh38.105.chr.gtf"
# gadi:
# human_annotation <- "/g/data/lf10/as7425/nanograd/analysis/2022-05-11_literature_metaplots/metaplot_input/GRCh38.104.metadata.txt"

####################################################################################################
####################################################################################################

# load libraries
library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)

# process the annotation and make a reference table with transcript structure information

# read in reference transcripts
gtf <- makeTxDbFromGFF(file=human_annotation, format = "gtf")

# make an exon database from the reference transcripts
exons <- exonsBy(gtf, "tx", use.names=TRUE)

# convert the transcripts to a tibble
exons_tib <- as_tibble(as(exons, "data.frame"))

# fetch the length of each transcript segment from the gtf
txlen <- transcriptLengths(gtf, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE) %>%
  as_tibble() %>%
  mutate(diff = tx_len - cds_len - utr5_len - utr3_len) %>%
  dplyr::rename(transcript_id = tx_name) %>%
  dplyr::select(-tx_id, -nexon, -gene_id, -diff)

# the last command doesn't store biotype, so we read in the gtf again using another package
tx_biotype <- rtracklayer::import(human_annotation) %>%
  as_tibble() %>%
  dplyr::select(transcript_id, transcript_biotype, gene_name, gene_id, transcript_version) %>%
  na.omit() %>%
  dplyr::distinct()

# merge the biotypes with the transcript segment lengths
merged_metadata <- inner_join(tx_biotype, txlen, by = "transcript_id") %>%
  unite("transcript", transcript_id, transcript_version, sep = ".")

write_tsv(merged_metadata, "/Users/asethi/localGadiData/2022-05-11_nanograd-metaplots-R/GRCh38.104.metadata.txt", col_names = T)

rm(exons, exons_tib, tx_biotype, txlen)

# read back in
merged_metadata <- read_tsv("/Users/asethi/localGadiData/2022-05-11_nanograd-metaplots-R/GRCh38.104.metadata.txt", col_names = T, col_types = "ffffnnnn")

# gadi
# merged_metadata <- read_tsv("/g/data/lf10/as7425/nanograd/analysis/2022-05-11_2022-05-11_literature_metaplots/metaplot_input/GRCh38.104.metadata.txt", col_names = T, col_types = "ffffnnnn")


##################################################

# import samtools depth
depth_data <- read_tsv(sdepth_out, col_names = T, col_types = "fnnnnnnnn") %>%
  dplyr::rename(transcript = 1, pos = 2, RNA_DRR178487 = 3, RNA_DRR178490 = 4, RNA_ERR4352443 = 5, RNA_ERR4352444 = 6, RNA_NA12878 = 7, RNA_SRR9646141 = 8, RNA_SRR9646143 = 9)

depth_data

##################################################

# attach the transcript lengths and biotype to output using the transcriptome coordinate
merge_out <- inner_join(depth_data, merged_metadata %>% dplyr::rename(transcript = transcript_id), by = "transcript") %>%
  filter(transcript_biotype == "protein_coding") %>%
  dplyr::rename(tx_coord = pos) %>%
  mutate(cds_start = utr5_len,
         cds_end = utr5_len + cds_len,
         tx_end = cds_end + utr3_len) %>%
  mutate(rel_pos = ifelse(tx_coord < cds_start, # if the site is in the 5' utr
                          ((tx_coord)/(cds_start)), # the relative position is simply the position of the site in the UTR
                          ifelse(tx_coord < cds_end, # else, if the site is in the CDS
                                 (1 + (tx_coord - utr5_len)/(cds_len)), # the relative position is betwee 1 and 2 and represents the fraction of the cds the site is in
                                 (2 + (tx_coord - utr5_len - cds_len) / utr3_len))),  # no final condition, the site must be in the 3' utr, similar as before but the rel_pos is between 2 and 3
         abs_cds_start = tx_coord - cds_start, # absolute pos relative to CDS start
         abs_cds_end = tx_coord - cds_end) %>% # absolute pos relative to CDS end
  dplyr::select(1, 2:7, rel_pos)

#%>%
  pivot_longer(cols = starts_with("RNA"), names_to = "library", values_to = "count") %>%
  filter(count > 0) %>%
  dplyr::select(library, rel_pos, count) %>%
  type_convert(col_type = "ffn") %>%
  filter(count > 20) %>%
  uncount(count)

ggplot(merge_out %>% slice_sample(n = 100000000), aes(x = rel_pos, group = library, color = library)) +
  geom_density(aes(y=0.05 * ..count..))

ggplot(merge_out %>% slice_sample(n = 10000000), aes(x = rel_pos, color = library)) +
  geom_histogram(aes(bindwidth = 0.05))



mutate(interval = cut(rel_pos, breaks, include.lowest = TRUE, right = TRUE)) %>%
  group_by(library, interval) %>%
  summarise(sum = sum(count)) %>%
  ungroup()

ggplot(merge_out, aes(x = interval, y = sum)) +
  geom_bar(aes(color = library))












merge_out
##################################################
