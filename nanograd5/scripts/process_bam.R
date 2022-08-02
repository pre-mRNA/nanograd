# make exploratory plots to see if we should use sequence read length or alignment length
# input is a 4-col vector consisting of:
# read ID, transcript, read sequence, cigar string

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("\nUsage: Rscript process_bam.R /path/to/bar.gtf /path/to/bamdata.txt /path/to/output.txt", call.=FALSE)
}

suppressPackageStartupMessages({
library(tidyverse)
library(GenomicFeatures)
library(rtracklayer)
})

input <- read_tsv(args[2], col_names = F, col_types = "ffcc") %>%
  dplyr::rename(read = 1, tx = 2, seq = 4, cig = 3) %>%
  mutate(seq_len = nchar(seq)) %>%
  dplyr::select(-seq)

############################################################
############################################################
############################################################

# parse CIGAR
matcher <- function(pattern, x) {
  ind = gregexpr(pattern, x)[[1]]
  start = as.numeric(ind)
  end = start + attr(ind, "match.length")- 2
  apply(cbind(start,end), 1, function(y) substr(x, start=y[1], stop=y[2]));
}

doone <- function(c, cigar) {
  pat <- paste("\\d+", c , sep="")
  sum(as.numeric(matcher(pat, cigar)), na.rm=T)
}

## takes a cigar string and parses it, not very fast but...
cigarsums <- function(cigar, chars=c("M")) {
  sapply (chars, doone, cigar)
}

# calculate matches
addcig <- input %>%
  rowwise %>% mutate(map_len = cigarsums(cig)) %>%
  dplyr::select(-cig)

rm(input)

# calculate transcript abundance
txabund <- addcig %>%
  group_by(tx) %>%
  mutate(nread = n()) %>%
  arrange(desc(nread)) %>%
  ungroup()

# calculate decay factor
decay <- txabund %>%
  dplyr::select(-read, -nread) %>%
  group_by(tx) %>%
  summarise(df = median(map_len), cov = n())

############################################################
############################################################
############################################################

# process the annotation
anno <- args[1]

# read in reference transcripts
gtf <- makeTxDbFromGFF(file=anno, format = "gtf")

# make an exon database from the reference transcripts
exons <- exonsBy(gtf, "tx", use.names=TRUE)

# convert the transcripts to a tibble
exons_tib <- as_tibble(as(exons, "data.frame"))

# fetch the lengthe of each transcript segment from the gtf
txlen <- transcriptLengths(gtf, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE) %>%
  as_tibble() %>%
  mutate(diff = tx_len - cds_len - utr5_len - utr3_len) %>%
  dplyr::rename(transcript_id = tx_name) %>%
  dplyr::select(-tx_id, -nexon, -gene_id, -diff)

# the last command doesn't store biotype, so we read in the gtf again using another package
tx_biotype <- rtracklayer::import(anno) %>%
  as_tibble() %>%
  dplyr::select(transcript_id, transcript_biotype, gene_name, gene_id) %>%
  na.omit() %>%
  distinct()

# merge the biotypes with the transcript segment lengths
merged_metadata <- inner_join(tx_biotype, txlen, by = "transcript_id")

rm(gtf, exons, exons_tib, txlen, tx_biotype, tmp)

# join data to decay
decay2 <- inner_join(decay %>% separate(tx, into=c("transcript_id", "transcript_version"), sep = "([.])", extra = "merge") %>% dplyr::select(-transcript_version),
                     merged_metadata, by = "transcript_id") %>%
  mutate(df_ratio = df / tx_len)

# clean
rm(decay, merged_metadata)

write_tsv(decay2, args[3], col_names = T)
