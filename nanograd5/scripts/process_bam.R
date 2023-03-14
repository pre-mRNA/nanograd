#!/usr/bin/env Rscript

# about

# make exploratory plots to see if we should use sequence read length or alignment length
# input is a 5-col vector consisting of:
# read ID, transcript, alignment 5' start, cigar string, sequence length

############################################################

# housekeeping

# import positional arguments
args = commandArgs(trailingOnly=TRUE)

### for testing, args <- c("/Users/asethi/localGadiData/2021-12-24_test-custom-metaplot/Mus_musculus.GRCm39.104.chr.gtf", "/Users/asethi/localGadiData/2022-03-15_develop-nanograd4/test_out.txt.temp", "/Users/asethi/localGadiData/2022-03-15_develop-nanograd4/nanograd5_Rout.txt")

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("\nUsage: Rscript process_bam.R /path/to/bar.gtf /path/to/bamdata.txt /path/to/output.txt", call.=FALSE)
}

# load packages quietly
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicFeatures)
  library(rtracklayer)
})

############################################################

# read in the trimmed BAM data

input <- read_delim(args[2], delim = " ", col_names = F, col_types = "ffdcd") %>%
  dplyr::rename(read = 1, tx = 2, start_coord = 3, cig = 4, seq_len = 5)

############################################################

# calculate the aligned length of read by parsing CIGAR matches

# parse CIGAR
matcher <- function(pattern, x) {
  ind = gregexpr(pattern, x)[[1]]
  start = as.numeric(ind)
  end = start + attr(ind, "match.length")- 2
  apply(cbind(start,end), 1, function(y) substr(x, start=y[1], stop=y[2]));
}

# The matcher function takes a regular expression pattern and a string x as input. It uses gregexpr to find all matches of the pattern in x. It then extracts the starting index of each match, and computes the ending index by adding the length of the matched pattern (which is stored as an attribute of ind) and subtracting 2 (since CIGAR strings use a 1-based coordinate system and include a type code at the end of each match). Finally, it uses apply to iterate over each matched substring and extract it from x using substr.

doone <- function(c, cigar) {
  pat <- paste("\\d+", c , sep="")
  sum(as.numeric(matcher(pat, cigar)), na.rm=T)
}

# The doone function takes a character code c (e.g. "M" for match) and a CIGAR string cigar as input. It creates a regular expression pattern for the code by concatenating the string \\d+ (which matches one or more digits) with c. It then calls matcher to extract all matches of the pattern in cigar, and converts the resulting vector of substrings to numeric values. The function returns the sum of these values, ignoring any missing values (indicated by NA) using the na.rm argument of sum.

## takes a cigar string and parses it, not very fast but...
cigarsums <- function(cigar, chars=c("M")) {
  sapply (chars, doone, cigar)
}

# The cigarsums function takes a CIGAR string and a vector of character codes (defaulting to "M") as input. It uses sapply to apply the doone function to each element of the chars vector, passing cigar as the second argument. The function returns a vector of numeric values, giving the sum of all matches of each code in chars in cigar.

# calculate matches
addcig <- input %>%
  rowwise %>% mutate(map_len = cigarsums(cig)) %>%
  dplyr::select(-cig)

# we apply the cigarsums funuction to each row of the input to calculate the total mapped length 

# memory management
rm(input)

# exploratory plots
# ggplot(addcig, aes(x = seq_len, y = map_len, color = tx)) + geom_point() + scale_x_log10() + scale_y_log10() + geom_abline(xintercept = 0)
# ggplot(addcig, aes(x = map_len, y = start_coord, color = tx)) + geom_point(alpha = 0.01) + facet_wrap(~tx)

############################################################

# calculate transcript abundance
decay <- addcig %>%
  group_by(tx) %>%
  mutate(nread = n()) %>%
  ungroup() %>%
  dplyr::select(-read, -nread) %>%
  group_by(tx) %>%
  summarise(median_align = median(map_len), cov = n())

# memory management
rm(addcig)

############################################################

# parse the annotation and merge transcript length and biotype information with the previous length calls

# process the annotation
anno <- args[1]

# read in reference transcripts
gtf <- makeTxDbFromGFF(file=anno, format = "gtf")

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

rm(gtf, txlen, tx_biotype)

############################################################

# merge transcriptome data with the original decay data

# join data to decay
decay2 <- inner_join(decay %>%
                       separate(tx, into=c("transcript_id", "transcript_version"), sep = "([.])", extra = "merge") %>%
                       dplyr::select(-transcript_version),
                     merged_metadata, by = "transcript_id") %>%
  mutate(df_ratio = median_align / tx_len)

# clean
rm(decay, merged_metadata)

write_tsv(decay2, paste(args[3], ".tmp", sep = ""), col_names = T)

############################################################

# calculate DIN and print

DIN <- decay2  %>%
  mutate(df_ratio = median_align / tx_len) %>%
  summarise(sum = sum(df_ratio), n = n()) %>%
  mutate(DIN = 30 * sum/n) %>%
  dplyr::select(DIN)

print(paste("Sample DIN is equal to ", round(DIN, 3)))

a <- c("Sample DIN:", DIN) %>% as.data.frame()

write_tsv(a, args[3], col_names = F)
