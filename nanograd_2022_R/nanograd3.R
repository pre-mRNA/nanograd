#!/usr/bin/env Rscript 

# import shell arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("\nUsage: Rscript liftoverCustomTranscriptome.R /path/to/annotation.gtf /path/to/alignment.bam /path/to/output.txt", call.=FALSE)
}

library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)
library(rsamtools)

# interactive args for testing (B iMac)
# args <- c("/Users/asethi/localGadiData/2022-01-31_make-annotate-2/annotate_out.bed.temp.gtf", "/Users/asethi/localGadiData/2022-01-31_make-annotate-2/annotate_out.bed.tempbed", "/Users/asethi/localGadiData/2022-01-31_make-annotate-2/annotate_out.bed", "pval")

################################################################################
################################################################################
################################################################################

# start

tstart <- print(paste("start time is", as.character(Sys.time())))

##################################################

# process the annotation and make a reference table with transcript structure information

# read in reference transcripts
gtf <- makeTxDbFromGFF(file=args[1], format = "gtf")

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
tx_biotype <- rtracklayer::import(args[1]) %>%
  as_tibble() %>%
  dplyr::select(transcript_id, transcript_biotype, gene_name, gene_id) %>%
  na.omit() %>%
  dplyr::distinct()

# merge the biotypes with the transcript segment lengths
merged_metadata <- inner_join(tx_biotype, txlen, by = "transcript_id")

##################################################

# calculate bam coverage, assuming each reference sequence in the bam represents a given transcript 

# read in the bam file
bam <- scanBam(bamfile)[[1]] # the result comes in nested lists
# filter reads without match position
ind <- ! is.na(bam$pos)
## remove non-matches, they are not relevant to us
bam <- lapply(bam, function(x) x[ind])
ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
## names of the bam data frame:
## "qname"  "flag"   "rname"  "strand" "pos"    "qwidth"
## "mapq"   "cigar"  "mrnm"   "mpos"   "isize"  "seq"    "qual"
## construc: genomic ranges object containing all reads
ranges <- GRanges(seqnames=Rle(bam$rname), ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname )
## returns a coverage for each reference sequence (aka. chromosome) in the bam file
return (mean(coverage(ranges)))    

