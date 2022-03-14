# make exploratory plots to see if we should use sequence read length or alignment length 
# input is a 4-col vector consisting of: 
# read ID, transcript, read sequence, cigar string 

# test data: 6m mouse reads aligned to a coding + noncoding transcriptome with no pseudoalignments 
# /g/data/lf10/as7425/2020-11_mouseBrain/analysis/2021-09-24_testAlignmentTranscriptome/allReadsTranscriptome_Ensembl.bam

library(tidyverse)

input <- read_tsv("/Users/asethi/localGadiData/2022-03-15_develop-nanograd4/mouse_mil_key.txt.gz", col_names = F, col_types = "ffcc") %>%
  rename(read = 1, tx = 2, seq = 4, cig = 3) %>% 
  mutate(seq_len = nchar(seq)) %>% 
  select(-seq) 

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
  select(-cig)

############################################################
############################################################
############################################################

# plot seq length vs map length 
ggplot(addcig, aes(x = seq_len, y = map_len)) + 
  geom_bin2d() +  
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Sequence length vs alignment length") + 
  xlab("Sequence length (nt)") + ylab("Alignment M length (nt)") + 
  scale_x_log10() + 
  scale_y_log10()

# check pearson correlation 
cor.test(addcig$seq_len, addcig$map_len, method="pearson")

############################################################
############################################################
############################################################

rm(input)
rm(addcig)

# calculate transcript abundance 
txabund <- addcig %>% 
  group_by(tx) %>% 
  mutate(nread = n()) %>% 
  arrange(desc(nread)) %>% 
  ungroup()

# most abundant read 
# ENSMUST00000082390.1; 35990 reads 

# make a scatterplot of seq len vs map len for most abundant read 
ggplot(txabund %>% select(-read) %>% filter(tx == "ENSMUST00000082390.1"), aes(x = seq_len, y = map_len)) + 
  geom_point() +  
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Sequence length vs alignment length for mt-Rnr2-201") + 
  xlab("Sequence length (nt)") + ylab("Alignment M length (nt)") + 
  scale_x_log10() + 
  scale_y_log10()

# make a scatterplot of seq len vs map len for second most abundant read 
ggplot(txabund %>% select(-read) %>% filter(tx == "ENSMUST00000097014.7"), aes(x = seq_len, y = map_len)) + 
  geom_point() +  
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Sequence length vs alignment length for Tuba1a") + 
  xlab("Sequence length (nt)") + ylab("Alignment M length (nt)") + 
  scale_x_log10() + 
  scale_y_log10()

############################################################
############################################################
############################################################

# calculate decay factor 
decay <- txabund %>% 
  select(-read, -nread) %>% 
  group_by(tx) %>% 
  summarise(df = median(map_len), cov = n())

# plot a histogram of coverage per transcript 
ggplot(decay, aes(x = cov)) +
  stat_ecdf() + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Distribution of counts per transcript from 1m reads") + 
  xlab("Counts per transcript") + ylab("Percentile") + 
  scale_x_log10()

# plot a histogram of decay factor 
ggplot(decay %>% filter(cov > 50), aes(x = df)) +
  geom_histogram(bins = 100) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("D50 distribution for transcripts with >50 reads") + 
  xlab("D50") + ylab("Count") 

# calculate mean and std.dev df at 50 reads 
tmp <- decay %>% filter(cov > 49)
mean(tmp$df)  
sqrt(var(tmp$df))

# plot coverage vs d50
ggplot(decay %>% filter(cov > 49), aes(x = cov, y = df)) +
  geom_point(alpha=0.25) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Total coverage vs d50") + 
  xlab("Transcript coverage") + ylab("Transcript d50")+ 
  scale_x_log10() + scale_y_log10()

# check correlation 
cor.test(tmp$cov, tmp$df)

############################################################
############################################################
############################################################

# import transcript attributes 
library(GenomicFeatures)
library(rtracklayer)

# process the annotation 
anno <- "/Users/asethi/localGadiData/2022-03-15_develop-nanograd4/Mus_musculus.GRCm39.104.chr.gtf"

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

