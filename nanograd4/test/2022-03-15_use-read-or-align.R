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

