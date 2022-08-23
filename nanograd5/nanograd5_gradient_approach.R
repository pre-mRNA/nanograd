# written by AJ Setih on 2022-08-02
# compare the gradient approach to transcript integrity factor 
# import the nanograd5 temporary output

# setup
library(tidyverse)
library(corrplot)
setwd("/Users/asethi/localGadiData/2022-08-04_test-nanograd5-decayFirst4")

# write a function to import the nanograd5 tmp output 
readTemp <- function(path, mySample, myCondition){
  n5temp <- read_tsv(path, col_names = T, col_types = "fddfffddddd") %>% 
    mutate(sample = mySample, condition = myCondition)
}

# read in data 
allRows <- bind_rows(readTemp("all.undegraded_hek293_pass1.fastq.gz.sorted_nanograd5_out.txt.tmp", "wt_rep1", "wt"),
          readTemp("all.undegraded_hek293_pass2.fastq.gz.sorted_nanograd5_out.txt.tmp", "wt_rep2", "wt"),
          readTemp("all.5mM_MgCl_degrdation_pass1.fastq.gz.sorted_nanograd5_out.txt.tmp", "deg_rep1", "deg"),
          readTemp("all.5mM_MgCl_degrdation_pass2.fastq.gz.sorted_nanograd5_out.txt.tmp", "deg_rep2", "deg")) %>% 
  mutate(condition = factor(condition, levels = c("wt", "deg")))

# note that df_ratio is equal to TIN 
# comparing df_ratio between conditions 
compare_tin <- allRows %>%
  filter(cov > 20) %>% # 53.6k transcripts without omitting NAs and without a coverage threshold, 6.3k transcripts with a coverage thredhold
  filter(transcript_biotype == "protein_coding") %>% 
  select(transcript_id, df_ratio, sample) %>% 
  pivot_wider(names_from="sample", values_from="df_ratio") %>% 
  na.omit() # 6.3k transcripts with coverage threshold to 682 transcripts without 

# make correlation plot 
compare_tin %>% 
  dplyr::select(c(2:5)) %>% 
  cor(., method = "pearson") %>% 
  corrplot(., 
           method = "pie", 
           type = "upper", 
           tl.col = "black", 
           tl.srt = 45, 
           diag = T, 
           order = c("hclust"),
           hclust.method = c("centroid"))

# write code to plot TIN for a specific transcript 
plotTin <- function(transcript){ 
  allRows %>%
  filter(cov > 20) %>% # 53.6k transcripts without omitting NAs and without a coverage threshold, 6.3k transcripts with a coverage thredhold
  filter(transcript_biotype == "protein_coding") %>% 
  select(transcript_id, df_ratio, sample, condition, gene_name) %>% 
  filter(transcript_id == transcript) %>% 
  arrange(., by = condition) %>% 
  ggplot(., aes(x = sample, y = df_ratio, fill = condition)) + 
    geom_bar(stat='identity')
}

plotTin("ENST00000650309")

# t-test between transcripts 
# t-test between col 2/3 and 4/5; rowwise t-test 
