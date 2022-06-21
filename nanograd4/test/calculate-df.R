

d1 <- read_tsv("/Users/asethi/localGadiData/2022-04-05_test-nanograd4-decay/all.5mM_MgCl_degrdation_pass1.fastq.gz.sorted_nanograd4_out.txt") %>% mutate(sample = "deg_1")
d2 <- read_tsv("/Users/asethi/localGadiData/2022-04-05_test-nanograd4-decay/all.5mM_MgCl_degrdation_pass2.fastq.gz.sorted_nanograd4_out.txt") %>% mutate(sample = "deg_2")
w1 <- read_tsv("/Users/asethi/localGadiData/2022-04-05_test-nanograd4-decay/all.undegraded_hek293_pass1.fastq.gz.sorted_nanograd4_out.txt") %>% mutate(sample = "wt_1")
w2 <- read_tsv("/Users/asethi/localGadiData/2022-04-05_test-nanograd4-decay/all.undegraded_hek293_pass2.fastq.gz.sorted_nanograd4_out.txt") %>% mutate(sample = "wt_2")

merge <- bind_rows(d1, d2, w1, w2)

# sample
merge %>% group_by(sample) %>% 
  summarise(sum = sum(df_ratio), n = n()) %>% 
  mutate(DIN = 30 * sum/n)
