library(tidyverse)

# read in nanopolish eventalign file 
a <- read_tsv("/Users/asethi/Documents/asethi-jcsmr/toLocal/nanopolish_test.txt", col_names = T) %>% 
  mutate(current = samples)

# plot event level means only 
ggplot(a, aes(x = position, y = event_level_mean)) + 
  geom_point() + 
  theme(text=element_text(size=24)) + 
  theme(plot.margin = margin(1,1,1,1, "cm"))

# split data so that each signal (col 14) has a new line 
b <- separate_rows(a, current, convert = TRUE) 
b <- b %>% mutate(id = seq(1:nrow(b))) %>%  # add time data 
  mutate(time = id / 3200) %>% # convert each observation to a time value by dividing by 3200 
  mutate(event_status = case_when(
    model_kmer == "NNNNN" ~ "Missed event", # mutate to add 'missed_event' label to certain kmers 
    model_kmer != "NNNNN" ~ "Aligned event")) # otherwise, label 'aligned event'

ggplot(b %>% arrange(event_status), aes(x = time, y = current)) + 
  geom_point(aes(color = event_status, size = event_status, alpha = event_status)) + 
  scale_colour_manual(values = c("azure4", "red")) +
  scale_size_manual(values = c("Missed event" = 1, "Aligned event"=0.6)) + 
  scale_alpha_manual(values = c("Missed event" = 0.6, "Aligned event"=0.5)) + 
  theme(text=element_text(size=18)) + 
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")) + 
  geom_line(aes(x = time, y = current), alpha = 0.2) + 
  xlab("time (seconds)") + 
  ylab("current (mA)") + 
  theme(legend.position = "none") + 
  xlim(50,60)
  
# bash command to find N-rich sequences from nanopolish eventalign output 
# time head -n 50000000 *txt | grep "NNNNN" | cut -f4 | tail -n +2 | uniq -c | sort -k1 -nr | head -n 20
# cat <(cat *txt | head -n 1) <(head -n 5000000 *txt | grep "92b303ed-bdcf-41e1-913c-5b7f99248a24") > ~/toLocal/nanopolish_test.txt


