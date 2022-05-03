library(tidyverse)

# read in squigglepull output 

myskip <- 99

# normal time 
a <- read_tsv("/Users/asethi/localGadiData/2022-05-03_nanograd4-squigglekit-test/read_1.tsv", 
              col_names = F, skip = myskip, n_max = 1) %>%
  head(n = 1) %>% 
  rename(library = X1, read = X2) %>% 
  select(-library) %>% 
  pivot_longer(cols = starts_with("X"), names_to = "time2", 
               values_to = "current") %>% 
  mutate(time = row_number()) %>% 
  arrange(desc(row_number())) %>% 
  mutate(reverse_time = row_number()) %>% 
  select(read, reverse_time, time, current) %>% 
  mutate(time = time / 3200, reverse_time = reverse_time / 3200) 

ggplot(a, aes(x = time, y = current)) + 
  geom_point(size = 0.2) + 
  geom_line(aes(x = time, y = current), alpha = 0.2) 

# reverse time
read_tsv("/Users/asethi/localGadiData/2022-05-03_nanograd4-squigglekit-test/read_1.tsv", 
         col_names = F, skip = myskip, n_max = 1) %>%
  head(n = 1) %>% 
  rename(library = X1, read = X2) %>% 
  select(-library) %>% 
  pivot_longer(cols = starts_with("X"), names_to = "time2", 
               values_to = "current") %>% 
  mutate(time = row_number()) %>% 
  arrange(desc(row_number())) %>% 
  mutate(reverse_time = row_number()) %>% 
  select(read, reverse_time, time, current) %>% 
  mutate(time = time / 3200, reverse_time = reverse_time / 3200) %>% 
  ggplot(., aes(x = reverse_time, y = current)) + 
  geom_point(size = 0.2) + 
  geom_line(aes(x = reverse_time, y = current), alpha = 0.2) +
  xlim(0,0.5)
