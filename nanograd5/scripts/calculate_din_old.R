

library(tidyverse)

# read in nanograd4 summary data
a <- read_tsv(file, col_names = T)

# calculate din
din <- a %>%
  summarise(sum = sum(df_ratio), n = n()) %>%
  mutate(DIN = 30 * sum/n) %>%
  select(DIN)

print(DIN)
