# import tidyverse
library(tidyverse)

# import tidy data 
wt_1 <- read_tsv("/Users/AJlocal/Downloads/BK_LIQA/undegraded_hek293_pass2_primary.txt", col_names = T)
deg_1 <- read_tsv("/Users/AJlocal/Downloads/BK_LIQA/5mM_MgCl_degrdation_pass1_primary.txt", col_names = T)

# like so, import the replicates 
merged <- inner_join(wt_1, deg_1, by = "GeneName", suffix = c("_wt1", "_deg_1"))

# plot isoform abundance 
ggplot(merged, aes(x = ReadPerGene_corrected_wt1, y = ReadPerGene_corrected_deg_1)) + 
  geom_point()
