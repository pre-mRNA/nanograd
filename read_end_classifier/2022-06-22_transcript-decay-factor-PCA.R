#!/usr/bin/env Rscript

# analyse merged sequencing summary and degradation data 
# written by AS on 2022-06-21

#################################################

library(tidyverse)

data_path <- "/Users/AJlocal/localGadiData/2022-06-21_degradation-all-data/2022-06-21_all_degradation_combined.txt.gz"

# import data 
input <- read_tsv(data_path, col_names = T, col_types = "ffddfddfffddf")

# prepare data for PCA 
pca_data_all <- input %>% 
  select(condition, transcript_id, df) %>% 
  group_by(condition, transcript_id) %>% 
  summarise(dfmedian = median(df), read_count = n()) %>% 
  group_by(transcript_id) %>% 
  mutate(tx_count = sum(read_count)) %>% 
  ungroup() %>% 
  select(-read_count) %>% 
  pivot_wider(names_from = condition, values_from = dfmedian) %>% 
  replace(is.na(.), 0)

rm(input)

# plot histogram of read_counts 
ggplot(pca_data_all, aes(x = tx_count)) + 
  stat_ecdf() + 
  scale_y_log10() + 
  scale_x_log10()

############################################################
############################################################

library(ggfortify)
library(ggrepel)

pca_data_temp <- pca_data_all %>% 
  filter(tx_count > 20) %>% 
  select(transcript_id, wt_rep1, wt_rep2, deg_rep1, deg_rep2) %>% 
  select(-transcript_id) %>% 
  as.matrix()

project.pca <- prcomp(t(pca_data_temp))
summary(project.pca)
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100

# scree plot 
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

# extract data 
data <- as_tibble(project.pca$x) %>% 
  mutate(Condition = c("Control", "Control", "Degraded", "Degraded")) %>% 
  mutate(sample = c("Control_rep1", "Control_rep2", "Deg_rep1", "Deg_rep2"))
data

dev.new(noRStudioGD = TRUE)
ggplot(data, aes(x = PC1, y = PC2, color = Condition)) + 
  geom_point(size = 3, alpha = 1) + 
  xlab(paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%")) + 
  ylab(paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%")) + 
  theme(text = element_text(size=24)) + 
  ggtitle(str_wrap("PCA of nanograd transcript decay factors", 60)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  geom_text_repel(data = data,
                  inherit.aes = TRUE,
                  aes(label = sample, x = PC1, y = PC2), 
                  box.padding = 1, 
                  size = 6)
                 