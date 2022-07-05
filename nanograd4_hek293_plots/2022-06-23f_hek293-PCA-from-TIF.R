#!/usr/bin/env Rscript

# Make PCA of transcript integrity factors 

#################################################

library(tidyverse)
library(ggrepel)
library(ggtext)

data_path <- "~/localGadiData/2022-06-23_hek293-nanograd4-analysis-TIF"

# import data 
(files <- fs::dir_ls(data_path, glob="*.txt"))
dat <- read_tsv(files, id="path") %>% 
  mutate(sample = sub('\\_nanograd4_out.txt$', '', basename(path))) %>% 
  select(-path)

# create study metadata
metadata <- tibble(sample = c("hassan_SRR13261194", "hassan_SRR13261195", "hassan_SRR13261196", "lorenz_SRR9646141", "lorenz_SRR9646142", "sethi_deg_rep1", "sethi_deg_rep2", "sethi_wt_rep1", "sethi_wt_rep2", "soneson_ERR3218376", "soneson_ERR3218377", "soneson_ERR3218378", "soneson_ERR3218379", "soneson_ERR3218380"), 
                study = c("Hassan_2021", "Hassan_2021", "Hassan_2021", "Lorenz_2020", "Lorenz_2020", "This_study", "This_study", "This_study", "This_study", "Soneson_2019", "Soneson_2019", "Soneson_2019", "Soneson_2019", "Soneson_2019"),
                condition = c("NA", "NA", "NA", "NA", "NA", "Degraded", "Degraded", "WT", "WT", "NA", "NA", "NA", "NA", "NA"))

# Join metadata
anno_data <- inner_join(dat, metadata, by = "sample")

# define factor levels for sample, study and condition
anno_data <- anno_data %>% 
  mutate(condition = factor(condition, levels = c("WT", "Degraded", "NA")), 
         study = factor(study, levels = c("Soneson_2019", "Lorenz_2020", "Hassan_2021", "This_study")), 
         sample = factor(sample, levels = c("soneson_ERR3218376", "soneson_ERR3218377", "soneson_ERR3218378", "soneson_ERR3218379", "soneson_ERR3218380", "lorenz_SRR9646141", "lorenz_SRR9646142", "hassan_SRR13261194", "hassan_SRR13261195", "hassan_SRR13261196", "sethi_wt_rep1", "sethi_wt_rep2", "sethi_deg_rep1", "sethi_deg_rep2")))

# deselect "hassan_SRR13261194", "hassan_SRR13261196"
filt_anno_data <- subset(anno_data, !(sample %in% c("hassan_SRR13261194", "hassan_SRR13261196")))

# plot histograms of DF 
ggplot(filt_anno_data, aes(x = df_ratio, fill = study)) + geom_density() + xlim(0,1) + facet_wrap(~sample, nrow = 3) + 
  xlab("Transcript Integrity Factor distribution") + 
  ggtitle("Transcript Integrity Factor distributions in HEK293 poly(A) sequencing")
  

# filter for transcripts observed in this study
tx_this_study <- filt_anno_data %>% filter(study == "This_study")
protein_coding_tx_this_study <- filt_anno_data %>% filter(transcript_id %in% tx_this_study$transcript_id) %>% filter(transcript_biotype == "protein_coding")
ggplot(protein_coding_tx_this_study, aes(x = df_ratio, fill = study)) + geom_histogram() + xlim(0,1) + facet_wrap(~sample, nrow = 3)

# prepare data for PCA 
pca_data_all <- filt_anno_data %>% 
  rename(dfmedian = df_ratio, tx_count = cov) %>% 
  filter(tx_count > 5) %>% 
  select(transcript_id, sample, dfmedian) %>% 
  pivot_wider(names_from = sample, values_from = dfmedian) %>% 
  na.omit()
         
# make correlation plots 
library(corrplot)

pca_data_all %>% select(-transcript_id) %>% cor(., method = "pearson") %>% 
  corrplot(., method = "pie", type = "upper", order = "hclust", tl.col = "black", tl.srt = 45, col.lim = c(0,1), diag = T)

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
  select(-transcript_id) %>% 
  as.matrix()

project.pca <- prcomp(t(pca_data_temp))
summary(project.pca)
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100

# scree plot 
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

# extract data 
data <- as_tibble(project.pca$x) %>% 
  mutate(sample = rownames(project.pca$x)) %>% 
  inner_join(., metadata, by = "sample")
data

dev.new(noRStudioGD = TRUE)

# Fix study names 
key <- tibble(tmp_study = c("Hassan_2021", "Lorenz_2020", "Soneson_2019", "This_study"),
              Study = c("Hassan et al, 2021", "Lorenz et al, 2020", "Soneson et al, 2019", "This study"))


data2 <- inner_join(data %>% rename(tmp_study = study), key, by = "tmp_study")

data2 <- data2 %>% mutate(Study = factor(Study, levels = c("Soneson et al, 2019", "Lorenz et al, 2020", "Hassan et al, 2021", "This study")))
 
dev.new(noRStudioGD = TRUE)

ggplot(data2, aes(x = PC1, y = PC2, color = Study)) + 
  geom_point(size = 3, alpha = 1) + 
  xlab(paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%")) + 
  ylab(paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%")) + 
  theme(text = element_text(size=16)) + 
  ggtitle(str_wrap("PCA of Transcript Integrity Factors", 50)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  geom_text_repel(data = data2 %>% filter(Study == "This study"),
                  inherit.aes = TRUE,
                  aes(label = condition, x = PC1, y = PC2), 
                  box.padding = 1, 
                  size = 6, 
                  show.legend = F)  

##################################################
##################################################
##################################################

# plot TIF for hek293 data 
this_study <- filt_anno_data %>% filter(study == "This_study") %>% 
  mutate(sample = str_remove(sample, 'sethi_')) %>% 
  mutate(condition = factor(condition, levels = c("Degraded", "WT")))

# plot transcript integrity factors 
ggplot(this_study, aes(x = df_ratio, fill = condition)) + 
  geom_density() + 
  facet_wrap(~sample, nrow = 2) + 
  xlim(0,1) + 
  ggtitle("HEK293 Transcript Integrity Factors") + 
  xlab("Transcript Integrity Factor") + 
  ylab("Density") + 
  theme(text = element_text(size=16)) + 
  theme(plot.title = element_text(hjust=0.5))

a <- tibble(sample = c("wt_rep1", "wt_rep2", "deg_rep1", "deg_rep2"), 
            condition = c("wt", "wt", "deg", "deg"),
            DIN = c("11.8", "11.8", "8.24", "8.11")) %>% 
  mutate(condition = factor(condition, levels = c("wt", "deg"))) %>% 
  mutate(sample = factor(sample, levels = c("wt_rep1", "wt_rep2", "deg_rep1", "deg_rep2"))) %>% 
  mutate(DIN = as.numeric(DIN))

ggplot(a, aes(x = sample, y = DIN, fill= condition)) + geom_bar(stat='identity') + 
  ggtitle("HEK293 Direct Integrity Numbers") + 
  xlab("Sample") + 
  ylab("Direct Integrity Number") + 
  theme(text = element_text(size=24)) + 
  theme(plot.title = element_text(hjust=0.5))            
  

##############################

a <- filt_anno_data %>% filter(transcript_id == "ENST00000376588") %>% arrange(desc(sample)) %>% 
  select(sample, df_ratio)
a

##########

a <- read_tsv("/Users/asethi/Documents/HEK293_DIN.txt", col_names = T, col_types = "fcd") %>% 
  mutate_all(~replace(., is.na(.), 10.7)) %>% 
  mutate(sample = "sample")

ggplot(a, aes(y = Study, x = DIN, color = Study)) + geom_point(size = 3, alpha = 1) + 
  ggtitle("HEK293 Direct Integrity Number distribution") + 
  xlab("Sample DIN") + 
  ylab("Study") + 
  theme(text = element_text(size=16)) + 
  theme(plot.title = element_text(hjust=0.5))  

  

    