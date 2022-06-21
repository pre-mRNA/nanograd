#!/bin/user/env/R

# analyse liqa output from degradation experiment
# written by BK and AS on 10-05-2022


# download data from gadi using a shell command
# scp bk9031@gadi.nci.org.au:/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/LIQA/isoform_expression/* /home/bhavika_kumar/localGadiData/isoform_expression/

###########################

# load tidyverse
library(tidyverse)

# Importing all 4 datasets
wt_1 <- read_tsv("/Users/asethi/localGadiData/2022-05-12_liqa/undegraded_hek293_pass1_primary.txt", col_names=T)
wt_2 <- read_tsv("/Users/asethi/localGadiData/2022-05-12_liqa/undegraded_hek293_pass2_primary.txt", col_names=T)
deg_1 <- read_tsv("/Users/asethi/localGadiData/2022-05-12_liqa/5mM_MgCl_degrdation_pass1_primary.txt", col_names=T)
deg_2 <- read_tsv("/Users/asethi/localGadiData/2022-05-12_liqa/5mM_MgCl_degrdation_pass2_primary.txt", col_names=T)

# Full join wt_1 and wt_2
sum_wt <- full_join(wt_1,wt_2, suffix = c("_wt1", "_wt2"), by = "IsoformName") %>%
  dplyr::mutate(transcript = IsoformName,
                wt_1_count = ReadPerGene_corrected_wt1,
                wt_2_count = ReadPerGene_corrected_wt2) %>%
  mutate(wt_sum = wt_1_count + wt_2_count) %>%
  select(transcript, wt_sum) %>%
  mutate_at(vars(wt_sum), ~ as.integer(round(.x)))

sum_wt

# For deg
sum_deg <- full_join(deg_1,deg_2, suffix = c("_deg1", "_deg2"), by = "IsoformName") %>%
  dplyr::mutate(transcript = IsoformName,
                deg_1_count = ReadPerGene_corrected_deg1,
                deg_2_count = ReadPerGene_corrected_deg2) %>%
  mutate(deg_sum = deg_1_count + deg_2_count) %>%
  select(transcript, deg_sum) %>%
  mutate_at(vars(deg_sum), ~ as.integer(round(.x))) 

sum_deg

# merge sum deg and sum wt
all_sum <- full_join(sum_deg, sum_wt, by = "transcript")
all_sum

# check correlation
cor.test(all_sum$deg_sum, all_sum$wt_sum, method = "spearman")

### plot
g <- ggplot(all_sum, aes(x = wt_sum, y = deg_sum)) +
  geom_point() +
  theme_light() +
  ggtitle(str_wrap("Normalised HEK293 transcript counts, WT vs degraded", 40)) +
  xlab("Normalised trancript count in WT RNA") +
  ylab("Normalised trancript count in degraded RNA") +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title= element_text(size = 16, hjust = 0.5),
        legend.title = element_text(colour="black", size=0, face="bold"),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.position=c(2,2),
        legend.background=element_blank()) +
  scale_colour_manual(values = c("#7474E8", "#9E3954", "#ADA81D")) + # correspond to the level of factor
  scale_x_log10() + scale_y_log10() + 

# plot tpm vs tpm

sum_wt <- sum(all_sum$wt_sum %>% na.omit())
sum_deg <- sum(all_sum$deg_sum %>% na.omit())

tpm <- all_sum %>%
  replace(is.na(.), 0) %>%
  mutate(wt_sum = 1000000 * wt_sum / sum_wt,
         deg_sum = 1000000 * deg_sum / sum_deg)

cor.test(tpm$deg_sum, tpm$wt_sum)

g <- ggplot(tpm, aes(x = wt_sum, y = deg_sum)) +
  geom_point(color="cornflowerblue", size = 0.5) +
  geom_smooth(method = "lm") + 
  theme_light() +
  ggtitle(str_wrap("Normalised HEK293 transcript counts, WT vs degraded", 20)) +
  xlab("Normalised trancript count in WT RNA") +
  ylab("Normalised trancript count in degraded RNA") +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title= element_text(size = 16, hjust = 0.5),
        legend.title = element_text(colour="black", size=0, face="bold"),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.position=c(2,2),
        legend.background=element_blank()) +
  scale_colour_manual(values = c("#7474E8", "#9E3954", "#ADA81D"))

# saving the plot
ggsave("~/Desktop/tpm.png", plot=g, scale=1, width=4.5, height=4.5, unit="in", dpi=300, bg=NULL)


# plot tpm vs tpm

sum_wt <- sum(all_sum$wt_sum %>% na.omit())
sum_deg <- sum(all_sum$deg_sum %>% na.omit())

tpm <- all_sum %>%
  replace(is.na(.), 0) %>%
  mutate(wt_sum = 1000000 * wt_sum / sum_wt,
         deg_sum = 1000000 * deg_sum / sum_deg)

ggplot(tpm, aes(x = wt_sum, y = deg_sum)) +
  geom_point() +
  theme_light() +
  ggtitle(str_wrap("Normalised HEK293 transcript counts, WT vs degraded", 40)) +
  xlab("Normalised trancript count in WT RNA") +
  ylab("Normalised trancript count in degraded RNA") +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title= element_text(size = 16, hjust = 0.5),
        legend.title = element_text(colour="black", size=0, face="bold"),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.position=c(2,2),
        legend.background=element_blank()) +
  scale_colour_manual(values = c("#7474E8", "#9E3954", "#ADA81D"))
