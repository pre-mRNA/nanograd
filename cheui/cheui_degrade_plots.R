# make plots for CHEUI annotate2 output 

library(tidyverse)
library(scales)
library(ggrepel)

####################################################################################################
####################################################################################################
####################################################################################################

# import annotate2 liftover using a function 

# write a function to import each condition's cheui output 
importAnnotateTwo <- function(cFile, sName, sModel) {
  read_tsv(cFile, col_names = T, skip = 0, col_types = "fnncccncnnnffcccnnnnnnndnnnncc") %>% 
    mutate(sample = sName, model = sModel)
}

# import cheui output for all conditions
wt <- importAnnotateTwo("~/localGadiData/2022-05-14_cheui-degradation/nanopolish_HEK293_undegraded_pass2_liftover.txt", "WT", "m6a")
deg <- importAnnotateTwo("~/localGadiData/2022-05-14_cheui-degradation/nanopolish_5mM_MgCl_degradation_pass2_liftover.txt", "deg", "m6a")

# merge the rows 
compiledRows <- bind_rows(wt, deg) %>% 
  mutate(sample = factor(sample, levels = c("WT", "deg"))) %>% 
  mutate(model = factor(model, levels = c("m6a")))

# add a filter for significant sites, and count how many times a site is tested and significnat 
comp <- compiledRows %>% 
  mutate(filter = ifelse(prob > 0.9999, "sig", "ns")) %>% 
  group_by(transcript, start, strand, model) %>% 
  mutate(n_tested = n()) %>% 
  mutate(n_sig = sum(prob > 0.9999)) %>% 
  ungroup() 

# test to ensure that no site is captured more than 5 times 
a <- comp %>% slice_sample(n = 100000)

ggplot(a, aes(x = n_tested, y = n_sig)) + 
  geom_point()

# remove the old variables from memory 
rm(list = ls()[grep("^compiled", ls())])

# save the output 
# write_tsv(comp, "~/localGadiData/2022-01-31_run-annotate2-all-pval/comp.txt.gz", col_names = T)

# read the output back in 
# comp <- read_tsv("~/localGadiData/2022-01-31_run-annotate2-all-pval/comp.txt.gz", col_names = T)

####################################################################################################
####################################################################################################

# process the annotation metadata for downstream filtering 

library(GenomicFeatures)
library(rtracklayer)

annotation <- "/Users/asethi/localGadiData/2022-05-11_nanograd-metaplots-R/Homo_sapiens.GRCh38.105.chr.gtf"

# process the annotation 

# read in reference transcripts
gtf <- makeTxDbFromGFF(file=annotation, format = "gtf")

# fetch the lengthe of each transcript segment from the gtf 
txlen <- transcriptLengths(gtf, with.cds_len=TRUE, with.utr5_len=TRUE, with.utr3_len=TRUE) %>% 
  as_tibble() %>% 
  mutate(diff = tx_len - cds_len - utr5_len - utr3_len) %>% 
  dplyr::rename(transcript_id = tx_name) %>% 
  dplyr::select(-tx_id, -nexon, -gene_id, -diff)

# the last command doesn't store biotype, so we read in the gtf again using another package 
tx_biotype <- rtracklayer::import(annotation) %>% 
  as_tibble() %>% 
  dplyr::select(transcript_id, transcript_biotype, gene_name, gene_id) %>% 
  na.omit() %>% 
  distinct()

# merge the biotypes with the transcript segment lengths 
merged_metadata <- inner_join(tx_biotype, txlen, by = "transcript_id")

####################################################################################################
####################################################################################################

# make metaplots for all genes (testing only, not in powerpoint )

# straight up histogram for all tested sits 
ggplot(comp, aes(x = rel_pos, fill = model)) +
  geom_histogram(binwidth = 0.01) + 
  geom_vline(xintercept = 1:2, col = "grey") + 
  theme_bw() + 
  theme(text = element_text(size=24)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Metaplot for all tested sites") + 
  xlab("Relative location") + ylab("Count") + 
  scale_y_log10() + 
  facet_grid(~model) 

# normalized abundance plot
breaks <- seq(0,3,0.01)

out_ratio <- comp %>% 
  mutate(interval = cut(rel_pos, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)) %>%
  group_by(model, interval, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))

ggplot(out_ratio, aes(x = interval, y = ratio, color = model)) + 
  geom_point(alpha = 0.2) + geom_smooth() +  
  geom_vline(xintercept = c(100,200), col = "black") + 
  theme_bw() + 
  theme(text = element_text(size=24)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Ratio of significant sites") + 
  xlab("Relative location") + ylab("Count") + 
  scale_y_log10() + 
  facet_grid(~model) 

####################################################################################################
####################################################################################################

# plot the redundancy of tested sites 

red_data <- comp %>% group_by(chr, start, sample, model) %>% 
  summarise(n = n())

ggplot(red_data, aes(x = n, fill = interaction(model,sample))) + 
  geom_histogram(alpha = 0.8, binwidth = 1) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Redundancy of tested sites") + 
  xlab("Transcriptomic sites per genomic site") + ylab("Count") + 
  scale_y_log10() + 
  facet_grid(sample ~ model) + 
  xlim(0,8) # add limit for mutant E18, which goes up to 12 sites (while other conditions are limited at 8)

# sepearate by metagene location 
red_data_metagene <- comp %>% 
  rowwise() %>% 
  mutate(metagene = case_when(rel_pos < 1 ~ "5UTR",
                              (rel_pos > 1 & rel_pos < 2) ~ "CDS",
                              rel_pos > 2 ~ "3UTR")) %>% 
  group_by(chr, start, metagene, sample, model) %>% 
  summarise(n = n())

# omit NA observations
# takes us from 18.29m observations down to 18.28 M
red_data_metagene <- red_data_metagene %>% 
  na.omit() %>% 
  mutate(metagene= factor(metagene, levels = c("5UTR", "CDS", "3UTR")))

ggplot(red_data_metagene, aes(x = n, fill = metagene)) + 
  geom_histogram(alpha = 1, binwidth = 1, position = 'identity') + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Redundancy of tested sites") + 
  xlab("Transcriptomic sites per genomic site") + ylab("Count") + 
  scale_y_log10() + 
  facet_grid(sample ~ model ~ metagene) + 
  xlim(0,8) # add limit for mutant E18, which goes up to 12 sites (while other conditions are limited at 8)


# remove the old variables from memory 
rm(list = ls()[grep("^red_data", ls())])

####################################################################################################
####################################################################################################

# plot the number and proportion of significant sites by metagene location and by sample 

Sys.time()

# start by adding metagene location to the compiled data 
# then summarise how many sites are assigned to each condition 
comp_meta <- comp %>% 
  mutate(metagene = case_when(rel_pos < 1 ~ "5UTR",
                              (rel_pos > 1 & rel_pos < 2) ~ "CDS",
                              rel_pos > 2 ~ "3UTR")) %>% 
  mutate(metagene = factor(metagene, levels = c("5UTR", "CDS", "3UTR"))) %>% 
  group_by(sample, model, metagene) %>% 
  summarise(n_site = n(),
            n_sig_tx_site = sum(prob > 0.9999),
            percent_sig_tx_site = n_sig_tx_site / n_site) 

Sys.time()

# fucnction for plotting with scientific labels
# taken from: https://stackoverflow.com/a/18530540/6889619 
scientific_10 <- function(x) {
  gsub("e", " x 10^", scientific_format()(x))
}

# plot number of tested sites per region 
ggplot(comp_meta %>% na.omit(), aes(x = metagene, y = n_site, fill = metagene)) + 
  geom_bar(alpha = 1, stat = "identity") + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Number of tested sites in tx space") + 
  xlab("Metatranscriptomic region") + ylab("Number of tested sites") + 
  scale_y_continuous(label=scientific_10) + 
  facet_grid(sample ~ model) # + xlim(0,8) # add limit for mutant E18, which goes up to 12 sites (while other conditions are limited at 8)

# make the plot 
ggplot(comp_meta %>% na.omit(), aes(x = metagene, y = percent_sig_tx_site, fill = metagene)) + 
  geom_bar(alpha = 1, stat = "identity") + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Percentage significant sites in tx space") + 
  xlab("Metatranscriptomic region") + ylab("Proportion of significant sites") + 
  facet_grid(sample ~ model) # + xlim(0,8) # add limit for mutant E18, which goes up to 12 sites (while other conditions are limited at 8)

# plot the relationship between tested and signifiant sites 
ggplot(comp_meta %>% na.omit(), aes(x = n_site, y = percent_sig_tx_site, color = metagene, shape = sample)) + 
  geom_point(size = 10) + 
  theme_bw() + 
  theme(text = element_text(size=24)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Percentage significant sites in tx space") + 
  xlab("Number of tested sites") + ylab("Proportion of significant sites") +
  scale_x_log10() + scale_y_log10() 

comp_meta2 <- comp_meta %>% na.omit()
cor.test(comp_meta2$n_site, comp_meta2$percent_sig_tx_site)


####################################################################################################
####################################################################################################

# make metaplots for all samples, actual version 

# straight up histogram for all tested sits 
ggplot(comp, aes(x = rel_pos, fill = model)) +
  geom_histogram(binwidth = 0.01) + 
  geom_vline(xintercept = 1:2, col = "grey") + 
  theme_bw() + 
  theme(text = element_text(size=24)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Metatranscript density of all tested sites") + 
  xlab("Metatranscript location") + ylab("Number of tested sites") + 
  facet_grid(sample ~ model) 


# straight up histogram for all significant sites 
ggplot(comp %>% filter(prob > 0.9999), aes(x = rel_pos, fill = model)) +
  geom_histogram(binwidth = 0.1) + 
  geom_vline(xintercept = 1:2, col = "grey") + 
  theme_bw() + 
  theme(text = element_text(size=24)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Metatranscript density of all significant sites") + 
  xlab("Metatranscript location") + ylab("Number of significant sites") + 
  facet_grid(sample ~ model) 

# normalized abundance plot
breaks <- seq(0,3,0.025) # iterate over break width 

Sys.time()
out_ratio <- comp %>% 
  mutate(interval = cut(rel_pos, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)) %>%
  group_by(sample, model, interval, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))
Sys.time()

ggplot(out_ratio, aes(x = interval, y = ratio, color = model)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(span = 0.2) + # iterate over loess span 
  geom_vline(xintercept = c(80,40), col = "black") + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Ratio of significant sites in all transcripts") + 
  xlab("Relative metatranscriptomic location") + ylab("Proportion of significant sites") + 
  facet_grid(sample ~ model) # +
ylim(0,0.02)

# check how many transcripts there are 
comp %>% dplyr::select(transcript) %>% unique()

####################################################################################################
####################################################################################################

# make metaplots for single-isoform protein-coding transcripts 

# make a list of genes with single biotypes 
single_biotype_gene_keys <- tx_biotype %>% 
  dplyr::select(gene_id, transcript_biotype) %>% 
  group_by(gene_id) %>% 
  summarise(n = n()) %>% 
  filter(n == 1) %>% 
  dplyr::select(gene_id)
single_biotype_gene_keys %>% nrow() # 32886

# make a list of genes with only protein coding biotypes 
protein_coding_gene_keys <- tx_biotype %>% 
  filter(transcript_biotype == "protein_coding") %>% 
  dplyr::select(gene_id) %>% 
  distinct()
protein_coding_gene_keys %>% nrow() # 21761

# make a list of genes with a single transcript 
single_tx_keys <- tx_biotype %>% 
  dplyr::select(transcript_id, gene_id) %>% 
  group_by(gene_id) %>% 
  summarise(n = n()) %>% 
  filter(n == 1)

# now find the intersection of this list 
# genes with one transcript only that is protein-coding 
target_gene_id <- single_tx_keys %>% filter(gene_id %in% single_biotype_gene_keys$gene_id) %>% filter(gene_id %in% protein_coding_gene_keys$gene_id)

target_gene_id %>% nrow() # genes meet our criteria 
# 4260 genes 

# calculate the proportion of significant sites in each bin 
Sys.time()
out_ratio <- comp %>% 
  filter(gene_id %in% target_gene_id$gene_id) %>% 
  mutate(interval = cut(rel_pos, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)) %>%
  group_by(sample, model, interval, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))
Sys.time()

# plot 
ggplot(out_ratio, aes(x = interval, y = ratio, color = model)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(span = 0.2) + # iterate over loess span 
  geom_vline(xintercept = c(80,40), col = "black") + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Ratio of significant sites in single-isoform protein-coding genes") + 
  xlab("Relative metatranscriptomic location") + ylab("Proportion of significant sites") + 
  facet_grid(sample ~ model) + 
  ylim(0,0.01)

##################################################

# filter by heuristic 

# calculate peak metagene coverage and filter for those transcripts where the peak coverage is at 2.9 or higher 
hstic_tx <- comp %>%
  dplyr::select(transcript, transcript_biotype, coverage, rel_pos) %>% 
  filter(transcript_biotype == "protein_coding") %>% 
  group_by(transcript) %>%
  dplyr::slice(which.max(coverage)) %>% 
  filter(rel_pos > 2.7 ) %>% 
  dplyr::select(transcript)

hstic_tx %>% dim() # 4870 transcripts 

# normalized abundance plot
breaks <- seq(0,3,0.025) # iterate over break width 

Sys.time()
out_ratio <- comp %>% 
  filter(transcript %in% hstic_tx$transcript) %>% 
  mutate(interval = cut(rel_pos, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)) %>%
  group_by(sample, model, interval, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))
Sys.time()  

# plot 
ggplot(out_ratio, aes(x = interval, y = ratio, color = model)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(span = 0.2) + # iterate over loess span 
  geom_vline(xintercept = c(80,40), col = "black") + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Ratio of significant sites in heuristic-selected protein-coding genes") + 
  xlab("Relative metatranscriptomic location") + ylab("Proportion of significant sites") + 
  facet_grid(sample ~ model) + 
  ylim(0,0.01)

# add protein-coding genes with single transcripts 
pre_join <- left_join(pre, target_gene_id %>% mutate(single_transcript = "single_transcript"), by = "gene_id")
pre_join$single_transcript <- pre_join$single_transcript %>% replace_na("multi_transcript")

ggplot(pre_join %>% dplyr::rename(gene_type = single_transcript), aes(x = rel_pos, fill = gene_type)) + 
  geom_histogram(alpha = 0.8, bins = 40) + 
  theme_bw() + 
  theme(text = element_text(size=16)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Peak metagenic coverage, 1 site per gene") + 
  xlab("Metagene location") + ylab("Number of genes") + 
  facet_grid(~gene_type)

pre_join %>% filter(single_transcript == "single_transcript") %>% arrange(desc(coverage)) %>% select(gene_id, gene_name) %>% distinct()
meta %>% filter(gene_id == "ENSMUSG00000072235") %>% ggplot(aes(x = -rel_pos, y = coverage)) + geom_point()
meta %>% filter(gene_id == "ENSMUSG00000072235") %>% ggplot(aes(x = -abs_cds_start, y = coverage)) + geom_point()


####################################################################################################
####################################################################################################

# make metaplots for consensus-only sites 
cons_sites <- comp %>%
  mutate(filter = case_when(n_sig <= 1 ~ "ns",
                            n_sig > 1 ~ "sig"))

# normalized abundance plot
breaks <- seq(0,3,0.025) # iterate over break width 

Sys.time()
out_ratio <- comp %>% 
  mutate(interval = cut(rel_pos, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)) %>%
  group_by(sample, interval, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))
Sys.time()  

# plot 
g <- ggplot(out_ratio, aes(x = interval, y = ratio, color = sample)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(span = 0.5) + # iterate over loess span 
  geom_vline(xintercept = c(80,40), col = "black") + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  xlab("Relative metatranscriptomic location") + ylab("Proportion of m6a bases") + 
  ylim(0,0.01)+ 
  theme_light() +
  ggtitle(str_wrap("Discovery of transcriptome-wide m6a sites using CHEUI", 80)) +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title= element_text(size = 16, hjust = 0.5),
        legend.title = element_text(colour="black", size=0, face="bold"),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.15, 'cm'),
        legend.background=element_blank()) +
  scale_colour_manual(values = c("#7474E8", "#9E3954", "#ADA81D"))

# saving the plot
ggsave("~/Desktop/c.png", plot=g, scale=1, width=8, height=4.5, unit="in", dpi=300, bg=NULL)

####################################################################################################
####################################################################################################

# plot around splice juntions 
sj_data <- comp %>% 
  mutate(up_junc_dist = -up_junc_dist) %>% 
  pivot_longer(cols = c(up_junc_dist, down_junc_dist), names_to = "type_dist", values_to = "junc_dist") %>% 
  dplyr::select(-type_dist) %>% 
  dplyr::filter(junc_dist < 55 & junc_dist > -55) %>% 
  group_by(sample, model, junc_dist, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))

# plot 
ggplot(sj_data, aes(x = junc_dist, y = ratio, color = model)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(span = 0.2) + # iterate over loess span 
  geom_vline(xintercept = 0, col = "black") + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Methylation around splice junctions in all samples") + 
  xlab("Absolute distance to splice junction (NT)") + ylab("Proportion of significant sites") + 
  facet_grid(sample ~ model) 


##################################################


# plot around splice juntions, but increasing view from -150 to 150  
sj_data <- comp %>% 
  mutate(up_junc_dist = -up_junc_dist) %>% 
  pivot_longer(cols = c(up_junc_dist, down_junc_dist), names_to = "type_dist", values_to = "junc_dist") %>% 
  dplyr::select(-type_dist) %>% 
  dplyr::filter(junc_dist < 355 & junc_dist > -355) %>% 
  group_by(sample, model, junc_dist, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))

# plot 
ggplot(sj_data, aes(x = junc_dist, y = ratio, color = model)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(span = 0.2) + # iterate over loess span 
  geom_vline(xintercept = 0, col = "black") + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Methylation around splice junctions in all samples") + 
  xlab("Absolute distance to splice junction (NT)") + ylab("Proportion of significant sites") + 
  facet_grid(sample ~ model) 


##################################################

# plot around splice juntions but merge data by model 
sj_data <- comp %>% 
  mutate(up_junc_dist = -up_junc_dist) %>% 
  pivot_longer(cols = c(up_junc_dist, down_junc_dist), names_to = "type_dist", values_to = "junc_dist") %>% 
  dplyr::select(-type_dist) %>% 
  dplyr::filter(junc_dist < 55 & junc_dist > -55) %>% 
  group_by(model, junc_dist, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))

# plot 
ggplot(sj_data, aes(x = junc_dist, y = ratio, color = model)) + 
  geom_point(alpha = 0.8) + 
  geom_smooth(span = 0.2) + # iterate over loess span 
  geom_vline(xintercept = 0, col = "black") + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Methylation around splice junctions, merging samples") + 
  xlab("Absolute distance to splice junction (NT)") + ylab("Proportion of significant sites") + 
  facet_grid(~ model) 

##################################################

# check consensus sites 
cons_sj_sites <- comp %>% 
  mutate(up_junc_dist = -up_junc_dist) %>% 
  pivot_longer(cols = c(up_junc_dist, down_junc_dist), names_to = "type_dist", values_to = "junc_dist") %>% 
  dplyr::select(-type_dist) %>% 
  dplyr::filter(junc_dist < 55 & junc_dist > -55) %>% 
  mutate(filter = case_when(n_sig <= 1 ~ "ns",
                            n_sig > 1 ~ "sig"))


out_ratio <- cons_sj_sites %>% 
  group_by(model, junc_dist, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))


# plot 
ggplot(out_ratio, aes(x = junc_dist, y = ratio, color = model, alpha = ratio)) + 
  geom_point(alpha = 1) + 
  geom_smooth(span = 0.2) + # iterate over loess span 
  geom_vline(xintercept = c(0,0), col = "black") + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Consensus methylation around splice junctions (n >= 2) across all genes") + 
  xlab("Absolute distance from splice junction (NT)") + ylab("Proportion of significant sites") + 
  facet_grid(~ model) 


#########

# check which genes have consensus m5c at the +5 position 
targs <- cons_sj_sites %>% 
  filter(junc_dist == 5) %>% 
  filter(model == "m5c") %>% 
  filter(filter == "sig") %>% 
  dplyr::select(gene_name) %>% 
  unique()

clipr::write_clip(targs)

##################################################

# check consensus sites in our single_tranascript protein_coding genes 
cons_sj_sites <- comp %>% 
  filter(transcript %in% target_gene_id$transcript) %>% 
  mutate(up_junc_dist = -up_junc_dist) %>% 
  pivot_longer(cols = c(up_junc_dist, down_junc_dist), names_to = "type_dist", values_to = "junc_dist") %>% 
  dplyr::select(-type_dist) %>% 
  dplyr::filter(junc_dist < 55 & junc_dist > -55) %>% 
  mutate(filter = case_when(n_sig <= 1 ~ "ns",
                            n_sig > 1 ~ "sig"))

Sys.time()
out_ratio <- cons_sj_sites %>% 
  mutate(interval = cut(rel_pos, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)) %>%
  group_by(model, junc_dist, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))
Sys.time()  

# plot 
ggplot(out_ratio, aes(x = junc_dist, y = ratio, color = model)) + 
  geom_point(alpha = 0.8) + 
  geom_smooth(span = 0.2) + # iterate over loess span 
  geom_vline(xintercept = c(0,0), col = "black") + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Consensus methylation around splice junctions (n >= 2) across single-isoform protein-coding genes") + 
  xlab("Relative metatranscriptomic location") + ylab("Proportion of significant sites") + 
  facet_grid(~ model) 

##################################################

# check consensus sites in our single_tranascript protein_coding genes ONLY IN CDS 
cons_sj_sites <- comp %>% 
  filter(gene_id %in% target_gene_id$gene_id) %>% 
  mutate(up_junc_dist = -up_junc_dist) %>% 
  pivot_longer(cols = c(up_junc_dist, down_junc_dist), names_to = "type_dist", values_to = "junc_dist") %>% 
  dplyr::select(-type_dist) %>% 
  dplyr::filter(junc_dist < 55 & junc_dist > -55) %>% 
  mutate(filter = case_when(n_sig <= 1 ~ "ns",
                            n_sig > 1 ~ "sig"))

Sys.time()
out_ratio <- cons_sj_sites %>% 
  mutate(interval = cut(rel_pos, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)) %>%
  group_by(model, junc_dist, filter) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "filter", values_from = "n") %>% 
  mutate(ratio = sig / (sig + ns))
Sys.time()  

# plot 
ggplot(out_ratio, aes(x = junc_dist, y = ratio, color = model)) + 
  geom_point(alpha = 0.8) + 
  geom_smooth(span = 0.2) + # iterate over loess span 
  geom_vline(xintercept = c(0,0), col = "black") + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Consensus methylation around splice junctions (n >= 2) across single-isoform protein-coding genes in CDS only") + 
  xlab("Absolute distance from splice junction") + ylab("Proportion of significant sites") + 
  facet_grid(~ model) 


####################################################################################################

# plot mean methylation stoichiometry around splice juction

# unweighted 
sj_stoich <- comp %>% 
  mutate(up_junc_dist = -up_junc_dist) %>% 
  pivot_longer(cols = c(up_junc_dist, down_junc_dist), names_to = "type_dist", values_to = "junc_dist") %>% 
  dplyr::select(-type_dist) %>% 
  dplyr::filter(junc_dist < 55 & junc_dist > -55) %>% 
  group_by(model, junc_dist, filter) %>% 
  summarise(n = n(), mean_stoich = mean(stoich)) 

# plot 
ggplot(sj_stoich, aes(x = junc_dist, y = mean_stoich, color = model)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(span = 0.2) + # iterate over loess span 
  geom_vline(xintercept = 0, col = "black") + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Mean methylation stoichiometry around splice junctions in all samples") + 
  xlab("Absolute distance to splice junction (NT)") + ylab("Mean position stoichiometry") + 
  facet_grid(~ model) + 
  scale_y_log10()

#####################################################

# weight by probability 
sj_stoich <- comp %>% 
  mutate(up_junc_dist = -up_junc_dist) %>% 
  pivot_longer(cols = c(up_junc_dist, down_junc_dist), names_to = "type_dist", values_to = "junc_dist") %>% 
  dplyr::select(-type_dist) %>% 
  dplyr::filter(junc_dist < 55 & junc_dist > -55) %>% 
  group_by(model, junc_dist, filter) %>% 
  summarise(n = n(), mean_stoich = weighted.mean(stoich, log(prob))) 

# plot 
ggplot(sj_stoich, aes(x = junc_dist, y = mean_stoich, color = model)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(span = 0.2) + # iterate over loess span 
  geom_vline(xintercept = 0, col = "black") + # iterate vertical lines to match the breaks 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Log-probability-weighted methylation stoichiometry around splice junctions in all samples") + 
  xlab("Absolute distance to splice junction (NT)") + ylab("Mean log-prob weighted position stoichiometry") + 
  facet_grid(~ model) + 
  scale_y_log10()

####################################################################################################
####################################################################################################

# isoform-specific modification

# checkpoint
library(tidyverse)
library(scales)
library(ggrepel)
comp <- read_tsv("~/localGadiData/2022-01-31_run-annotate2-all-pval/comp.txt.gz", col_names = T)
merged_metadata <- read_tsv("~/localGadiData/2022-01-31_run-annotate2-all-pval/merged_metadata.txt.gz", col_names = T)

##########################################################

# add a score, prob * stoich 
scored <- comp %>% 
  mutate(score = prob * stoich) %>% 
  mutate(metagene = case_when(rel_pos < 1 ~ "5UTR",
                              (rel_pos > 1 & rel_pos < 2) ~ "CDS",
                              rel_pos > 2 ~ "3UTR")) %>% 
  mutate(metagene = factor(metagene, levels = c("5UTR", "CDS", "3UTR"))) %>% 
  group_by(chr, start, model, sample) %>% 
  mutate(delta_score = max(score) - min(score)) %>% 
  ungroup()

# plot the relationship between score and delta score 
ggplot(scored %>% dplyr::select(score, delta_score, model, metagene, sample) %>% na.omit(), aes(x = score, y = delta_score)) + 
  geom_bin2d(bins = 100) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Relationship between score and delta-score for all sites") + 
  xlab("Site score in transcriptome space (prob * stoich)") + ylab("Max delta score at genomic position between isoforms") + 
  facet_grid(sample ~ model) 

# plot the distribution of delta score by metagene location
ggplot(scored %>% dplyr::select(delta_score, model, metagene, sample) %>% na.omit(), aes(x = delta_score, fill = metagene)) + 
  geom_histogram(binwidth = 0.01, alpha = 0.6) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Delta score between isoforms within conditions") + 
  xlab("delta score") + ylab("abundance") + 
  scale_y_log10() + 
  facet_grid(sample ~ model) 

##########################################################


# remake the same plots except for genomic positions considered significant in at least one condition 
sig_genome <- comp %>% 
  dplyr::select(chr, start, sample, model, prob) %>% 
  filter(prob > 0.9999) %>% 
  dplyr::select(-prob) %>% 
  unique()

sig_all <- inner_join(sig_genome, scored, by = c("chr", "start", "sample", "model"))
# 44329 site-level calls 

sig_scored <- sig_all %>% 
  mutate(score = prob * stoich) %>% 
  mutate(metagene = case_when(rel_pos < 1 ~ "5UTR",
                              (rel_pos > 1 & rel_pos < 2) ~ "CDS",
                              rel_pos > 2 ~ "3UTR")) %>% 
  mutate(metagene = factor(metagene, levels = c("5UTR", "CDS", "3UTR"))) %>% 
  group_by(chr, start, sample, model) %>% 
  mutate(delta_score = max(score) - min(score)) %>% 
  ungroup()

# plot the relationship between score and delta score 
ggplot(sig_scored %>% dplyr::select(score, delta_score, model, metagene, sample) %>% na.omit(), aes(x = score, y = delta_score)) + 
  geom_bin2d(bins = 40) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Relationship between score and delta-score for significant genomic sites within conditions") + 
  xlab("Site score in transcriptome space (prob * stoich)") + ylab("Max delta score at genomic position between isoforms") + 
  facet_grid(sample ~ model) 

# plot the distribution of delta score by metagene location
ggplot(sig_scored %>% dplyr::select(delta_score, model, metagene, sample) %>% na.omit(), aes(x = delta_score, fill = metagene)) + 
  geom_histogram(binwidth = 0.03, alpha = 0.6) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Delta score between isoforms across conditions for significant genomic sites") + 
  xlab("delta score") + ylab("abundance") + 
  scale_y_log10() + 
  facet_grid(sample~model) 

##########################################################

# calculate the max delta score for sig sites by gene 
sig_gene_list <- sig_scored %>% 
  arrange(desc(delta_score)) %>% 
  dplyr::select(delta_score, gene_name, model, sample) %>% 
  group_by(gene_name, model, sample) %>% 
  summarise(max_delta_score = max(delta_score)) 

# make a histogram of max delta score by site by gene 
ggplot(sig_gene_list, aes(x = max_delta_score, fill = model)) + 
  geom_histogram(binwidth = 0.03, alpha = 0.6) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Maximum delta score per gene on significant genomic sites") + 
  xlab("Max delta score for gene") + ylab("abundance") + 
  scale_y_log10() + 
  facet_grid(sample~model) 

# plot the relationship between m6a and m5c max gene delta score on significant sites 
sig_gene_list %>% pivot_wider(names_from = "model", values_from = "max_delta_score") %>% 
  ggplot(., aes(x = m6a, y = m5c)) + 
  geom_point(color = "red", alpha = 0.5) + 
  geom_smooth(span = 0.2) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Relationship between score and delta-score for significant genomic sites") + 
  xlab("Maximum m6a delta_score for gene on significant site") + 
  ylab("Maximum m5c delta_score for gene on significant site") + 
  facet_grid(~sample)

# check the correlation of the above data
cor_in <- sig_gene_list %>% pivot_wider(names_from = "model", values_from = "max_delta_score") 
cor.test(cor_in$m6a, cor_in$m5c)

# check genes with max delta score over 0.25
cor_in %>% filter(m5c > 0.25) %>% dplyr::select(gene_name) %>% unique() %>% clipr::write_clip()
cor_in %>% filter(m6a > 0.25) %>% dplyr::select(gene_name) %>% unique() %>% clipr::write_clip()


####################################################################################################
####################################################################################################

# repeat the isoform-specific modification analysis but ignore condition, i.e. consider methylation changes on the same isoform between conditions 

# checkpoint
library(tidyverse)
library(scales)
library(ggrepel)
comp <- read_tsv("~/localGadiData/2022-01-31_run-annotate2-all-pval/comp.txt.gz", col_names = T)
merged_metadata <- read_tsv("~/localGadiData/2022-01-31_run-annotate2-all-pval/merged_metadata.txt.gz", col_names = T)

# add a score, prob * stoich 
scored <- comp %>% 
  mutate(score = prob * stoich) %>% 
  mutate(metagene = case_when(rel_pos < 1 ~ "5UTR",
                              (rel_pos > 1 & rel_pos < 2) ~ "CDS",
                              rel_pos > 2 ~ "3UTR")) %>% 
  mutate(metagene = factor(metagene, levels = c("5UTR", "CDS", "3UTR"))) %>% 
  group_by(chr, start, model) %>% 
  mutate(delta_score = max(score) - min(score)) %>% 
  ungroup()

# plot the relationship between score and delta score 
ggplot(scored %>% dplyr::select(score, delta_score, model, metagene) %>% na.omit(), aes(x = score, y = delta_score)) + 
  geom_bin2d(bins = 100) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Relationship between score and delta-score for all sites") + 
  xlab("Site score in transcriptome space (prob * stoich)") + ylab("Max delta score at genomic position between isoforms and conditions") + 
  facet_grid(metagene ~ model) 

# plot the distribution of delta score by metagene location
ggplot(scored %>% dplyr::select(delta_score, model, metagene) %>% na.omit(), aes(x = delta_score, fill = metagene)) + 
  geom_histogram(binwidth = 0.01, alpha = 0.6) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Delta score between isoforms across conditions") + 
  xlab("delta score") + ylab("abundance") + 
  scale_y_log10() + 
  facet_grid(metagene~model) 

##########################################################

# remake the same plots except for genomic positions considered significant in at least one condition 
sig_genome <- comp %>% 
  dplyr::select(chr, start, sample, model, prob) %>% 
  filter(prob > 0.9999) %>% 
  dplyr::select(-prob) %>% 
  unique()

sig_all <- inner_join(sig_genome, comp, by = c("chr", "start", "sample", "model"))

sig_scored <- sig_all %>% 
  mutate(score = prob * stoich) %>% 
  mutate(metagene = case_when(rel_pos < 1 ~ "5UTR",
                              (rel_pos > 1 & rel_pos < 2) ~ "CDS",
                              rel_pos > 2 ~ "3UTR")) %>% 
  mutate(metagene = factor(metagene, levels = c("5UTR", "CDS", "3UTR"))) %>% 
  group_by(chr, start, model) %>% 
  mutate(delta_score = max(score) - min(score)) %>% 
  ungroup()

# plot the relationship between score and delta score 
ggplot(sig_scored %>% dplyr::select(score, delta_score, model, metagene) %>% na.omit(), aes(x = score, y = delta_score)) + 
  geom_bin2d(bins = 40) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Relationship between score and delta-score for significant genomic sites") + 
  xlab("Site score in transcriptome space (prob * stoich)") + ylab("Max delta score at genomic position between isoforms and conditions") + 
  facet_grid(metagene ~ model) 

# plot the distribution of delta score by metagene location
ggplot(sig_scored %>% dplyr::select(delta_score, model, metagene) %>% na.omit(), aes(x = delta_score, fill = metagene)) + 
  geom_histogram(binwidth = 0.03, alpha = 0.6) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Delta score between isoforms across conditions for significant genomic sites") + 
  xlab("delta score") + ylab("abundance") + 
  scale_y_log10() + 
  facet_grid(metagene~model) 

##########################################################

# calculate the max delta score for sig sites by gene 
sig_gene_list <- sig_scored %>% 
  arrange(desc(delta_score)) %>% 
  dplyr::select(delta_score, gene_name, model) %>% 
  group_by(gene_name, model) %>% 
  summarise(max_delta_score = max(delta_score)) 

# make a histogram of max delta score by site by gene 
ggplot(sig_gene_list, aes(x = max_delta_score, fill = model)) + 
  geom_histogram(binwidth = 0.03, alpha = 0.6) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Maximum delta score per gene on significant genomic sites") + 
  xlab("Max delta score for gene") + ylab("abundance") + 
  scale_y_log10() + 
  facet_grid(~model) 

# plot the relationship between m6a and m5c max gene delta score on significant sites 
sig_gene_list %>% pivot_wider(names_from = "model", values_from = "max_delta_score") %>% 
  ggplot(., aes(x = m6a, y = m5c)) + 
  geom_point(color = "red", alpha = 0.5) + 
  geom_smooth(span = 0.2) + 
  theme_bw() + 
  theme(text = element_text(size=12)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Relationship between score and delta-score for significant genomic sites") + 
  xlab("Maximum m6a delta_score for gene on significant site") + 
  ylab("Maximum m5c delta_score for gene on significant site") 

# check the correlation of the above data
cor_in <- sig_gene_list %>% pivot_wider(names_from = "model", values_from = "max_delta_score") 
cor.test(cor_in$m6a, cor_in$m5c)

# check genes with max delta score over 0.4
cor_in %>% filter(m5c > 0.4) %>% dplyr::select(gene_name) %>% unique() %>% clipr::write_clip()
cor_in %>% filter(m6a > 0.30) %>% dplyr::select(gene_name) %>% unique() %>% clipr::write_clip()






