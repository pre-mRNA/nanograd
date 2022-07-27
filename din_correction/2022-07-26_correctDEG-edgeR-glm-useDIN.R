<<<<<<< HEAD
#!/bin/bash

# written by AJ abd BK on 2022-07-25

###########################################################

# load libraries 
library(tidyverse)
library(edgeR)
library(EnhancedVolcano)

# JCSMR iMac: 
#counts_file <- "/Users/AJlocal/localGadiData/2022-06-22_HEK293-degradation-first4-AR_liqa-genome-alignments_BK/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

# Bhavika's laptop: 
 counts_file <- "d:/Users/Sujata Kumar/Desktop/Project/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

########## liqa transcript abundance

# JCSMR iMac
#liqa_transcript_counts <- "/Users/AJlocal/localGadiData/2022-07-17_LIQA_isoform_counts/undegraded_hek293_pass1_primary.txt"

 liqa_transcript_counts <- "//wsl.localhost/Ubuntu/home/bhavika_kumar/localGadiData/isoform_expression/undegraded_hek293_pass1_primary.txt"

########## biomart transcript lengths

# JCSMR iMac
#bm_tx_lengths <- "/Users/AJlocal/localGadiData/2022-07-17_LIQA_isoform_counts/2022-07-17_biomart-human-transcript-lengths.txt"

 bm_tx_lengths <- "d:/Users/Sujata Kumar/Desktop/Project/2022-07-17_biomart-human-transcript-lengths.txt"

# import and preprocess the counts data from featureCounts + uLTRA 

# read in the counts file 
raw_counts <- read_tsv(counts_file, col_names = T, skip = 1, col_types = "fccccdfdddd") %>%
  dplyr::rename(gene_id = Geneid, deg_rep1 = 8, deg_rep2 = 9, wt_rep1 = 10, wt_rep2 = 11) %>%
  mutate(across(c(Start, End, Chr, Strand), gsub, pattern = ";.*", replacement = "")) %>% 
  mutate(chr = Chr) %>% 
  select(gene_id, gene_name, chr, wt_rep1, wt_rep2, deg_rep1, deg_rep2)


####################################################################################################

# use LIQA output to calculate the most likely isoform and transcript length for each gene 

liqa <- read_tsv(liqa_transcript_counts, col_names = T, col_types = "ffddd") %>% 
  rename(transcript_id = IsoformName) %>% 
  rename(transcript_abundance = ReadPerGene_corrected) %>% 
  select(transcript_id, transcript_abundance)

# from bioMart, for humans, download a table of gene id, transcript_id, and transcirpt length (including UTRs and CDS)
# uploaded to OneDrive as 2022-07-17_biomart-transcript-lengths.txt

bm <- read_tsv(bm_tx_lengths, col_names = T) %>% 
  rename(gene_id = 1, transcript_id = 3, transcript_length = 5) %>% 
  select(gene_id, transcript_id, transcript_length)

# merge the data 
merged <- left_join(bm, liqa, by = "transcript_id") %>% 
  group_by(gene_id) %>% 
  slice(which.max(transcript_abundance))

# figure out which genes are not in merged
`%ni%` <- Negate(`%in%`)
outside_genes <- bm %>% filter(gene_id %ni% merged$gene_id) %>% 
  group_by(gene_id) %>% 
  summarise(transcript_length = mean(transcript_length))

# merge merged with outside_genes to have a length for each gene 
final_length_key <- bind_rows(merged %>% select(gene_id, transcript_length), outside_genes)

####################################################################################################
# attach gene length to raw_counts 

counts <- inner_join(raw_counts, final_length_key, by = "gene_id") %>% 
  select(gene_id, wt_rep1, wt_rep2, deg_rep1, deg_rep2)

name_key <- inner_join(raw_counts, final_length_key, by = "gene_id") %>% select(gene_id, gene_name)

####################################################################################################

# convert raw counts to matrix format 

  
counts_import_matrix <- as.data.frame(counts)
data_clean <- counts_import_matrix[,-1]
row.names(data_clean) <- counts_import_matrix[,1]
  
# log transform the counts 
cpm_log <- cpm(data_clean, log = TRUE)
  
# get median count
median_log2_cpm <- apply(cpm_log, 1, median)
  
  
####################################################################################################

# first, do edgeR GLM without including DIN 

expr_cutoff <- 5

# filter for genes where the median log2 expression is greater than the cutoff 
data_clean <- data_clean[median_log2_cpm > expr_cutoff, ]

# print how many genes we still have in the analysis 
rows <- nrow(data_clean)
print(paste("you have ", rows, " genes remaining in the analysis"))

# define the groups for comparison 
group <- c("control","control","degraded", "degraded")

# design the experiment 
design <- model.matrix(~group)
design

# make DGE object 
y <- DGEList(counts=data_clean,group=factor(group))

# normalize counts using edgeR 
y <- calcNormFactors(y)
y$samples

# estimate dispersion 
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

# extract the differential expression data and convert it to a tibble and add gene names 
DEtib <- lrt[[14]] %>% 
  as_tibble(rownames = "gene_id") %>% 
  dplyr::rename(p.raw = PValue) %>% 
  right_join(name_key, ., by = "gene_id")

# correct the p-value
DEtib$p.adj <- p.adjust(DEtib$p.raw, method = "fdr")

total_genes <- nrow(DEtib)
sig_genes <- nrow(DEtib %>% filter(p.adj < 0.05 & abs(logFC) > 1))
specificity <- total_genes/(total_genes + sig_genes)
print(paste(total_genes, " total genes; ", sig_genes, " significant genes;", specificity, " specificity"))

# assign to a variable, uncorrected 
uncorrected_DEgenes <- DEtib

##############################################################################################################

# now repeat, but incorporate DIN 

# take the expression cutoff from the input 
expr_cutoff <- 5

# filter for genes where the median log2 expression is greater than the cutoff 
data_clean <- data_clean[median_log2_cpm > expr_cutoff, ]

# print how many genes we still have in the analysis 
rows <- nrow(data_clean)
print(paste("you have ", rows, " genes remaining in the analysis"))

# define the groups for comparison 
group <- c("control","control","degraded", "degraded")

# add DIN 
din <- c(11.8, 11.8, 8.24, 8.11)

# design the experiment 
design <- model.matrix(~group+din)
design

# make DGE object 
y <- DGEList(counts=data_clean,group=factor(group))

# normalize counts using edgeR 
y <- calcNormFactors(y)
y$samples

# estimate dispersion 
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

# extract the differential expression data and convert it to a tibble and add gene names 
DEtib <- lrt[[14]] %>% 
  as_tibble(rownames = "gene_id") %>% 
  dplyr::rename(p.raw = PValue) %>% 
  right_join(name_key, ., by = "gene_id")

# correct the p-value
DEtib$p.adj <- p.adjust(DEtib$p.raw, method = "fdr")


total_genes <- nrow(DEtib)
sig_genes <- nrow(DEtib %>% filter(p.adj < 0.05 & abs(logFC) > 1))
specificity <- total_genes/(total_genes + sig_genes)
print(paste(total_genes, " total genes; ", sig_genes, " significant genes;", specificity, " specificity"))

# assign to a variable, corrected 
corrected_DEgenes <- DEtib

##############################################################################################################

# join the corrected and uncorrected analyses 
joined_DEgenes <- full_join(uncorrected_DEgenes, corrected_DEgenes, by = "gene_id", suffix=c("_uncorrected", "_corrected")) %>% 
  mutate(deltaCPM = logCPM_corrected - logCPM_uncorrected, 
         deltapadj = p.adj_corrected -  p.adj_uncorrected,
         deltalogFC = logFC_corrected - logFC_uncorrected)

ggplot(joined_DEgenes, aes(x = deltaCPM)) + geom_histogram() + scale_y_log10()
ggplot(joined_DEgenes, aes(x = deltaCPM)) + stat_ecdf()

ggplot(joined_DEgenes, aes(x = deltalogFC)) + geom_histogram()
ggplot(joined_DEgenes, aes(x = deltalogFC)) + stat_ecdf()

ggplot(joined_DEgenes, aes(x = deltapadj)) + geom_histogram()
ggplot(joined_DEgenes, aes(x = deltapadj)) + stat_ecdf()

ggplot(joined_DEgenes, aes(x = deltalogFC, y = deltapadj)) + geom_point() 

ggplot(joined_DEgenes, aes(x = deltaCPM, y = deltapadj)) + geom_point()

##############################################################################################################

# adding gene length to joined_DEgenes

merged_DEgenes <- full_join(joined_DEgenes,final_length_key, by="gene_id")

ggplot(merged_DEgenes, aes(x=deltaCPM, y=transcript_length))+ geom_point() 

ggplot(merged_DEgenes, aes(x=deltapadj, y=transcript_length))+ geom_point() 

ggplot(merged_DEgenes, aes(x=deltalogFC, y=transcript_length))+ geom_point()+ scale_y_log10()

cor.test(merged_DEgenes$deltalogFC,merged_DEgenes$transcript_length)
cor.test(merged_DEgenes$deltaCPM, merged_DEgenes$transcript_length)
cor.test(merged_DEgenes$deltapadj,merged_DEgenes$transcript_length)

##############################
EnhancedVolcano(merged_DEgenes,
                lab = merged_DEgenes$gene_name_corrected,
                x = "deltalogFC",
                y = "deltapadj",
                title = '',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.1,
                pointSize = 2.0,
                labSize = 4,
                colAlpha = 0.5,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = FALSE,
                widthConnectors = 0.5, 
                typeConnectors = "closed",
                endsConnectors = "first",
                lengthConnectors = unit(0.01, "npc"),
                colConnectors = "grey10")
            

EnhancedVolcano(uncorrected_DEgenes,
                lab = uncorrected_DEgenes$gene_name,
                x = "logFC",
                y = "p.adj",
                title = '',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.1,
                pointSize = 2.0,
                labSize = 4,
                colAlpha = 0.5,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = FALSE,
                widthConnectors = 0.5, 
                typeConnectors = "closed",
                endsConnectors = "first",
                lengthConnectors = unit(0.01, "npc"),
                colConnectors = "grey10")

######################################################

# save volcano plot for edgeR GLM degradation experiment 

# open a plotting device
dev.new(width = 9, height = 4.5, noRStudioGD = TRUE, unit = "cm")

# make the uncorrected volcano plot 
EnhancedVolcano(uncorrected_DEgenes,
                lab = uncorrected_DEgenes$gene_name,
                x = "logFC",
                y = "p.adj",
                title = '',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.1,
                pointSize = 2.0,
                labSize = 4,
                colAlpha = 0.5,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 4.0,
                drawConnectors = FALSE,
                widthConnectors = 0.5, 
                typeConnectors = "closed",
                endsConnectors = "first",
                lengthConnectors = unit(0.01, "npc"),
                colConnectors = "grey10")

# export the plot 
ggsave(filename = "<name>", plot = last_plot(), width = 5, height = 5, units = "cm", dpi = 300, scale = 2)
ggsave(filename = "<name", plot = last_plot(), width = 5, height = 5, units = "cm", dpi = 300, scale = 2)
=======
#!/bin/bash

# written by AJ abd BK on 2022-07-25

###########################################################

# load libraries 
library(tidyverse)
library(edgeR)
library(EnhancedVolcano)

# JCSMR iMac: 
counts_file <- "/Users/AJlocal/localGadiData/2022-06-22_HEK293-degradation-first4-AR_liqa-genome-alignments_BK/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

# Bhavika's laptop: 
# counts_file <- "d:/Users/Sujata Kumar/Desktop/Project/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

########## liqa transcript abundance

# JCSMR iMac
liqa_transcript_counts <- "/Users/AJlocal/localGadiData/2022-07-17_LIQA_isoform_counts/undegraded_hek293_pass1_primary.txt"

# liqa_transcript_counts <- "//wsl.localhost/Ubuntu/home/bhavika_kumar/localGadiData/isoform_expression/undegraded_hek293_pass1_primary.txt"

########## biomart transcript lengths

# JCSMR iMac
bm_tx_lengths <- "/Users/AJlocal/localGadiData/2022-07-17_LIQA_isoform_counts/2022-07-17_biomart-human-transcript-lengths.txt"

# bm_tx_lengths <- "d:/Users/Sujata Kumar/Desktop/Project/2022-07-17_biomart-human-transcript-lengths.txt"

# import and preprocess the counts data from featureCounts + uLTRA 

# read in the counts file 
raw_counts <- read_tsv(counts_file, col_names = T, skip = 1, col_types = "fccccdfdddd") %>%
  dplyr::rename(gene_id = Geneid, deg_rep1 = 8, deg_rep2 = 9, wt_rep1 = 10, wt_rep2 = 11) %>%
  mutate(across(c(Start, End, Chr, Strand), gsub, pattern = ";.*", replacement = "")) %>% 
  mutate(chr = Chr) %>% 
  select(gene_id, gene_name, chr, wt_rep1, wt_rep2, deg_rep1, deg_rep2)


####################################################################################################

# use LIQA output to calculate the most likely isoform and transcript length for each gene 

liqa <- read_tsv(liqa_transcript_counts, col_names = T, col_types = "ffddd") %>% 
  rename(transcript_id = IsoformName) %>% 
  rename(transcript_abundance = ReadPerGene_corrected) %>% 
  select(transcript_id, transcript_abundance)

# from bioMart, for humans, download a table of gene id, transcript_id, and transcirpt length (including UTRs and CDS)
# uploaded to OneDrive as 2022-07-17_biomart-transcript-lengths.txt

bm <- read_tsv(bm_tx_lengths, col_names = T) %>% 
  rename(gene_id = 1, transcript_id = 3, transcript_length = 5) %>% 
  select(gene_id, transcript_id, transcript_length)

# merge the data 
merged <- left_join(bm, liqa, by = "transcript_id") %>% 
  group_by(gene_id) %>% 
  slice(which.max(transcript_abundance))

# figure out which genes are not in merged
`%ni%` <- Negate(`%in%`)
outside_genes <- bm %>% filter(gene_id %ni% merged$gene_id) %>% 
  group_by(gene_id) %>% 
  summarise(transcript_length = mean(transcript_length))

# merge merged with outside_genes to have a length for each gene 
final_length_key <- bind_rows(merged %>% select(gene_id, transcript_length), outside_genes)

####################################################################################################
# attach gene length to raw_counts 

counts <- inner_join(raw_counts, final_length_key, by = "gene_id") %>% 
  select(gene_id, wt_rep1, wt_rep2, deg_rep1, deg_rep2)

name_key <- inner_join(raw_counts, final_length_key, by = "gene_id") %>% select(gene_id, gene_name)

####################################################################################################

# convert raw counts to matrix format 

  
counts_import_matrix <- as.data.frame(counts)
data_clean <- counts_import_matrix[,-1]
row.names(data_clean) <- counts_import_matrix[,1]
  
# log transform the counts 
cpm_log <- cpm(data_clean, log = TRUE)
  
# get median count
median_log2_cpm <- apply(cpm_log, 1, median)
  
  
####################################################################################################

# first, do edgeR GLM without including DIN 

expr_cutoff <- 5

# filter for genes where the median log2 expression is greater than the cutoff 
data_clean <- data_clean[median_log2_cpm > expr_cutoff, ]

# print how many genes we still have in the analysis 
rows <- nrow(data_clean)
print(paste("you have ", rows, " genes remaining in the analysis"))

# define the groups for comparison 
group <- c("control","control","degraded", "degraded")

# design the experiment 
design <- model.matrix(~group)
design

# make DGE object 
y <- DGEList(counts=data_clean,group=factor(group))

# normalize counts using edgeR 
y <- calcNormFactors(y)
y$samples

# estimate dispersion 
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

# extract the differential expression data and convert it to a tibble and add gene names 
DEtib <- lrt[[14]] %>% 
  as_tibble(rownames = "gene_id") %>% 
  dplyr::rename(p.raw = PValue) %>% 
  right_join(name_key, ., by = "gene_id")

# correct the p-value
DEtib$p.adj <- p.adjust(DEtib$p.raw, method = "fdr")

total_genes <- nrow(DEtib)
sig_genes <- nrow(DEtib %>% filter(p.adj < 0.05 & abs(logFC) > 1))
specificity <- total_genes/(total_genes + sig_genes)
print(paste(total_genes, " total genes; ", sig_genes, " significant genes;", specificity, " specificity"))

# assign to a variable, uncorrected 
uncorrected_DEgenes <- DEtib

##############################################################################################################

# now repeat, but incorporate DIN 

# take the expression cutoff from the input 
expr_cutoff <- 5

# filter for genes where the median log2 expression is greater than the cutoff 
data_clean <- data_clean[median_log2_cpm > expr_cutoff, ]

# print how many genes we still have in the analysis 
rows <- nrow(data_clean)
print(paste("you have ", rows, " genes remaining in the analysis"))

# define the groups for comparison 
group <- c("control","control","degraded", "degraded")

# add DIN 
din <- c(11.8, 11.8, 8.24, 8.11)

# design the experiment 
design <- model.matrix(~group+din)
design

# make DGE object 
y <- DGEList(counts=data_clean,group=factor(group))

# normalize counts using edgeR 
y <- calcNormFactors(y)
y$samples

# estimate dispersion 
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

# extract the differential expression data and convert it to a tibble and add gene names 
DEtib <- lrt[[14]] %>% 
  as_tibble(rownames = "gene_id") %>% 
  dplyr::rename(p.raw = PValue) %>% 
  right_join(name_key, ., by = "gene_id")

# correct the p-value
DEtib$p.adj <- p.adjust(DEtib$p.raw, method = "fdr")


total_genes <- nrow(DEtib)
sig_genes <- nrow(DEtib %>% filter(p.adj < 0.05 & abs(logFC) > 1))
specificity <- total_genes/(total_genes + sig_genes)
print(paste(total_genes, " total genes; ", sig_genes, " significant genes;", specificity, " specificity"))

# assign to a variable, corrected 
corrected_DEgenes <- DEtib

##############################################################################################################

# join the corrected and uncorrected analyses 
joined_DEgenes <- full_join(uncorrected_DEgenes, corrected_DEgenes, by = "gene_id", suffix=c("_uncorrected", "_corrected")) %>% 
  mutate(deltaCPM = logCPM_corrected - logCPM_uncorrected, 
         deltapadj = p.adj_corrected -  p.adj_uncorrected,
         deltalogFC = logFC_corrected - logFC_uncorrected)

ggplot(joined_DEgenes, aes(x = deltaCPM)) + geom_histogram() + scale_y_log10()
ggplot(joined_DEgenes, aes(x = deltaCPM)) + stat_ecdf()

ggplot(joined_DEgenes, aes(x = deltalogFC)) + geom_histogram()
ggplot(joined_DEgenes, aes(x = deltalogFC)) + stat_ecdf()

ggplot(joined_DEgenes, aes(x = deltapadj)) + geom_histogram()
ggplot(joined_DEgenes, aes(x = deltapadj)) + stat_ecdf()

ggplot(joined_DEgenes, aes(x = deltalogFC, y = deltapadj)) + geom_point() 
>>>>>>> 82d55c9103a0ca344523916949e7e0e352063970
