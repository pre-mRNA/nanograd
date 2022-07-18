#!/bin/bash

# written by AJ Sethi on 2022-07-17
# aim: 

####################################################################################################
####################################################################################################

# modules
library(tidyverse)
library(edgeR)

####################################################################################################

# import and preprocess the counts data from featureCounts + uLTRA 

# path to counts file 

# JCSMR iMac: 
counts_file <- "/Users/AJlocal/localGadiData/2022-06-22_HEK293-degradation-first4-AR_liqa-genome-alignments_BK/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

# Bhavika's laptop: 
# counts_file <- "d:/Users/Sujata Kumar/Desktop/Project/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

# read in the counts file 
raw_counts <- read_tsv(counts_file, col_names = T, skip = 1, col_types = "fccccdfdddd") %>%
  dplyr::rename(gene_id = Geneid, deg_rep1 = 8, deg_rep2 = 9, wt_rep1 = 10, wt_rep2 = 11) %>%
  mutate(across(c(Start, End, Chr, Strand), gsub, pattern = ";.*", replacement = "")) %>% 
  mutate(chr = Chr) %>% 
  select(gene_id, gene_name, chr, wt_rep1, wt_rep2, deg_rep1, deg_rep2)

####################################################################################################

# use LIQA output to calculate the most likely isoform and transcript length for each gene 

liqa <- read_tsv("/Users/AJlocal/localGadiData/2022-07-17_LIQA_isoform_counts/undegraded_hek293_pass1_primary.txt", col_names = T, col_types = "ffddd") %>% 
  rename(transcript_id = IsoformName) %>% 
  rename(transcript_abundance = ReadPerGene_corrected) %>% 
  select(transcript_id, transcript_abundance)

# from bioMart, for humans, download a table of gene id, transcript_id, and transcirpt length (including UTRs and CDS)
# uploaded to OneDrive as 2022-07-17_biomart-transcript-lengths.txt

bm <- read_tsv("/Users/AJlocal/localGadiData/2022-07-17_LIQA_isoform_counts/2022-07-17_biomart-human-transcript-lengths.txt", col_names = T) %>% 
  rename(gene_id = 1, transcript_id = 3, transcript_length = 5) %>% 
  select(gene_id, transcript_id, transcript_length)

# merge the data 
merged <- left_join(bm, liqa, by = "transcript_id") %>% 
  group_by(gene_id) %>% 
  slice(which.max(transcript_abundance))

# figre out which genes are not in merged ]
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

# convert to matrix 
counts_import_matrix <- as.data.frame(counts)
counts_matrix <- counts_import_matrix[,-1]
row.names(counts_matrix) <- counts_import_matrix[,1]

# create list and calculate library sizes 
group <- c("control","control","degraded", "degraded")

# make DGE object 
d <- DGEList(counts=counts_matrix,group=factor(group))
d

# Normalization
dt <- calcNormFactors(d,method="TMM")

# make a PCA plot 
plotMDS(dt, gene.selection="common")

# Estimate dispersion
d1 <- estimateCommonDisp(dt, verbose=T)
d1 <- estimateTagwiseDisp(d1)

# run edgeR
design.mat <- model.matrix(~ 0 + dt$samples$group)
colnames(design.mat) <- levels(dt$samples$group)
d2 <- estimateGLMCommonDisp(dt,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="auto")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2,design.mat)


# do edgeR using GLM
fit <- glmFit(d2,design.mat)

# plot BCV 
plotBCV(d2)

# extract results
lrt12_raw <- glmLRT(fit,contrast=c(1,-1))

lrt12 <- glmLRT(fit,contrast=c(-1,1))[[14]] %>% 
  as_tibble(rownames = "geneName")


##############################################################################################################

# note: Figure out how to filter for expression level, we need a better method 

# add a filter for expression level 

filt <- lrt12 %>% filter(logCPM > 2.75)

ggplot(filt, aes(x = PValue)) + geom_histogram()

####################################################################################################

# add back gene names, where they exist 

named <- inner_join(filt %>% rename(gene_id = geneName), name_key, by = "gene_id")
named$p.adj <- p.adjust(named$PValue, method = "fdr")

named %>% filter(p.adj < 0.1)

##############################################################################################################

# convert the output to data.table for enhanced volcano
lrt12_dt <- as.data.frame(lrt12_names[,-1])
names2 <- make.unique(lrt12_dt$geneName)
rownames(lrt12_dt) <- names2

library(EnhancedVolcano)

EnhancedVolcano(named,
                lab = named$gene_name,
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
                colConnectors = "grey10",
                ylim = c(0,30),
                xlim = c(-3,3))

##############################################################################################################





