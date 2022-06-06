# load tidyverse 
library(tidyverse)
library(EnhancedVolcano)
library(biomaRt)

# read in the data (metafeature)
import_data <- read_tsv("/Users/asethi/localGadiData/2022-05-11_degradation-genomic-alignments/featureCounts_primary.txt", 
                        col_names = T, 
                        comment = "#", 
                        col_types = "fcnncnffnnnn") %>% 
  dplyr::rename(geneID = 1, chr = 2, start = 3, end = 4, strand = 5, length = 6, biotype = 7, name = 8, deg_1 = 9, deg_2 = 10, wt_1 = 11, wt_2 = 12) %>% 
  dplyr::select(geneID, wt_1, wt_2, deg_1, deg_2)

# convert tibble to data.frame for DESeq2
import_df <- data.frame(import_data)
rownames(import_df) <- import_df$geneID
import_df <- import_df[,-1]

# trim extra data and convert to matrix 
countdata <- import_df[]
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition
condition <- factor(c("WT",
                      "WT",
                      "deg",
                      "deg"))

# Analysis with DESeq2 ----------------------------------------------------

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Plot dispersions
plotDispEsts(dds, main="Dispersion plot")


# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# make sample distance heatmap 
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(20, 20), main="Sample Distance Matrix", 
          cexRow = 2, cexCol = 2)


# Principal components analysis

# make sample distance heatmap 
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# plot PCA
DESeq2::plotPCA(rld, intgroup="condition")+
  theme_bw() +
  theme(text = element_text(size=24)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  geom_label(aes(label = name)) # to add names


# preparing data for volcano plots 
coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds


dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized=TRUE)
cntnap2_counts <- normalized_counts %>% as_tibble(rownames = "geneID") %>% filter(grepl("ENSMUSG00000039419", geneID)) 
clipr::write_clip(cntnap2_counts)

# check the difference between raw and normalized counts for cntnap2
cntnap2_raw <- import_data %>% filter(grepl("ENSMUSG00000039419", geneID))

cntnap2_merged <- bind_rows(cntnap2_counts, cntnap2_raw) %>% dplyr::select(-chr, -start, -end, -strand, -biotype, -geneID, -length) %>% rownames_to_column %>%
  gather(variable, value, -rowname) %>% 
  spread(rowname, value) %>% 
  dplyr::rename(library = 1, count = 3, normCount = 2) %>% 
  type_convert(col_types = "cdd") %>% 
  dplyr::mutate(
    mutant = case_when(
      str_detect(library, "M") ~ TRUE, 
      TRUE ~ FALSE))

cor.test(cntnap2_merged$count, cntnap2_merged$normCount)
ggplot(cntnap2_merged, aes(x = count, y = normCount, color = mutant)) + geom_point(size = 5) +  theme_bw() +
  theme(text = element_text(size=24)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Raw and norm counts for Cntnap2") + 
  xlab("Counts in genomic alignment") + ylab("DESeq2 normalised counts") + 
  geom_smooth(method='lm', formula= y~x)

# also impport counts from the latest transcriptome alignments 
# cd /g/data/xc17/akanksha/mouse_brains/minimap2_nopseudo 
# for i in $(find $(pwd) -name "*.bam"); do printf "${i##*/}\t$(samtools view -F 2308 $i | grep -E 'ENSMUST00000207647.2|ENSMUST00000114641.8|ENSMUST00000196561.2|ENSMUST00000150737.2|ENSMUST00000199100.5|ENSMUST00000060839.8' | wc -l)\n" >> ${scratch}/counts.txt & done 
# manually repair names and save to /data/2021-09-29_cntnap2-transcriptomeCounts.txt

cntnap2_transcriptome_counts <- read_tsv("./data/2021-09-29_cntnap2-transcriptomeCounts.txt", col_names = F) %>% dplyr::rename(lib = 1, count = 2)

cntnap2_merged <- inner_join(cntnap2_raw %>% dplyr::select(-chr, -start, -end, -strand, -biotype, -geneID, -length) %>% rownames_to_column %>% gather(variable, value, -rowname) %>% spread(rowname, value) %>% dplyr::rename(lib = 1, geneCount = 2), 
                             cntnap2_transcriptome_counts %>% dplyr::rename(lib = 1, transcriptomeCount = 2)) %>% 
  dplyr::mutate(
    mutant = case_when(
      str_detect(lib, "M") ~ TRUE, 
      TRUE ~ FALSE))


cor.test(cntnap2_merged$geneCount, cntnap2_merged$transcriptomeCount)
ggplot(cntnap2_merged, aes(x = geneCount, y = transcriptomeCount, color = mutant)) + geom_point(size = 5) + theme_bw() +
  theme(text = element_text(size=24)) + 
  theme(plot.title = element_text(hjust=0.5)) +  
  ggtitle("Genomic and transcriptomic counts for Cntnap2") + 
  xlab("Counts in genomic alignment") + ylab("Counts in transcriptomic alignment")


##########################################################################################
##########################################################################################
##########################################################################################

# first contrast, WT_E15 vs WT_E18, as WT15-18

# put second condition first in the contrast such that upregulated genes are upregulated in the first condition 
WT15_18 <- as_tibble(results(dds,contrast=c("condition", "WT_E18", "WT_E15")), rownames = "geneID") %>% 
  dplyr::rename(LFC = 3) %>% 
  na.omit()


WT15_18

WT15_18 %>% nrow() # 1469 
WT15_18 %>% filter(padj < 0.1) %>% nrow() # 782
WT15_18 %>% filter(padj < 0.05) %>% nrow() # 782
WT15_18 %>% filter(padj < 0.1) %>% filter(abs(LFC) > 1) %>% nrow() # 499
WT15_18 %>% filter(padj < 0.05) %>% filter(abs(LFC) > 1) %>% nrow() # 390 

# add geneIDs to contrast a (WT_E15 vs WT_E18)
ens <- WT15_18$geneID
symbols <- mapIds(org.Mm.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')

symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(WT15_18$geneID, names(symbols))]
WT15_18$geneName <- symbols

WT15_18 <- WT15_18 %>% relocate(geneName, .after = "geneID")
write_tsv(WT15_18, "./data/2022-04-26_WT15_18.txt.gz")

# make p-value histogram 
WT15_18 %>% 
  dplyr::select(pvalue, padj) %>% 
  pivot_longer(cols = starts_with("p"), names_to = "stat", values_to = "value") %>% 
  mutate(stat = factor(stat, levels = c("pvalue", "padj"))) %>% 
  ggplot(aes(x = value, fill = stat)) + 
  geom_histogram(bins = 100, alpha = 1, position = "identity") + 
  facet_wrap(~stat, nrow = 2) +
  theme_bw() +
  theme(text = element_text(size=24)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle("Raw and adjusted p for WT_E15 vs _E18") + 
  xlab("p") + ylab("count")


# determine significant genes 
sigGenes <- WT15_18 %>% filter(padj < 0.1) %>% filter(abs(LFC) > 1) %>% dplyr::select(geneName)
clipr::write_clip(sigGenes)

sigGenes2 <- WT15_18 %>% filter(padj < 0.05) %>% filter(abs(LFC) > 1) %>% dplyr::select(geneName)
clipr::write_clip(sigGenes)

upGenes <- WT15_18 %>% filter(padj < 0.1) %>% filter(LFC > 1) %>% dplyr::select(geneName)
clipr::write_clip(upGenes)

downGenes <- WT15_18 %>% filter(padj < 0.1) %>% filter(LFC < -1) %>% dplyr::select(geneName)
clipr::write_clip(downGenes)

allGenes <- WT15_18 %>% dplyr::select(geneName)
clipr::write_clip(allGenes)

# make a volcano plot for WT15_18
EnhancedVolcano(WT15_18,
                lab = WT15_18$geneName,
                x = "LFC",
                y = "padj",
                title = 'DEG between WT15 vs WT18, padj < 0.1',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
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