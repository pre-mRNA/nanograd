###DESEQ2 on degraded vs non-degraded RNA

library(DESeq2)
library(tidyverse)

# set path to counts file

# kat's laptop
# counts_file <- "C:/Users/User/Documents/ANU/Year 4/Semester 2/PhD/Nanograd/Featurecounts output/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

# AJ's JCSMR iMac 
# counts_file <- "/Users/AJlocal/localGadiData/2022-06-22_HEK293-degradation-first4-AR_liqa-genome-alignments_BK/2022-07-05_ARdegradation-genomic-featurecounts-AS.txt"

# import counts as a tibble 
raw_counts <- read_tsv(counts_file, col_names = T, skip = 1, col_types = "fccccdfdddd") %>% 
  dplyr::rename(gene_id = Geneid, deg_rep1 = 8, deg_rep2 = 9, wt_rep1 = 10, wt_rep2 = 11) 

# make a counts table 
counts_matrix <- raw_counts %>% 
  dplyr::select(gene_id, wt_rep1, wt_rep2, deg_rep1, deg_rep2) %>% 
  remove_rownames() %>% column_to_rownames(var="gene_id") %>% 
  as.matrix()

# define condition 
conds <- factor(c("Control", "Control", "Degraded", "Degraded"))
# summary(conds)

# add column data 
colData=data.frame(condition=conds)

# run DESEq2
dds <- DESeqDataSetFromMatrix(counts_matrix, colData, design = ~condition)

##PCA
rld <- rlog(dds)
head(assay(rld))
plotPCA(rld,intgroup = "condition")

#DE analysis
dds <- DESeq(dds)
res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)

# extract DESeq2 contrast data as tibble 
results_tibble <- as_tibble(results(dds,contrast=c("condition", "Degraded", "Control")), rownames = "geneID") %>% 
  dplyr::rename(LFC = 3) %>% 
  na.omit()

# get a list of geneIDs and geneNames to annotate DESeq2 results 
gene_key <- raw_counts %>% select(gene_id, gene_name) %>% unique() %>% rename(geneID = gene_id)

# ane merge 
results_annotated <- inner_join(gene_key, results_tibble, by = "geneID")

# plot for KW 

##########################################################################################
##########################################################################################
##########################################################################################

# KW old code


# first contrast, WT_E15 vs WT_E18, as WT15-18

# put second condition first in the contrast such that upregulated genes are upregulated in the first condition 
WT_KO <- as_tibble(results(dds,contrast=c("condition", "KO", "WT")), rownames = "geneID") %>% 
  dplyr::rename(LFC = 3) %>% 
  na.omit()

##########################################################################################
##########################################################################################
##########################################################################################


##ENhanced Volcano

library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "padj",
                title = 'DEG between Decayed vs Non-decayed',
                xlab = "log2FoldChange",
                ylim = c(0, -log10(10e-12)),
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
