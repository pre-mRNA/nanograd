#!/usr/bin/env Rscript

# written by Kat Woodward on 2022-05-02 and edited by AJ Sethi on 2022-05-05 

# Aim: To identify differential gene expression and differential transcript expression in Agin's degradation data, 
# taking salmon transcript quantification of Agin's HEK293 DRS reads as input 

# to run, direct salmon_output to the correct path for salmon output on your machine 
salmon_output <- ("~/localGadiData/2022-05-05_KW-Salmon-DTU")

####################################################################################################
####################################################################################################

# start 

# load libraries 
library(tidyverse)
library(tximport)

# library(dplyr) 
# don't need to load this if you import tidyverse 

# library(ggplot2) 
# don't need to load this if you import tidyverse 

# library(BiocManager) 
# don't need to load this unless you're installing packages from BioConductor 

# list all files
# setwd("") 
# don't need to do this because you can specify path in list.files() - working directory should stay same as Rproj 

# specify path (in variable salmon_output) to list.files() directly 
files <- list.files(salmon_output, pattern = "quant.sf", full.names = TRUE, recursive = TRUE)
names(files) <- files
all(file.exists(files))

# make tx2gene file

# tx2gene is a dataframe not a file

# library(ensembldb)
# library(AnnotationHub)

ah <- AnnotationHub()
# query(ah, "EnsDb.Hsapiens")
# edb <- ah[["AH98047"]] # how did you get this? It's not in the query
# txs <- transcripts(edb, return.type = "DataFrame")
# txs
# tx2gene <- data.frame(txs$tx_id_version, txs$gene_id)

# import files 
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# load DESeq2
library(DESeq2)

#make colData
conds <- factor(c("Degraded", "Degraded", "Undegraded", "Undegraded"))
colData=data.frame(condition=conds)

dds <- DESeqDataSetFromTximport(txi, colData, ~condition)

##PCA plot
rld <- rlog(dds)
head(assay(rld))
plotPCA(rld,intgroup = "condition")

png("PCA_rlog_DDS.png")
svg("PCA_rlog_dds.svg")
dev.off()

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
plotPCA(vsd, intgroup="condition")

png("PCA_vst_DDS.png")
svg("PCA_vst_dds.svg")
dev.off()

##DESEQ2
dds <- DESeq(dds)
res <- results(dds)
which(res$padj < 0.05)

# extract a results from DESEQ
res_1 <- as_tibble(results(dds,contrast=c("condition", "Degraded", "Undegraded")), rownames = "geneID") %>%
  dplyr::rename(LFC = 3) %>%
  na.omit()

# use biomart to add geneIDs to output
library(biomaRt)
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)


# add geneIDs to contrast
ens <- res_1$geneID
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(res_1$geneID, names(symbols))]
res_1$geneName <- symbols
res_1 <- res_1 %>% relocate(geneName, .after = "geneID")
write_tsv(res_1, "2022-05-01_DEG_vs_NONDEG.txt.gz")

library(EnhancedVolcano)
EnhancedVolcano(res_1,
                lab = res_1$geneName,
                x = "LFC",
                y = "padj",
                title = 'DEG between Degraded vs Undegraded',
                ylim = c(0, -log10(10e-12)),
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

png("Volcano_plot_standard_DESEQ.png", width = 800, height = 600)
svg("Volcano_plot_standard_DESEQ.svg")
dev.off()


