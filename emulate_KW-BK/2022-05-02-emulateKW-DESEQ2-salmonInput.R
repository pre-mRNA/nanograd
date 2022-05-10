#!/usr/bin/env Rscript

# written by Kat Woodward on 2022-05-02 and edited by AJ Sethi on 2022-05-05 

# Aim: To identify differential gene expression and differential transcript expression in Agin's degradation data, 
# taking salmon transcript quantification of Agin's HEK293 DRS reads as input 

# to run, direct salmon_output to the correct path for salmon output on your machine 
salmon_output <- ("~/localGadiData/2022-05-05_KW-Salmon-DTU") # AJ's iMac

# also provide the human annotation in GTF format 
human_annotation <- "~/localGadiData/2022-05-05_KW-Salmon-DTU/Homo_sapiens_transcriptsOnly_GRCh38.104.chr.gtf"

# AS: I filtered the annotation for 'transcript' categories before downloading it to save space on my Mac 
# The filtering command I used is: 
# $ cat all.gtf | awk '($3 == "transcript"){print}' > transcripts.gtf


####################################################################################################
####################################################################################################

# start 

# load libraries 
library(tidyverse)
library(tximport)

# library(dplyr) 
# don't need to load this if you import tidyverse 

# library(ggplot2) 
# AS: don't need to load this if you import tidyverse 

# library(BiocManager) 
# AS: don't need to load this unless you're installing packages from BioConductor 

# list all files
# setwd("") 
# AS: don't need to do this because you can specify path in list.files() - working directory should stay same as Rproj 

# AS: specify path (in variable salmon_output) to list.files() directly 
files <- list.files(salmon_output, pattern = "quant.sf", full.names = TRUE, recursive = TRUE)
names(files) <- files
all(file.exists(files))

##################################################
##################################################

# make tx2gene file
# AS: note that tx2gene is a dataframe not a file

# library(ensembldb)
# library(AnnotationHub)
# AS: I don't like this approach because I dont' want to download this to some random location in my local machine 
# AS: It's also harder to version control because annotationhub is continously updated, 
# e.g. because the version of R on my computer was slightly old, I couldn't access ah[["AH98047"]]

# ah <- AnnotationHub()
# query(ah, "EnsDb.Hsapiens")
# edb <- ah[["AH98047"]] # how did you get this? It's not in the query
# txs <- transcripts(edb, return.type = "DataFrame")
# txs
# tx2gene <- data.frame(txs$tx_id_version, txs$gene_id)

##################################################
##################################################

# AS: I use a different approach to make the tx2gene dataframe 
# AS: Imort rtracklayer to be able to read in a GTF as an R object 
library(rtracklayer)

# AS: Import annotation and combine transcript id with transcript version 
tx_gene_map <- rtracklayer::import(human_annotation) %>% 
  as_tibble() %>% # AS: Convert to tibble format as soon as possible 
  dplyr::select(transcript_id, transcript_version, transcript_biotype, gene_name, gene_id) %>% 
  unite(col = "transcript", c(transcript_id, transcript_version), sep = ".") %>% # AS: split transcript version and transcript ID 
  dplyr::rename(gene = gene_name) %>% 
  dplyr::select(transcript, gene) %>% 
  na.omit() %>% 
  distinct() 

##################################################
##################################################

# import salmon files into R using the transcript annotation we just made 
txi <- tximport(files, type = "salmon", tx2gene = tx_gene_map %>% as_data_frame())

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

# png("PCA_rlog_DDS.png")
# svg("PCA_rlog_dds.svg")
# dev.off()

# DESeq2 tutorial 
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

# AS: make a p-value histogram 
contrast_tibble <- res %>% 
  as_tibble(rownames = "geneID") %>% 
  dplyr::rename(LFC = 3) %>% 
  na.omit()
contrast_tibble %>% ggplot(aes(x = pvalue)) + geom_histogram()

# adjust p-values 
contrast_tibble$padj <- p.adjust(contrast_tibble$pvalue, method = "fdr")

# plot a histogram of adjusted p-values 
contrast_tibble %>% ggplot(aes(x = padj)) + geom_histogram()

# check number of differentially-expressed genes 
contrast_tibble %>% filter(padj < 0.05) %>% filter(abs(LFC) > 1) %>% dplyr::select(geneName)
clipr::write_clip(sigGenes)


# extract a results from DESEQ
# res_1 <- as_tibble(results(dds), rownames = "geneID") %>%
#  dplyr::rename(LFC = 3) %>%
#  na.omit()

# use biomart to add geneIDs to output
 #library(biomaRt)
#library("AnnotationDbi")
#library("org.Hs.eg.db")
#columns(org.Hs.eg.db)


# add geneIDs to contrast
#ens <- res_1$geneID
#symbols <- mapIds(org.Hs.eg.db, keys = ens,
#                  column = c('SYMBOL'), keytype = 'ENSEMBL')
#symbols <- symbols[!is.na(symbols)]
#symbols <- symbols[match(res_1$geneID, names(symbols))]
#res_1$geneName <- symbols
#res_1 <- res_1 %>% relocate(geneName, .after = "geneID")
#write_tsv(res_1, "2022-05-01_DEG_vs_NONDEG.txt.gz")

library(EnhancedVolcano)
EnhancedVolcano(contrast_tibble,
                lab = contrast_tibble$geneID,
                x = "LFC",
                y = "padj",
                title = 'DEG between degraded vs undegraded',
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


