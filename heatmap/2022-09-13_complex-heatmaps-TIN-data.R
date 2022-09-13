#!/usr/bin/env Rscript 

# written by AJ and BK on 2022-09-13

# aim: To visualise RNA degradation. Making complex heat maps from TIN data.

###########################################################

# load libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# importing dataset

data_filt12 <- read_tsv("2022-08-24_degradationFirst6-poolTin-filt12.txt.gz")
data_filt5 <- read_tsv("2022-08-24_degradationFirst6-poolTin-filt5.txt.gz", col_types="fdd")
data_nocutoff <- read_tsv("2022-08-24_degradationFirst6-poolTin-noCutoff.txt.gz")

# making a single heatmap of data_filt5

matrix <- as.matrix(data_filt5[2:ncol(data_filt5)])
rownames(matrix) <- data_filt5$transcript_id
Heatmap(matrix)

# Apply clustering

Heatmap(matrix, name="matrix", cluster_rows=F)

Heatmap(matrix, name="matrix", column_dend_height = unit(2,"cm"), row_dend_width = unit(2, "cm"))  

# Applying Hierarchical clustering
Heatmap(matrix, name="matrix", clustering_distance_rows = "pearson", column_title = "Pearson")
Heatmap(matrix, name="matrix", clustering_distance_rows = "kendall", column_title = "kendall")

Heatmap(matrix, name="matrix", clustering_method_rows = "single")

# Splitting by k-means clustering
Heatmap(matrix, name="matrix", row_km=4)

# Adjusting colors in heatmap
col_matrix <- colorRamp2(c(0,1), c("white", "blue"))
col_matrix(seq(0, 1, 0.01))
Heatmap(matrix, name="matrix", row_km=4, col=col_matrix, show_row_names = F)

