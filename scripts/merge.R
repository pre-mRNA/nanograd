# written by A.J. Sethi on 2020-12-13
# use tidyverse to merge cluster information with decay info and length info

# first arg: cluster info
# second arg: decay info
# third arg: length info
# fourth arg: output information

library(tidyverse)

# parse the arguments
args <- commandArgs()
clusterImportF <- args[7]
clusterDecayF <- args[8]
clusterLengthF <- args[9]
outFileF <- args[10]

# import the files
clusterImport <- read_tsv(clusterImportF, col_names = T)
clusterDecay <- read_delim(clusterDecayF, col_names = F, delim = " ") %>% select(cluster = 1, decayCoefficient = 2)
clusterLength <- read_delim(clusterLengthF, col_names = F, delim = " ") %>% select(cluster = 1, clusterLength = 2)

# merge the files
firstMerge <- left_join(clusterDecay, clusterImport, by = "cluster")
secondMerge <- left_join(firstMerge, clusterLength, by = "cluster") %>% select(chromosome, start, end, cluster, count, strand, clusterLength, decayCoefficient)

# write output
write_tsv(secondMerge , outFileF, col_names = T)
