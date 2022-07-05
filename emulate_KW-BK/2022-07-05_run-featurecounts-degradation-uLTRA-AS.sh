#!/bin/bash

# written by AJ Sethi on 2022-05-11
# aim: Quantify gene expression in uLTRA genomic alignments of Agin's degradation data using subread featurecounts
# rus on AJ's iMac

############################################################
############################################################

export annotation="/Users/AJlocal/localGadiData/0_a_tmp/Homo_sapiens.GRCh38.104.chr.gtf"
# brew install brewsci/bio/subread
cd /Users/AJlocal/localGadiData/2022-06-22_HEK293-degradation-first4-AR_liqa-genome-alignments_BK

featureCounts --primary -L -s 1 -a ${annotation} --extraAttributes "gene_name" `ls *bam` -o 2022-07-05_ARdegradation-genomic-featurecounts-AS.txt
