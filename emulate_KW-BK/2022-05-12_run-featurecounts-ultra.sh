#!/bin/bash

# written by AJ Sethi on 2022-05-11
# aim: run subread featurecounts for Bhavika's uLTRA alignments

############################################################
############################################################

# count genes
export annotation="/g/data/xc17/bk9031/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"
export PATH="${PATH}:/g/data/lf10/as7425/apps/subread-2.0.1-Linux-x86_64/bin"

export wd="/g/data/lf10/as7425/nanograd/analysis/2022-05-12_uLTRA-featurecounts"
mkdir $wd
cd $wd

# run featurecounts
featureCounts --primary -L -T 56 -s 1 --extraAttributes "gene_biotype, gene_name" -a ${annotation} -o ${wd}/featureCounts_primary.txt `ls /g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/primaryAlignments/*.bam`

# zip results
cat ${wd}/featureCounts_primary.txt | gzip -c > ${wd}/2022-05-12_degradation-first4-featureCounts_primary.txt.gz
