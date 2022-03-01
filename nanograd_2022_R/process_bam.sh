#!/bin/bash

# some test data
# export bam="/g/data/lf10/as7425/2020-11_mouseBrain/analysis/2021-10-03_analyseEnsemblMappingFreezepoint/mergedH.bam"; export twd="/scratch/lf10/as7425/test_nanograd3"

# process a bam file to get genome coverage
export bam=$1
export twd=$2

# convert the bam to bed (smoothes indels)
bedtools bamtobed -split -i ${bam} > ${twd}/remove_indel.bed

# create a genome file for bedtools
samtools view -H ${bam} | grep @SQ | sed 's/@SQ\tSN:\|LN://g' >  ${twd}/genome.txt

# compute genome coverage
bedtools genomecov -i ${twd}/remove_indel.bed -g ${twd}/genome.txt -dz > ${twd}/coverage.txt
