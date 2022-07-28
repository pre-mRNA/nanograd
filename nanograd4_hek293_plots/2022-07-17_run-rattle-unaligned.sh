#!/bin/bash

# written by Aditya Sethi on 2022-07-17
# ran rattle on unaligned reads to figure out their identity
# start from transcriptome alignments

##################################################
##################################################

# run rattle on all hek293 data

reads="/g/data/lf10/as7425/nanograd/data/2022-06-23_hek293-nanograd4-analysis"

# make analysis directory
mkdir -p /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/

# merge all the reads
zcat ${reads}/*gz > /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/all_reads.fastq

# subset the reads
cd /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/
cat all_reads.fastq | head -n 2000000 > 500k_read_subset.fastq
cat all_reads.fastq | head -n 8000000 > 4m_read_subset.fastq
# run rattle

# download

# rattle step 1: cluster reads
/g/data/lf10/as7425/apps/RATTLE/rattle cluster -i /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/all_reads.fastq -t 96 --iso --rna -o /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/
# job ended after 24h without finishing

# rerun on the 500k read subset
mkdir /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/rattle-cluster-subset
/g/data/lf10/as7425/apps/RATTLE/rattle cluster -i /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/500k_read_subset.fastq -t 48 --iso --rna -o /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/rattle-cluster-subset

# RNA mode: true
# Reading fasta file...
# Done
# Gene clustering done
# 35788 gene clusters found
# finishes in 1:15

# rerun on the 2m read subset
# rerun on the 500k read subset
mkdir /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/rattle-cluster-subset-2m
/g/data/lf10/as7425/apps/RATTLE/rattle cluster -i /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/4m_read_subset.fastq -t 48 --iso --rna -o /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/rattle-cluster-subset-2m/

# 8h walltime for 2m reads

# RNA mode: true
# Reading fasta file...
# Done
# Gene clustering done
# 121755 gene clusters found

# step 2: cluster extraction

/g/data/lf10/as7425/apps/RATTLE/rattle extract_clusters -i /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/4m_read_subset.fastq -t 48 --iso --rna -o /g/data/lf10/as7425/nanograd/analysis/2022-07-17_run-rattle-hek293/rattle-cluster-subset-2m/
