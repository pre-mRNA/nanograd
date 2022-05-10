#!/bin/bash

# written by AJ Sethi on 2022-05-04

# emulate some of Bhavika's code, unzipping fastqs to see if it works

# original code snippet:

# run uLTRA
# conda activate ultra
# module load R
#
# export genome="/g/data/xc17/bk9031/genomes/human_genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
# export annotation="/g/data/xc17/bk9031/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"
# export allAlign="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/all_alignments"
# export allReads="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/fastq"
# export ultraIndex="/g/data/xc17/bk9031/genomes/human_genome/ultra_index"; mkdir -p ${ultraIndex}
#
# # make uLTRA index for human
# echo "$(date) .... starting indexing for ultra"
# uLTRA index "${genome}" "${annotation}" "${ultraIndex}" || echo "$(date) ..... ultra indexing failed"
# echo "$(date) .... done indexing for ultra"
#
# # loop over samples and start uLTRA
# for i in ${allReads}/*; do
#    sample=$(basename $i .fastq.gz | cut -c 5-)
#    echo "$(date) .... starting aligning for ${sample}"
#    time uLTRA align "${genome}" "${i}" "${allAlign}" --index "${ultraIndex}" --ont --t 48 --prefix "${sample}" ||  echo "$(date) .... failed to align for ${sample}"
#    echo "$(date) .... finished aligning for ${sample}"
#  done
#

### my code

# copy the fastqs to another directory and unzip them
my_reads="/scratch/lf10/as7425/test_ultra_nanograd/fastq"; mkdir -p ${my_reads}
bk_reads="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/fastq"
cp ${bk_reads}/* ${my_reads}
gunzip ${my_reads}/*

export genome="/g/data/xc17/bk9031/genomes/human_genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
export annotation="/g/data/xc17/bk9031/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"
export allAlign="/scratch/lf10/as7425/test_ultra_nanograd/alignments"; mkdir -p ${allAlign}
export allReads="/scratch/lf10/as7425/test_ultra_nanograd/fastq"
export ultraIndex="/g/data/xc17/bk9031/genomes/human_genome/ultra_index"; mkdir -p ${ultraIndex}

conda activate ultra
module load R

# loop over samples and start uLTRA
for i in ${allReads}/*; do
   sample=$(basename $i .fastq | cut -c 5-)
   echo "$(date) .... starting aligning for ${sample}"
   time uLTRA align "${genome}" "${i}" "${allAlign}" --index "${ultraIndex}" --ont --t 48 --prefix "${sample}" ||  echo "$(date) .... failed to align for ${sample}"
   echo "$(date) .... finished aligning for ${sample}"
 done
