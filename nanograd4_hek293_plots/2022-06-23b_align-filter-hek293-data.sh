#!/bin/bash

# Align and filter published direct RNA sequencing data from HEK293 prior to running nanograd4
# Written by A.J. Sethi on 2022-05-23

# refer to 'HEK293 Studies BK AS on Nanograd ANU OneDrive account'
# 2022-06-23 at https://anu365-my.sharepoint.com/:x:/r/personal/u6081208_anu_edu_au/_layouts/15/Doc.aspx?sourcedoc=%7BCFA6BAB3-538E-4999-869F-523EB1C6BDFD%7D&file=HEK293%20studies%20BK%20AS%20.xlsx&action=default&mobileredirect=true

##################################################
##################################################

module load samtools/1.10
module load minimap2

# set up data and analysis directories
export seq_data="/g/data/lf10/as7425/nanograd/data/2022-06-23_hek293-nanograd4-analysis"
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-06-23_hek293-nanograd4-analysis"; mkdir -p ${wd}

# set up directories for all alignments and primary alignments
export all_align="${wd}/all_alignments"; mkdir -p ${all_align}
export primary_align="${wd}/primary_alignments"; mkdir -p ${primary_align}

# set up reference for alignment
export genome="/g/data/lf10/as7425/genomes/human_genome/transcriptome/GRCh38_codingPlusNoncoding_noPsuedo.fa"

# index the genome
minimap2 -ax map-ont -k14 -d "${genome%.*}_mm2_index.mmi" "${genome}" || echo "cannot build index"

# export the index
export index="/g/data/lf10/as7425/genomes/human_genome/transcriptome/GRCh38_codingPlusNoncoding_noPsuedo_mm2_index.mmi"

# write a function to align all reads to the reference tranascriptome and filter out primary reads
# argument 1 is the path to the fastq
function align(){

  local reads=${1}
  local sample_name=$(basename $reads .fastq.gz)
  # echo $sample_name

  # align using minimap2
  minimap2 -t 48 -2 -I 16G -ax map-ont -k14 "${index}" "${reads}" | samtools view -@ 48 -m 3.8G -b > "${all_align}/all_${sample_name}.bam" || echo "Cannot align ${reads}"

  # sort for primary alignents
  samtools view -@ 48 -m 3.8G -F 2324 -b "${all_align}/all_${sample_name}.bam" > "${primary_align}/primary_${sample_name}.bam" || echo "cannot filter ${reads}"

}; export -f align

for i in ${seq_data}/*fastq*; do align $i || echo "cannot align ${i}"; done

# gather some data for onedrive
cd ${seq_data} && for i in *fastq.gz; do printf "${i}\t$(echo $(zcat $i|wc -l)/4|bc)\n" >> ${seq_data}/2022-06-23_hek293-reads-data.txt; done
cd ${all_align} && for i in *bam; do printf "${i}\t$(samtools view -@ 48 $i | wc -l)\n" >> ${all_align}/2022-06-23_hek293-all-alignments-data.txt; done
cd ${primary_align} && for i in *bam; do printf "${i}\t$(samtools view -@ 48 $i | wc -l)\n" >> ${primary_align}/2022-06-23_hek293-primary-alignments-data.txt; done
