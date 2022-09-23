#!/bin/bash

# written by BK on 2022-09-23
# aim: script for nanopolish 

#############################################################
#############################################################

# getting data
export data="/g/data/xc17/degradation_project/Mg_degraded/data"

#  data preprocessing
running_nanopolish() {
  local seq_sum=$1
  local fast5=$2
  local fastq=$3
  local alignment=$4
  local genome=$5
  echo "$(data)....starting indexing for nanopolish"
  nanopolish index -s ${seq_sum} -d ${fast5} ${fastq} || echo "$(date)....nanopolish indexing failed"
  echo "$(date)....done indexing for nanopolish" 
  
  echo "$(data)....starting alignment for nanopolish"
  nanopolish eventalign --reads ${fastq} --bam ${alignment} --genome ${genome} || echo "$(date)....nanopolish alignment failed"
  echo "$(date)....done alignment for nanopolish"
}; export -f running_nanopolish


"/g/data/xc17/degradation_project/Mg_degraded/data/5mM_MgCl_degrdation_pass2 $path2";
"/g/data/xc17/degradation_project/Mg_degraded/data/MgCl_degrdaded_mild_pass1
"/g/data/xc17/degradation_project/Mg_degraded/data/MgCl_degrdaded_mild_pass2
"/g/data/xc17/degradation_project/Mg_degraded/data/undegraded_hek293_pass1
"/g/data/xc17/degradation_project/Mg_degraded/data/undegraded_hek293_pass2

"/g/data/xc17/degradation_project/Mg_degraded/data/undegraded_hek293_pass2/sequencing_summary_FAP73818_93b719e9.txt $seq_sum";
"/g/data/xc17/degradation_project/Mg_degraded/data/5mM_MgCl_degrdation_pass1 $fast5";



running_nanopolish “seq_sum" “fast5" "fastq" "alignment" "genome"

