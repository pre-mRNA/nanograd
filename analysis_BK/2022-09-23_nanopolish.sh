#!/bin/bash

# written by BK on 2022-09-23
# aim: script for nanopolish

#############################################################
#############################################################

# getting data
export data="/g/data/xc17/degradation_project/Mg_degraded/data"

# installing nanopolish
# export install_nanopolish="/g/data/xc17/bk9031/apps"

# cd ${install_nanopolish}

# git clone --recursive https://github.com/jts/nanopolish.git
# cd nanopolish
# make

#  data preprocessing
runNano() {
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
}; export -f runNano

# defining path of genome
export genome="/g/data/xc17/bk9031/genomes/human_genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"

# example command
runNano ${sum_wt_1} ${fast5_wt_1} ${fastq_wtt_1} ${bam_wt_1} ${genome}


# genome
export genome="/g/data/xc17/bk9031/genomes/human_genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"

# for sample wt_1
sum_wt_1= "/g/data/xc17/degradation_project/Mg_degraded/data/undegraded_hek293_pass1/sequencing_summary_FAQ86281_15c37cc7.txt"
fast5_wt_1="/g/data/xc17/degradation_project/Mg_degraded/data/undegraded_hek293_pass1"
fastq_wt_1="/g/data/xc17/degradation_project/Mg_degraded/basecalled/undegraded_hek293_pass1/all.undegraded_hek293_pass1.fastq.gz"
bam_wt_1="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/primaryAlignments/undegraded_hek293_pass1_primary.bam"

# for sample wt_2
sum_wt_2="/g/data/xc17/degradation_project/Mg_degraded/data/undegraded_hek293_pass2/sequencing_summary_FAP73818_93b719e9.txt"
fast5_wt_2="/g/data/xc17/degradation_project/Mg_degraded/data/undegraded_hek293_pass2";
fastq_wt_2="/g/data/xc17/degradation_project/Mg_degraded/basecalled/undegraded_hek293_pass2/all.undegraded_hek293_pass2.fastq.gz"
bam_wt_2="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/primaryAlignments/undegraded_hek293_pass2_primary.bam"

# for MgCl mild degradation sample 1
sum_mild_deg1="/g/data/xc17/as7425/sharing/2022-09-26_nanograd-mildDegradation-transcriptomePrimaryAlignments/basecall_mild_degraded_rep1/sequencing_summary.txt"
fast5_mild_deg1="/g/data/xc17/degradation_project/Mg_degraded/data/MgCl_degrdaded_mild_pass1/fast5"
fastq_mild_deg1="/g/data/xc17/as7425/sharing/2022-09-26_nanograd-mildDegradation-transcriptomePrimaryAlignments/fastq/mild_degradataion_rep1.fastq"
bam_mild_deg1="/g/data/xc17/as7425/sharing/2022-09-26_nanograd-mildDegradation-transcriptomePrimaryAlignments/primary_mild_degradataion_rep1.fastq.bam"

# for MgCl mild degradation sample 2
sum_mild_deg2="/g/data/xc17/as7425/sharing/2022-09-26_nanograd-mildDegradation-transcriptomePrimaryAlignments/basecall_mild_degraded_rep2sequencing_summary.txt"
fast5_mild_deg2="/g/data/xc17/degradation_project/Mg_degraded/data/MgCl_degrdaded_mild_pass2/fast5"
fastq_mild_deg2="/g/data/xc17/as7425/sharing/2022-09-26_nanograd-mildDegradation-transcriptomePrimaryAlignments/fastq/mild_degradataion_rep2.fastq"
bam_mild_deg2="/g/data/xc17/as7425/sharing/2022-09-26_nanograd-mildDegradation-transcriptomePrimaryAlignments/primary_mild_degradataion_rep2.fastq.bam"

# for 5mM MgCl degraded sample 1
sum_deg1="/g/data/xc17/degradation_project/Mg_degraded/data/5mM_MgCl_degrdation_pass1/sequencing_summary_FAR08510_cd74bc1a.txt"
fast5_deg1="/g/data/xc17/degradation_project/Mg_degraded/data/5mM_MgCl_degrdation_pass1    "
fastq_deg1="/g/data/xc17/bk9031/2022_nanograd_bk/data/2021_HEK293-degradation-first4-AR/5mM_MgCl_degrdation_pass1/all.5mM_MgCl_degrdation_pass1.fastq.gz"
bam_deg1="/g/data/xc17/bk9031/2022_nanograd_bk/data/2021_HEK293-degradation-first4-AR/5mM_MgCl_degrdation_pass1/all.5mM_MgCl_degrdation_pass1.fastq.gz.sorted.bam"

# for 5mM MgCl degraded sample 2
sum_deg2="/g/data/xc17/degradation_project/Mg_degraded/data/5mM_MgCl_degrdation_pass2/sequencing_summary_FAP73846_c4e33c9d.txt"
fast5_deg2="/g/data/xc17/degradation_project/Mg_degraded/data/5mM_MgCl_degrdation_pass2        "
fastq_deg2="/g/data/xc17/bk9031/2022_nanograd_bk/data/2021_HEK293-degradation-first4-AR/5mM_MgCl_degrdation_pass2/all.5mM_MgCl_degrdation_pass2.fastq.gz"
bam_deg2="/g/data/xc17/bk9031/2022_nanograd_bk/data/2021_HEK293-degradation-first4-AR/5mM_MgCl_degrdation_pass2/all.5mM_MgCl_degrdation_pass2.fastq.gz.sorted.bam"
