#!/bin/bash
#PBS -P xc17
#PBS -l walltime=48:00:00
#PBS -l mem=190GB
#PBS -l ncpus=48
#PBS -q normal
#PBS -M bhavika.kumar@anu.edu.au
#PBS -l storage=scratch/xc17+gdata/xc17

# written by BK on 2022-09-23
# aim: script for nanopolish

#############################################################
#############################################################



# path to data
export data="/g/data/xc17/degradation_project/Mg_degraded/data"

# path to reference genome
export genome="/g/data/xc17/as7425/sharing/genomes/human_genome/transcriptome/GRCh38_codingPlusNoncoding_noPsuedo.fa"

# defining output path
export outputDir="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/nanopolish"
mkdir -p ${outputDir}

# nanopolish path
export PATH="${PATH}:/g/data/xc17/as7425/sharing/apps/nanopolish"

##################################################
##################################################

# for sample wt_1
sum_wt_1="/g/data/xc17/degradation_project/Mg_degraded/data/undegraded_hek293_pass1/sequencing_summary_FAQ86281_15c37cc7.txt"
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

##################################################
##################################################

# write a function to run nanopolish for each library
runNano() {

  # get positional arguments
  local seq_sum=$1
  local fast5=$2
  local fastq=$3
  local alignment=$4


  # get sample name
  local sampleLong=${alignment##*/}
  local sample=${sampleLong%.*}

  # make a nanopolish index
  echo "$(date)....starting indexing for nanopolish"
  nanopolish index -s ${seq_sum} -d ${fast5} ${fastq} || echo "$(date)....nanopolish indexing failed"
  echo "$(date)....done indexing for nanopolish"

  # align events
  echo "$(date)....starting alignment for nanopolish"
  nanopolish eventalign --reads ${fastq} --bam ${alignment} --genome ${genome} -v -t 48 --scale-events --print-read-names --samples > "${outputDir}/${sample}_events.tsv" || echo "$(date)....nanopolish alignment failed"
  echo "$(date)....done alignment for nanopolish"
}

# export function
export -f runNano

##################################################
##################################################

# test function for wt_rep1
time runNano "${sum_wt_1}" "${fast5_wt_1}" "${fastq_wt_1}" "${bam_wt_1}" || echo "failed for wt 1"
time runNano "${sum_wt_2}" "${fast5_wt_2}" "${fastq_wt_2}" "${bam_wt_2}" || echo "failed for wt 2"
time runNano "${sum_mild_deg1}" "${fast5_mild_deg1}" "${fastq_mild_deg1}" "${bam_mild_deg1}" || echo "failed for mild 1"
time runNano "${sum_mild_deg2}" "${fast5_mild_deg2}" "${fastq_mild_deg2}" "${bam_mild_deg2}" || echo "failed for mild 2"
time runNano "${sum_deg1}" "${fast5_deg1}" "${fastq_deg1}" "${bam_deg1}" || echo "failed for deg 1"
time runNano "${sum_deg2}" "${fast5_deg2}" "${fastq_deg2}" "${bam_deg2}" || echo "failed for deg 2"
