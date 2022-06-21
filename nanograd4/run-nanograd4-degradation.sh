#!/bin/bash

# run nanograd4 for Agin's degraded samples

# written by AJ Sethi on 2022-04-05

# define variables
export data="/g/data/xc17/degradation_project/Mg_degraded/basecalled/"
export nanograd4="/home/150/as7425/nanograd/nanograd4/nanograd4.sh"
export annotation="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"
# export annotation="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens_transcriptsOnly_GRCh38.104.chr.gtf" # transcripts only

# make working directory
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-04-05_test-nanograd4-decay"

# write a function to run nanograd for a given bam
function runNano(){
  local bam=$1

  # get the bamname
  local bamLong=${bam##*/}
  local name=${bamLong%.*}

  time bash ${nanograd4} ${bam} ${annotation} ${wd}/${name}_nanograd4_out.txt
}; export -f runNano

# load R
module load R samtools

# find bams in data
for i in $(find $data -name "*bam" ); do runNano $i; done; wait && echo "done for all"

# for benchmarking, count reads in library
export data="/g/data/xc17/degradation_project/Mg_degraded/basecalled/"
for i in $(find $data -name "*bam" ); do printf "$i\t$(samtools view $i | wc -l)\n"; done; wait && echo "done for all"

# for benchmarking, sample 60% of reads from wt_rep1 and run nanograd on that in a temporary directory
export wd="/g/data/lf10/as7425/nanograd/analysis"
mkdir ${wd}
export full_bam="/g/data/xc17/degradation_project/Mg_degraded/basecalled/undegraded_hek293_pass1/all.undegraded_hek293_pass1.fastq.gz.sorted.bam"
samtools view -b -s 0.6 ${full_bam} > ${wd}/wt_rep1_60percent.bam
export bam="${wd}/wt_rep1_60percent.bam"
samtools view ${bam} | wc -l # count reads
export nanograd4="/home/150/as7425/nanograd/nanograd4/nanograd4.sh"
export annotation="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens_transcriptsOnly_GRCh38.104.chr.gtf" # transcripts only
time bash ${nanograd4} ${bam} ${annotation} ${wd}/${name}_nanograd4_out.txt

# for benchmarking, sample 60% of reads from wt_rep1 and run nanograd on that in a temporary directory
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-06-22_benchmark-nanograd4-runtime/"
mkdir ${wd}
export full_bam="/g/data/xc17/degradation_project/Mg_degraded/basecalled/undegraded_hek293_pass1/all.undegraded_hek293_pass1.fastq.gz.sorted.bam"
samtools view -@ 48 -b -s 0.8 ${full_bam} > ${wd}/wt_rep1_80percent.bam
export bam="${wd}/wt_rep1_80percent.bam"
samtools view -@ 48 ${bam} | wc -l # count reads
export nanograd4="/home/150/as7425/nanograd/nanograd4/nanograd4.sh"
export annotation="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens_transcriptsOnly_GRCh38.104.chr.gtf" # transcripts only
time bash ${nanograd4} ${bam} ${annotation} ${wd}/wt_rep1_80percent_nanograd4_out.txt
