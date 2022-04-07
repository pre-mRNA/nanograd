#!/bin/bash

# run nanograd4 for Agin's degraded samples

# written by AJ Sethi on 2022-04-05

# define variables
export data="/g/data/xc17/degradation_project/Mg_degraded/basecalled/"
export nanograd4="/home/150/as7425/nanograd/nanograd4/nanograd4.sh"
export annotation="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"

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
module load R

# find bams in data
for i in $(find $data -name "*bam" ); do runNano $i & done; wait && echo "done for all"
