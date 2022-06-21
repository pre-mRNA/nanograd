#!/bin/bash

# data is currently distributed in /g/data/xc17/degradation_project/Mg_degraded/data
all_data="/g/data/xc17/degradation_project/Mg_degraded/data"

# sequencing summary data to be collected in degradation project
ss_data="/g/data/lf10/as7425/nanograd/data/sequencing_summary_degradation"

# copy all sequencing summary files from all_data to ss_data
find ${all_data} -name "*sequencing_summary*" -exec cp "{}" ${ss_data} \;

# run nanograd for degradation data, retaining temporary outputs
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

# working directory is /g/data/lf10/as7425/nanograd/analysis/2022-04-05_test-nanograd4-decay
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-04-05_test-nanograd4-decay"
cd $wd
