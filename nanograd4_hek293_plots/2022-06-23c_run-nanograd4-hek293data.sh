#!/bin/bash

# run nanograd4 for all the hek293 data
# written by AJ Sethi on 2022-04-05

##################################################
##################################################

# define variables
export nanograd4="/home/150/as7425/nanograd/nanograd4/nanograd4.sh"
export primary_alignments="/g/data/lf10/as7425/nanograd/analysis/2022-06-23_hek293-nanograd4-analysis/primary_transcriptome_alignments"
export nanograd_out="/g/data/lf10/as7425/nanograd/analysis/2022-06-23_hek293-nanograd4-analysis/nanograd4_out"; mkdir -p "${nanograd_out}"
export annotation="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"

##################################################

# configure environment
module load R samtools
export R_LIBS=/g/data/lf10/as7425/apps/Rlibs
export TIMEFORMAT=%R # to print the time in seconds

##################################################

# write a function to run nanograd for a given bam
# bam is supplied as the first positional argument
function runNano(){

  local primary_alignment="${1}"

  # get the bamname
  local sample_name=$(basename ${primary_alignment} .bam | cut -c 9-)

  # test if the sample name is correct
  # echo ${sample_name}

  # run nanograd
  time bash "${nanograd4}" "${primary_alignment}" "${annotation}" "${nanograd_out}/${sample_name}_nanograd4_out.txt" || echo "nanograd4 failed for ${sample_name}"
  echo "$(date) ... done for ${sample_name}"

}; export -f runNano

# completed succefully 2022-06-23 AS on first try! 

for i in ${primary_alignments}/*bam; do runNano "${i}" || echo "cannot run nanograd for ${i}" & done
