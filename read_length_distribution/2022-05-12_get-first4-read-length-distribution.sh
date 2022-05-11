#!/bin/bash

# written by AJ Sethi on 2022-05-05
# Get empirical read length distribution for uLTRA alignments of first 4 degradation libraries

# use tutorial at https://www.biostars.org/p/65216/

####################################################################################################
####################################################################################################

# set directories
align_dir="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/primaryAlignments"

# make output directory
analysis="/g/data/lf10/as7425/nanograd/analysis/2022-05-12_read-length-distribution"; mkdir -p ${analysis}

# configure environments
module load samtools

# loop and calculate lengths
for i in ${align_dir}/*bam; do
  export sample=$(basename $i .primary.bam)
  echo "$(date) .... starting read count for ${sample}"

    samtools view $i | awk '{print length($10)"\t"ENVIRON["sample"]}' > ${analysis}/${sample}_alignment_length.txt 

  echo "$(date) .... finished aligning for ${sample}"
done; wait && echo "done calculating lengths for all"
