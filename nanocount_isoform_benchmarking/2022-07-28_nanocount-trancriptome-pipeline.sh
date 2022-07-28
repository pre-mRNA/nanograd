#!/bin/bash

# written by AJ Sethi on 2022-05-11
# aim: Write a pipeline to do efficient transcriptome alignment using nanocount + minimap2

############################################################

# housekeeping

conda activate NanoCount

# reference transcriptome
export transcripts="/g/data/lf10/as7425/genomes/human_genome/transcriptome/GRCh38_codingPlusNoncoding_noPsuedo.fa"

# minimap2 to path
export PATH="${PATH}:/g/data/lf10/as7425/apps/minimap2"
module load samtools

# read directory
export data="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/fastq"

# make an output directory
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-07-28_align-transcriptome-nanocount"


############################################################

# function to align using minimap2 and select alignment using nanocount
function alignNano(){

  # import reads
  local reads=${1}

  # get sample name
  sample=$(basename $reads .fastq | cut -c 5-)

  # # say that we're aligning
  # echo "$(date) .... aligning for ${sample}"
  #
  # # first, align the reads to the reference transcriptome using minimap2
  # minimap2 -t 48 -ax map-ont -N 30 "${transcripts}" "${reads}" | samtools view -b -@ 48 | samtools sort -@ 48 > "${wd}/allAlignments_sorted_${sample}.bam" || echo "alignment failed for ${sample}"
  # samtools index "${wd}/allAlignments_sorted_${sample}.bam"
  #
  # # say that we're aligning
  # echo "$(date) .... nanocount for ${sample}"

  # run nanocount
  NanoCount -i "${wd}/allAlignments_sorted_${sample}.bam" -o "${wd}/isoformCounts_${sample}.txt" -b "${wd}/filtAlignments_sorted_${sample}.bam" || echo "nanocount failed for ${sample}"


}; export -f alignNano

# run the pipeline for each fastq
for i in ${data}/*fastq; do alignNano "${i}" & done; wait && echo "done for all"
exit

# count gene expression
for i in filt*; do samtools view -F 2308 $i | cut -f2 | sort | uniq -c > "./counts_${i##*/}.txt" & done; wait && echo "done for all"
