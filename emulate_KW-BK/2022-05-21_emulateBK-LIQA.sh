#!/bin/bash

# written by Bhavika Kumar and AJ Sethi on 2022-05-05

# aim: Identify differential splicing in degradation data using LIQA on genomic alignments
# samples: HEK293 polyA RNA (2 x WT RNA, 2 x degraded RNA)

############################################################
############################################################

# start

# set variables
# export annotation="/g/data/xc17/bk9031/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"
# export LIQA_index="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/LIQA/index"; mkdir -p ${LIQA_index} 2>/dev/null
# export primaryAlign="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/primaryAlignments"
# export isoform_expression="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/LIQA/isoform_expression"; mkdir -p ${isoform_expression} 2>/dev/null
#
# # set environment
# conda activate liqa
# module load R/3.6.1

############################################################

# make liqa index
# liqa -task refgene -ref ${annotation} -format gtf -out ${LIQA_index}/GRCh38.104.refgene

# works succesfully

############################################################

# set variables
conda activate liqa

export LIQA_index="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/LIQA/index"; mkdir -p ${LIQA_index} 2>/dev/null
export primaryAlign="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/primaryAlignments"
export isoform_expression="/scratch/lf10/as7425/2021_HEK293-degradation-first4-AR/LIQA/isoform_expression"; mkdir -p ${isoform_expression} 2>/dev/null
export annotation="/g/data/xc17/bk9031/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"

# loop for samples and start LIQA
IsoQuant(){
   local i=$1
   sample=$(basename $i .bam)
   echo "$(date) .... starting LIQA for ${sample}"
   time liqa -task quantify -refgene "${LIQA_index}/GRCh38.104.refgene" -bam "${i}" -out ${isoform_expression}/${sample}.txt -max_distance 10 -f_weight 1 ||  echo "$(date) .... failed to run LIQA for ${sample}"
   echo "$(date) .... finished LIQA for ${sample}"
}; export -f IsoQuant

 for i in ${primaryAlign}/*bam; do IsoQuant ${i} & done
 wait && echo "done for all"
