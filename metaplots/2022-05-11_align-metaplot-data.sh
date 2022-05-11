#!/bin/bash

# written by AJ Sethi on 2022-05-11
# aim: Align metaplots studies to human transcriptome

############################################################
############################################################

# preprocessing

# data is in /g/data/xc17/bk9031/2022_nanograd_bk/data/2021_HEK293-degradation-first4-AR
export data="/g/data/lf10/as7425/nanograd/data/metaplot_data"

# make an output directory
export analysis="/g/data/lf10/as7425/nanograd/analysis/metaplot_data"
mkdir -p ${analysis} 2>/dev/null

# make a folder for mapping output
mkdir ${analysis}/all_alignments/

#######################################################

# run minimap2

conda activate ultra
module load R samtools

export genome="/g/data/xc17/degradation_project/Reference/GRCh38_codingPlusNoncoding_noPsuedo.fa"
export annotation="/g/data/xc17/bk9031/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"
export allReads="/g/data/lf10/as7425/nanograd/data/metaplot_data"
export analysis="/g/data/lf10/as7425/nanograd/analysis/literature_metaplots"
export allAlign="${analysis}/all_alignments"; mkdir -p ${allAlign}
export priAlign="${analysis}/primary_alignments"; mkdir -p ${priAlign}

# unzip all fastq
gunzip ${allReads}/*

# loop over samples and start uLTRA
for i in ${allReads}/*; do
   sample=$(basename $i .fastq)
   echo "$(date) .... starting aligning for ${sample}"
   #minimap2 -ax map-ont -k14 "${genome}" "${i}" -a -t 48 | samtools view -@ 48 -b > ${allAlign}/${sample}.bam ||  die "$(date) .... failed to align for ${sample}"
   samtools view -b -F 2308 -u -@ 48 ${allAlign}/${sample}.bam | samtools sort -@ 48 -m 4G  -l 9 > ${priAlign}/${sample}_primary.bam && samtools index ${priAlign}/${sample}_primary.bam &
   echo "$(date) .... finished aligning for ${sample}"
 done; wait && echo "done aliging for all"
