#!/bin/bash

# written by AJ Sethi on 2022-05-11
# aim: Preprocess metaplot studies and calculate metacoverage

############################################################
############################################################

# paths
export genome="/g/data/xc17/degradation_project/Reference/GRCh38_codingPlusNoncoding_noPsuedo.fa"
export analysis="/g/data/lf10/as7425/nanograd/analysis/literature_metaplots"
export priAlign="${analysis}/primary_alignments"; mkdir -p ${priAlign}
export metaInput="${analysis}/metaplot_input"; mkdir -p ${metaInput}

# environment
module load samtools

cd ${priAlign} && time samtools depth -H *bam | gzip -c > ${metaInput}/2022-11-05_samtools-depth.txt.gz
