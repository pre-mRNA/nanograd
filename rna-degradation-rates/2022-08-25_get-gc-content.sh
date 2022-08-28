#!/bin/bash

# get GC content for all transcripts in the human genome

# install seqkit
cd /g/data/lf10/as7425/apps
wget https://github.com/shenwei356/seqkit/releases/download/v2.3.0/seqkit_linux_amd64.tar.gz
tar -xvzf seqkit_linux_amd64.tar.gz
export PATH="${PATH}:/g/data/lf10/as7425/apps/"

export genome="/g/data/lf10/as7425/genomes/human_genome/transcriptome/GRCh38_codingPlusNoncoding_noPsuedo.fa"
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-08-25_human-gc-content"; mkdir -p ${wd} 2>/dev/null

seqkit fx2tab --name --gc -G -B -n -i ${genome} | gzip -c > ${wd}/GRCh38_gc.txt.gz
