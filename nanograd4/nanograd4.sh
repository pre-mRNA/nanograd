#!/bin/bash

# arg1 bam (transcriptome alignment)
# arg2 gtf
# arg3 output

##########
samtools view -F 256 ${bam} | cut -f1,3,6,10
