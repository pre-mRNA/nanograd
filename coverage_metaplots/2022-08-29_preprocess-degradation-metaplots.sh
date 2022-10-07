#!/bin/bash

# Aim: Generate data to produce coverage metaplots for Agin's degradation data
# Written by AJ Sethi on 2022-09-29
# Last modified on 2022-09-29

##############################

# assemble all 6 degradation primary alignments
export data="/g/data/lf10/as7425/nanograd/data/2022-08-29_degradation-first6-transcriptomePrimaryAlignments"
mkdir ${data} 2>/dev/null

# get the first 4 samples
cp /g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/transcriptome_alignments_fromAJ/*bam ${data}
cp /g/data/lf10/as7425/nanograd/analysis/2022-08-22_basecall-new-samples/primary_transcriptome_alignments/*bam ${data}

# sort data
mkdir ${data}/sorted/
module load samtools
for i in ${data}/*bam; do samtools sort -@ 48 ${i} > "${data}/sorted/${i##*/}.bam" & done

##############################

# make working directory
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-08-29_degradation-first6-metaplots"
mkdir ${wd} 2>/dev/null

# recall data dir
export data="/g/data/lf10/as7425/nanograd/data/2022-08-29_degradation-first6-transcriptomePrimaryAlignments/sorted/"

# load samtools
module load bedtools

##############################

# use bedtools genomecov to get coverage
for i in ${data}/*bam; do
  name=$(basename $i .bam)
  bedtools genomecov -ibam ${i} -dz -split | awk '{ print $1"\t"$2"\t"$2+1"\t.\t"$3"\t+"}' > "${wd}/${name}_coverage.txt" || echo "failed for ${name}" &
done
wait && echo "done for all samples"

for i in ${wd}/*txt; do sed -i '1i chr\tstart\tend\tname\tscore\tstrand' ${i}; done

##############################

# annotate samples

# make working directory
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-08-29_degradation-first6-metaplots"
mkdir ${wd}/annotated/
export annotated="${wd}/annotated/"

# use txannotate to annotate sites
module load R

# link GTF
export annotation="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"
export data="/g/data/lf10/as7425/nanograd/analysis/2022-08-29_degradation-first6-metaplots/primary_mild_degradataion_rep2.fastq.bam_coverage.txt"

for i in ${wd}/*txt; do
time Rscript /home/150/as7425/txannotate/annotate.R "${i}" "${annotation}" "${annotated}/annotated_${i##*/}" || echo "failed for ${i}" & done; wait && echo "done for all"

##############################

# filter out key data, i.e. rel_pos, score, transcript_biptype, transcript_length

# make working directory
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-08-29_degradation-first6-metaplots"
export annotated="${wd}/annotated/"
export trimmed="${wd}/trimmed"
mkdir ${trimmed}

for i in ${annotated}/*txt; do cat $i | cut -f5,7,10,17 | gzip -c > "${trimmed}/trimmed_${i##*/}.gz" & done
