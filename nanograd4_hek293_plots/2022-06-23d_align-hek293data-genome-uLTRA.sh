#!/bin/bash

# WARNING: Script is still being debugged

# written by Aditya Sethi on 2022-06-23
# aim: Align all hek293 data to the human genome using uLTRA

##################################################
##################################################

# configure environment

# loads conda in which the ultra enivronment is stored
source ~/.bashrc
conda activate ultra # activate the ultra environment

# load modules
module load R parallel

# define R library directory
export R_LIBS=/g/data/lf10/as7425/apps/Rlibs

# configure the timeformat for 'time'
export TIMEFORMAT=%R

##################################################

# configure run-specific variables

export ultraIndex="/g/data/lf10/as7425/genomes/human_genome/uLTRA_index_Homo_sapiens.GRCh38.dna_sm.primary_assembly/"
export genome="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
export annotation="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"
export allReads="/g/data/lf10/as7425/nanograd/data/2022-06-23_hek293-nanograd4-analysis"
export allAlign="/g/data/lf10/as7425/nanograd/analysis/2022-06-23_hek293-nanograd4-analysis/genome_all_alignments"; mkdir -p ${allAlign}
export primaryAlign="/g/data/lf10/as7425/nanograd/analysis/2022-06-23_hek293-nanograd4-analysis/genome_primary_alignments"; mkdir -p ${primaryAlign} 2>/dev/null

##################################################

# make uLTRA index for human
# we don't need to do this because we used Bhavika'a uLTRA index
# which we copied to "/g/data/lf10/as7425/genomes/human_genome/uLTRA_index_Homo_sapiens.GRCh38.dna_sm.primary_assembly/"

  # echo "$(date) .... starting indexing for ultra"
  # uLTRA index "${genome}" "${annotation}" "${ultraIndex}" || echo "$(date) ..... ultra indexing failed"
  # echo "$(date) .... done indexing for ultra"

# unzip all fastq
# gunzip "${allReads}/*"

# write a fuction to run uLTRA for each sample
function runUltra(){

  # get the bam name
  local reads="${1}"

  # get the sample name, omitting the extension '.fastq.gz'
  local sample=$(basename ${reads} .fastq.gz)

  # unzip the reads, only for running ultra
  echo "$(date) .... unzipping ${sample}"
  zcat ${reads} > ${allAlign}/tmp.${sample}.fastq

  # make the output in a subfolder because ultra temp files seem not to use the prefix?
  mkdir "${allAlign}/tmp_${sample}/"

  # run uLTRA
  echo "$(date) .... starting aligning for ${sample}"
  time uLTRA align "${genome}" "${allAlign}/tmp.${sample}.fastq" "${allAlign}/tmp_${sample}/" --index "${ultraIndex}" --ont --t 48 --prefix "${sample}" ||  echo "$(date) .... failed to align for ${sample}"
  echo "$(date) .... finished aligning for ${sample}"

  # remove the temporary unzipped fastq to save space
  # rm "${allAlign}/tmp.${sample}.fastq"

  # move everything out of the temp folder and remove it
  # mv "${allAlign}/tmp_${sample}/" "${allAlign}/" && rm -rf "${allAlign}/tmp_${sample}/"

  # convert the bam to sam and remove the original sam
  samtools view -b -@ 48 -u "<naming convention of the output>" | tee >(samtools sort -@ 48 -m 4G  -l 9 > "${allAlign}/all_genome_${sample}.bam") >(samtools view -b -@ 48 -u -F 2308 | samtools sort -@ 48 -m 4G -l 9 > "${primaryAlign}/primary_genome_${sample}.bam")

  # filter for primary alignments
  # samtools view -b -F 2308 -u -@ 48 ${i} | samtools sort -@ 48 -m 4G  -l 9 > ${primaryAlign}/${name}_primary.bam && samtools index ${primaryAlign}/${name}_primary.bam

}; export -f runUltra

rm -rf /g/data/lf10/as7425/nanograd/analysis/2022-06-23_hek293-nanograd4-analysis/genome_all_alignments/*
runUltra /scratch/lf10/as7425/test.fastq.gz


# use parallel to call 4 instances of uLTRA
# we will run the jobs on hugemem with 48 CPU cores and 1900 gb ram
# find ${allReads} -type f -name "*fastq.gz" | parallel -j 4 runUltra {} && echo "$(date) ... Ran uLTRA for all samples" || echo "$(date) ... failed to run uLTRA for one or more samples"
# remove the temporary files if the job finished succefully
# rm -rf ${allAlign}/tmp*

####################################################################################################
####################################################################################################

# runUltra /scratch/lf10/as7425/test.fastq.gz
# convert the sam to bam and filter by BAM FLAG
export allAlign="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/all_alignments"

# use loop structure to simplify
for i in ${allAlign}/*; do
name=$(basename $i .sam)
samtools view -b -F 2308 -u -@ 48 ${i} | samtools sort -@ 48 -m 4G -l 9 > ${primaryAlign}/${name}_primary.bam && samtools index ${primaryAlign}/${name}_primary.bam
done
