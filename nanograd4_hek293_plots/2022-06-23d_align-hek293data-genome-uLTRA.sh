#!/bin/bash

# written by Aditya Sethi on 2022-06-23
# aim: Align all hek293 data to the human genome using uLTRA

##################################################
##################################################

# configure environment
source ~/.bashrc # loads conda in which the ultra enivronment is stored
conda activate ultra
export R_LIBS=/g/data/lf10/as7425/apps/Rlibs
export TIMEFORMAT=%R # to print the time in seconds

##################################################

# configure run-specific variables

export ultraIndex="/g/data/lf10/as7425/genomes/human_genome/uLTRA_index_Homo_sapiens.GRCh38.dna_sm.primary_assembly/"
export genome="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
export annotation="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"
export allAlign="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/all_alignments"
export allReads="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/fastq"


##################################################

# make uLTRA index for human
echo "$(date) .... starting indexing for ultra"
uLTRA index "${genome}" "${annotation}" "${ultraIndex}" || echo "$(date) ..... ultra indexing failed"
echo "$(date) .... done indexing for ultra"

# unzip all fastq
gunzip "${allReads}/*"

# loop over samples and start uLTRA
<<<<<<< HEAD
# note: for 48 cores, > 190 gb is used.
=======
# note: for 48 cores, > 190 gb is used.
>>>>>>> bcf3dc3d3c6b43012c4ce3d74ec43c235482f46d
for i in ${allReads}/*; do
   sample=$(basename $i .fastq | cut -c 5-)
   echo "$(date) .... starting aligning for ${sample}"
   time uLTRA align "${genome}" "${i}" "${allAlign}" --index "${ultraIndex}" --ont --t 48 --prefix "${sample}" ||  echo "$(date) .... failed to align for ${sample}"
   echo "$(date) .... finished aligning for ${sample}"
 done

# convert the sam to bam and filter by BAM FLAG
export allAlign="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/all_alignments"
export primaryAlign="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/primaryAlignments"; mkdir -p ${primaryAlign} 2>/dev/null

# use loop structure to simplify
for i in ${allAlign}/*; do
name=$(basename $i .sam)
samtools view -b -F 2308 -u -@ 48 ${i} | samtools sort -@ 48 -m 4G  -l 9 > ${primaryAlign}/${name}_primary.bam && samtools index ${primaryAlign}/${name}_primary.bam
done
