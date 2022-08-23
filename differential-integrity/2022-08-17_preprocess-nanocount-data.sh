#!/bin/bash

# Written by AJ Sethi on 2022-07-18
# download and preprocess the nanocount data for nanograd5

##################################################
##################################################

# download the sequencing data

# set wd
cd /g/data/lf10/as7425/nanograd/data/2022-08-18_nanocount-differentiation-data

# load aspera

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/003/ERR4352443/ERR4352443.fastq.gz -o ERR4352443_MinION_sequencing.fastq.gz || echo "could not download 1" &
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/001/ERR4352441/ERR4352441.fastq.gz -o ERR4352441_MinION_sequencing.fastq.gz || echo "could not download 2" &
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/002/ERR4352442/ERR4352442.fastq.gz -o ERR4352442_MinION_sequencing.fastq.gz || echo "could not download 3" &
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR436/009/ERR4368409/ERR4368409.fastq.gz -o ERR4368409_MinION_sequencing.fastq.gz || echo "could not download 4" &
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/004/ERR4352444/ERR4352444.fastq.gz -o ERR4352444_MinION_sequencing.fastq.gz || echo "could not download 5" &
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR436/000/ERR4368410/ERR4368410.fastq.gz -o ERR4368410_MinION_sequencing.fastq.gz || echo "could not download 6" &
wait && echo "done for all"

##################################################
##################################################

# Align and filter the sequencing data

module load samtools/1.10
module load minimap2

# set up data and analysis directories
export seq_data="/g/data/lf10/as7425/nanograd/data/2022-08-18_nanocount-differentiation-data"
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-08-17_nanocount-differential-integrity-analysis"; mkdir -p ${wd}

# set up directories for all alignments and primary alignments
export all_align="${wd}/all_transcriptome_alignments"; mkdir -p ${all_align}
export primary_align="${wd}/primary_transcriptome_alignments"; mkdir -p ${primary_align}

# set up reference for alignment
export genome="/g/data/lf10/as7425/genomes/human_genome/transcriptome/GRCh38_codingPlusNoncoding_noPsuedo.fa"

# index the genome
# note: this has already been done
# minimap2 -ax map-ont -k14 -d "${genome%.*}_mm2_index.mmi" "${genome}" || echo "cannot build index"

# export the index
export index="/g/data/lf10/as7425/genomes/human_genome/transcriptome/GRCh38_codingPlusNoncoding_noPsuedo_mm2_index.mmi"

# write a function to align all reads to the reference tranascriptome and filter out primary reads
# argument 1 is the path to the fastq
function align(){

  local reads=${1}
  local sample_name=$(basename $reads .fastq.gz)
  # echo $sample_name

  # align using minimap2 and convert to bam without filtering
  minimap2 -t 48 -2 -I 16G -ax map-ont -k14 "${index}" "${reads}" | samtools view -@ 48 -m 3.8G -b > "${all_align}/all_${sample_name}.bam" || echo "Cannot align ${reads}"

  # filter for primary alignents
  samtools view -@ 48 -m 3.8G -F 2324 -b "${all_align}/all_${sample_name}.bam" > "${primary_align}/primary_${sample_name}.bam" || echo "cannot filter ${reads}"

}; export -f align

for i in ${seq_data}/*fastq*; do align $i || echo "cannot align ${i}"; done

# gather some data for onedrive
cd ${seq_data} && for i in *fastq.gz; do printf "${i}\t$(echo $(zcat $i|wc -l)/4|bc)\n" >> ${seq_data}/2022-08-17_nanocount-differential-integrity-reads-data.txt; done
cd ${all_align} && for i in *bam; do printf "${i}\t$(samtools view -@ 48 $i | wc -l)\n" >> ${all_align}/2022-08-17_nanocount-differential-integrity-all-alignments-data.txt; done
cd ${primary_align} && for i in *bam; do printf "${i}\t$(samtools view -@ 48 $i | wc -l)\n" >> ${primary_align}/2022-08-17_nanocount-differential-integrity-primary-alignments-data.txt; done

##################################################
##################################################

# run nanograd5

# define variables
export nanograd5="/home/150/as7425/nanograd/nanograd5/nanograd5.sh"
export primary_alignments="/g/data/lf10/as7425/nanograd/analysis/2022-08-17_nanocount-differential-integrity-analysis/primary_transcriptome_alignments"
export nanograd_out="/g/data/lf10/as7425/nanograd/analysis/2022-08-17_nanocount-differential-integrity-analysis/nanograd5_out"; mkdir -p "${nanograd_out}"
export annotation="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"

##################################################

# configure environment
module load R samtools
export R_LIBS=/g/data/lf10/as7425/apps/Rlibs
export TIMEFORMAT=%R # to print the time in seconds

##################################################

# write a function to run nanograd for a given bam
# bam is supplied as the first positional argument
function runNano(){

  local primary_alignment="${1}"

  # get the bamname
  local sample_name=$(basename ${primary_alignment} .bam | cut -c 9-)

  # test if the sample name is correct
  # echo ${sample_name}

  # run nanograd
  time bash "${nanograd5}" "${primary_alignment}" "${annotation}" "${nanograd_out}/${sample_name}_nanograd5_out.txt" || echo "nanograd5 failed for ${sample_name}"
  echo "$(date) ... done for ${sample_name}"

}; export -f runNano

# completed succefully 2022-06-23 AS on first try!

for i in ${primary_alignments}/*bam; do runNano "${i}" || echo "cannot run nanograd for ${i}" & done
