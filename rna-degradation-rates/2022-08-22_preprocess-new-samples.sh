#!/bin/bash

# Aim: Basecall the two intermediate degradation HEK293 samples, and align to transcriptome and then run nanograd5
# Last modified on 2022-08-22 by AS

############################################################
############################################################

# setting up a working directory
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-08-22_basecall-new-samples"
mkdir -p ${wd} 2>/dev/null

############################################################

# basecalling the samples

# write a function for basecalling
basecall(){

  # get fast5 from first positional argument
  local fast5="${1}"

  # get fastq from second positional arguments
  local fastq="${2}"

  echo "running module guppy for ${fast5}"

  # guppy parameters tuned via https://hackmd.io/@Miles/S12SKP115
  # and via https://gist.github.com/sirselim/2ebe2807112fae93809aa18f096dbb94#guppy-basecalling-benchmarking-on-a-titan-rtx

  # /g/data/lf10/as7425/apps/ont-guppy/bin/guppy_basecaller --input_path "${fast5}" --save_path "${fastq}" \
  # --flowcell FLO-MIN106 --kit SQK-RNA002 --disable_qscore_filtering --reverse_sequence --u_substitution --disable_pings \
  # --records_per_fastq 2000 --max_queued_reads 40000 --chunks_per_caller 4000 --num_callers 32 \
  # --device "cuda:all:100%" -–chunk_size 500 --chunk_per_runner 1024 -–gpu_runners_per_device 8

  /g/data/lf10/as7425/apps/ont-guppy/bin/guppy_basecaller \
  --input_path ${fast5} \
  --save_path ${fastq} \
  --flowcell FLO-MIN106 --kit SQK-RNA002 --compress_fastq --recursive --num_callers 48 --disable_qscore_filtering

# --device "cuda:all:100%"
}; export -f basecall

# define the data for the function
export rep1="/g/data/xc17/degradation_project/Mg_degraded/data/MgCl_degrdaded_mild_pass1/fast5/"
export rep2="/g/data/xc17/degradation_project/Mg_degraded/data/MgCl_degrdaded_mild_pass2/fast5/"

# call the function
basecall "${rep1}" "${wd}/basecall_mild_degraded_rep1" || echo "failed for rep1" & # process forking
basecall "${rep2}" "${wd}/basecall_mild_degraded_rep2" || echo "failed for rep2" &

############################################################

# Align and filter the sequencing data

# environment management, using gadi native modules
module load samtools/1.10
module load minimap2

# set up working directory
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-08-22_basecall-new-samples"; mkdir -p ${wd}

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
