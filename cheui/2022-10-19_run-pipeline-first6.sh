#!/bin/bash

# written by
# aim

# define global variables for all runs
export pipeline="/home/150/as7425/nanograd/cheui/2022-10-12_cheui-runnable-script.sh"
export output="/g/data/lf10/as7425/nanograd/analysis/2022-10-19_run-cheui-first6/"; mkdir -p ${output} 2>/dev/null
export annotation="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"

# get the data directory
export data_dir="/g/data/xc17/as7425/sharing/2022-10-10_nanopolish-first6"

# list the samples
export control_1="${data_dir}/all.undegraded_hek293_pass1.fastq.gz.sorted.bam_sorted_primary_events.tsv"
export control_2="${data_dir}/all.undegraded_hek293_pass2.fastq.gz.sorted.bam_sorted_primary_events.tsv"
export mild_1="${data_dir}/primary_mild_degradataion_rep1.fastq.bam_sorted_primary_events.tsv"
export mild_2="${data_dir}/primary_mild_degradataion_rep2.fastq.bam_sorted_primary_events.tsv"
export harsh_1="${data_dir}/all.5mM_MgCl_degrdation_pass1.fastq.gz.sorted.bam_sorted_primary_events.tsv"
export harsh_2="${data_dir}/all.5mM_MgCl_degrdation_pass2.fastq.gz.sorted.bam_sorted_primary_events.tsv"

# log directory
export logs="/home/150/as7425/logs"

# 12 jobs; 6 samples * 2 models
export submit="qsub -P lf10 -l walltime=24:00:00,mem=190GB,ncpus=48,ngpus=4,storage=scratch/lf10+gdata/lf10+gdata/xc17,jobfs=100GB -q gpuvolta -o ${logs} -e ${logs} ${pipeline}"

# setup positional arguments
${submit} "${control_1}" "${annotation}" "${output}" "control_1" "A"

# test for m6A
bash ${pipeline} "${data}" "${output}/m6A" "test_data" "A"

# test for m5C
bash ${pipeline} "${data}" "${output}/m5C" "test_data" "C"


# works
qsub -P lf10 -l walltime=24:00:00,mem=190GB,ncpus=48,ngpus=4,storage=scratch/lf10+gdata/lf10+gdata/xc17,jobfs=100GB -q gpuvolta -o ${logs} -e ${logs} -v a="${data}",b="${output}/m6A",c="test_data",d="A" /home/150/as7425/nanograd/cheui/2022-10-12_cheui-runnable-script.sh

# fails
qsub -P lf10 -l walltime=24:00:00,mem=190GB,ncpus=48,ngpus=4,storage=scratch/lf10+gdata/lf10+gdata/xc17,jobfs=100GB -q gpuvolta -o ${logs} -e ${logs} -v a="${data}",b="${output}/m6A",c="test_data",d="A" /home/150/as7425/nanograd/cheui/2022-10-12_cheui-runnable-script.sh

# works
qsub -P lf10 -l walltime=24:00:00,mem=190GB,ncpus=48,ngpus=4,storage=scratch/lf10+gdata/lf10+gdata/xc17,jobfs=100GB -q gpuvolta -o ${logs} -e ${logs} -v 1="${data}",2="${output}/m6A",3="test_data",4="A" /home/150/as7425/nanograd/cheui/2022-10-12_cheui-runnable-script.sh

qsub -P lf10 -l walltime=24:00:00,mem=190GB,ncpus=48,ngpus=4,storage=scratch/lf10+gdata/lf10+gdata/xc17,jobfs=100GB -q gpuvolta -o ${logs} -e ${logs} -v 1="${data}",2="${output}/m6A",3="test_data",4="A" /home/150/as7425/nanograd/cheui/2022-10-19_echo-positional-arguments.sh

# doesn't work
export logs="/home/150/as7425/logs"2
qsub -P lf10 -l walltime=1:00:00,mem=1GB,ncpus=1,storage=scratch/lf10+gdata/lf10+gdata/xc17,jobfs=1GB -q express -o ${logs} -e ${logs} -v 1="${data}",2="${output}/m6A",3="test_data",4="A" /home/150/as7425/nanograd/cheui/2022-10-19_echo-positional-arguments.sh


qsub -P lf10 -l walltime=1:00:00,mem=1GB,ncpus=1,storage=scratch/lf10+gdata/lf10+gdata/xc17,jobfs=1GB -q express -o ${logs} -e ${logs} -v "${data}","${output}/m6A","test_data","A" /home/150/as7425/nanograd/cheui/2022-10-19_echo-positional-arguments.sh


qsub -P lf10 -l walltime=1:00:00,mem=1GB,ncpus=1,storage=scratch/lf10+gdata/lf10+gdata/xc17,jobfs=1GB -q express -o ${logs} -e ${logs} -- bash /home/150/as7425/nanograd/cheui/2022-10-19_echo-positional-arguments.sh "${data}" "${output}/m6A" "test_data" "A"
