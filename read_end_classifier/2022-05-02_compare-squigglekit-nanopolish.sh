#!/bin/bash

# written by AJ Sethi on 2022-05-02
# aim: to compare squigglekit to nanopolish

# test data: random reads from Kat's sequencing for Bootes nep1 (that has already been nanopolished)

######################################################################
######################################################################

# make a directory to fetch read ID
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-05-02_test-squigglekit"
mkdir -p ${wd}/nanopolish_reads

# fetch 5 reads from Kat's mutant_rep1 (nep1) sequencing containing the modified site
cd /g/data/xc17/kat_nm/NanoPolish/NEP1_Replicate1

# get a list of 5 target reads
time head -n 50000000 *txt | awk '($2 == 2137){print}' | cut -f4 | uniq | head -n 5 > ${wd}/read_targets.txt

# get the name of an example read
head -n 1 ${wd}/read_targets.txt
# a3f8f3ad-594f-4260-a757-631cc0b40ab9

# for each target, extract the data from annopolish eventAlign output
for i in $(cat ${wd}/read_targets.txt); do cat <(cat *txt | head -n 1) <(head -n 5000000 *txt | grep $i)  > ${wd}/nanopolish_reads/${i}.txt; done

# useing squigglekit now;
conda env list
conda activate squigglekit

# fast5 path for mut_rep2
export mr2seq="/g/data/xc17/kat_nm/Sequencing_round2/NEP1_S2/NEP1_KD_sample2/20211029_1236_MN38051_FAP77920_e9bffbe8/fast5_gzip"

# make a directory
mkdir ${wd}/squigglekit

######################

# fetch a fast5

# print target name
echo "a3f8f3ad-594f-4260-a757-631cc0b40ab9" | gzip -c > ${wd}/squigglekit/target1.txt.gz

# zip sequencing summary file
cat /g/data/xc17/kat_nm/fastqs/NEP1_Mut_Replicate1/sequencing_summary.txt | gzip -c > ${wd}/squigglekit/sequencing_summary.txt.gz

# call fast5 fetcher
export export wd="/g/data/lf10/as7425/nanograd/analysis/2022-05-02_test-squigglekit"
export mr2seq="/g/data/xc17/kat_nm/Sequencing_round2/NEP1_S2/NEP1_KD_sample2/20211029_1236_MN38051_FAP77920_e9bffbe8/fast5_gzip"
fast5_fetcher="/g/data/lf10/as7425/apps/SquiggleKit/fast5_fetcher_multi.py"

python ${fast5_fetcher} -f "${wd}/squigglekit/target1.txt.gz" -s "${wd}/squigglekit/sequencing_summary.txt.gz" -m "${mr2seq}" -o "${wd}/squigglekit/fast5"
python ${fast5_fetcher} -s "${wd}/squigglekit/sequencing_summary.txt.gz" -m "/g/data/xc17/kat_nm/Sequencing_round2/NEP1_S2/NEP1_KD_sample2/20211029_1236_MN38051_FAP77920_e9bffbe8/fast5_gzip" -o "${wd}/squigglekit/fast5"

# squigglepull
SquigglePull="/g/data/lf10/as7425/apps/SquiggleKit/SquigglePull.py"

##########

python $SquigglePull -p ${mr2seq} -v -r > ${wd}/squigglekit/read_1.tsv

# download data and open in R

# search the original file for a read
# fetch 5 reads from Kat's mutant_rep1 (nep1) sequencing containing the modified site
eventAlign="/g/data/xc17/kat_nm/NanoPolish/NEP1_Replicate1/NEP_Replicate1_nanopolish.txt"
read="071928ed-64d3-4b5e-8108-a7b0a1ca07d4"
module load parallel

zcat ${eventAlign} | grep "${read}"

# final read of interest: 071928ed-64d3-4b5e-8108-a7b0a1ca07d4

eventAlign="/g/data/xc17/kat_nm/NanoPolish/NEP1_Replicate1/NEP_Replicate1_nanopolish.txt"
read="071928ed-64d3-4b5e-8108-a7b0a1ca07d4"
module load parallel
export export wd="/g/data/lf10/as7425/nanograd/analysis/2022-05-02_test-squigglekit"

time < "${eventAlign}" parallel --pipe --block 1000M grep -i -C 5 "${read}" >> ${wd}/out_071928ed-64d3-4b5e-8108-a7b0a1ca07d4_eventAlign.tsv
