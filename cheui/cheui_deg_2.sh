#!/bin/bash

# written by AJ Sethi on 2022-05-05
# Run CHEUI for sample AR1

####################################################################################################
####################################################################################################

# set up cheui run

# define nanopolish data
export events="/g/data/xc17/degradation_project/nanopolish/nanopolish_5mM_MgCl_degradation_pass2.txt" # nanopolish output file

# define library condition
export condition="degraded"

# define output directory for CHEUI
export ctd="/g/data/lf10/as7425/nanograd/analysis/2022-05-11_run-cheui"; mkdir -p ${ctd} # working directory

# define cheui paths
export cheui="/g/data/lf10/as7425/apps/CHEUI" # CHEUI install directory
export k_model="${cheui}/kmer_models/model_kmer.csv" # kmer model
export dm_m6a_1="${cheui}/CHEUI_trained_models/CHEUI_m6A_model1.h5" # m6a deep learning model 1
export dm_m6a_2="${cheui}/CHEUI_trained_models/CHEUI_m6A_model2.h5" # m6a deep learning model 2

# define library namea
export sample=$(basename ${events} .txt)

# establish environment

# for conda install
# source ~/.bashrc
# conda activate cheui

# for gadi install
export PYTHONPATH="${PYTHONPATH}:/g/data/xc17/pm1122/Epinano_IVT/scripts/"
module load tensorflow/2.3.0
module load cuda/10.1

##################################################

# run cheui

# run cheui preprocessing
printf "$(date) .... starting to preprocess ${sample}\n\n\n"

time python3 "${cheui}/scripts/CHEUI_preprocess_m6A.py" -i "${events}" -m ${k_model} -s "${sample}" -n 48 -o "${ctd}" || printf "$(date) ..... preprocssing failed\n"

printf "$(date) ..... done preprocessing ${sample}\n\n\n"

# run model I
printf "$(date) .... starting model 1 ${sample}\n\n\n"

time python3 "${cheui}/scripts/CHEUI_predict_model1.py" -i "${ctd}/${sample}_signals+IDS.p" -m "${dm_m6a_1}" -l "${condition}" -o "${ctd}/${sample}_model1.txt" || printf "$(date) ..... model 1 failed\n"

printf "$(date) ..... done with model 1${sample}\n\n\n"

# sort the output
sort -k1 --parallel=48 "${ctd}/${sample}_model1.txt" > "${ctd}/${sample}_model1_sorted.txt"

# launch model 2
printf "$(date) .... starting model 2 ${sample}\n\n\n"

time python3 "${cheui}/scripts/CHEUI_predict_model2.py" -i "${ctd}/${sample}_model1_sorted.txt" \
  -m "${dm_m6a_2}" \
  -o "${ctd}/${sample}_model2.txt" || printf "$(date) ..... model 2 failed\n"

printf "$(date) ..... finished model 2 ${sample}\n\n\n"

# model 2 completes succesfully

##################################################
