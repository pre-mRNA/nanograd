#!/bin/bash

# written by AJ Sethi on 2022-10-12
# executable script to run cheui

# test data
# export events="/g/data/xc17/degradation_project/nanopolish/nanopolish_5mM_MgCl_degradation_pass2.txt" # nanopolish output file
# export condition="degraded"
# export ctd="/g/data/lf10/as7425/nanograd/analysis/2022-05-11_run-cheui"

####################################################################################################
####################################################################################################

# housekeeping 'die' funvtion

# die function
die() { printf "$(date +%F)\t$(date +%T)\t[scriptDied] CHEUI_PIPELINE died because it $*\n"; exit 1; }; export -f die

# set up cheui run using positional arguments

# argument 1: eventAlign file
# argument 2: a folder in which the output is generated
# argument 3: a string, referring to the condition of the run
# argument 4: A or C, corresponding to detection of m6A or m5C (default A)

# define nanopolish data
export events="${1}" # nanopolish output file

# define library namea
export sample=$(basename ${events} .txt)

# define output directory for CHEUI
export ctd="${2}"
mkdir -p ${ctd} # working directory
cd ${ctd} || die "cannot access user-supplied working directory ${ctd}"

# define library condition
export condition=${3}

# test the condition
[ -z "$condition" ] && die "did not detect a valid condition in user-supplied ARG3"

# get the run mode
export mode=${4}
# test the mode
[ -z "$mode" ] && die "did not detect a valid condition in user-supplied ARG3"

# get script paths
[ -n "${SCRIPTPATH+set}" ] || SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )" || die "cannot get script path"
export SCRIPTPATH="${SCRIPTPATH}" # get the script path # taken from https://stackoverflow.com/questions/4774054/reliable-way-for-a-bash-script-to-get-the-full-path-to-itself
echo "scriptpath is ${SCRIPTPATH}"

# get cheui paths
export cheui="${SCRIPTPATH}/CHEUI" # CHEUI install directory
export k_model="${cheui}/kmer_models/model_kmer.csv" # kmer model

# get deep learning models for m6A
export dm_m6a_1="${cheui}/CHEUI_trained_models/CHEUI_m6A_model1.h5" # m6a deep learning model 1
export dm_m6a_2="${cheui}/CHEUI_trained_models/CHEUI_m6A_model2.h5" # m6a deep learning model 2

# get deep learning models for m5C
export dm_m5c_1="${cheui}/CHEUI_trained_models/CHEUI_m5C_model1.h5" # m6a deep learning model 1
export dm_m5c_2="${cheui}/CHEUI_trained_models/CHEUI_m5C_model2.h5" # m6a deep learning model 2

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
if [ ${mode} == "A" ]
then

  printf "$(date) .... Using m6A models\n\n\n"

  # preprocess the data
  printf "$(date) .... starting to preprocess ${sample}\n\n\n"

  python3 "${cheui}/scripts/CHEUI_preprocess_m6A.py" -i "${events}" -m ${k_model} -s "${sample}_m6A" -n 48 -o "${ctd}" || die "$(date) ..... preprocssing failed\n"

  printf "$(date) ..... done preprocessing ${sample}\n\n\n"

  # run model I
  printf "$(date) .... starting model 1 ${sample}\n\n\n"

  python3 "${cheui}/scripts/CHEUI_predict_model1.py" -i "${ctd}/${sample}_signals+IDS.p" -m "${dm_m6a_1}" -l "${condition}" -o "${ctd}/${sample}_model1_m6A.txt" || die "$(date) ..... model 1 failed\n"

  printf "$(date) ..... done with model 1${sample}\n\n\n"

  # sort the output
  (head -1 "${ctd}/${sample}_model1.txt" && tail -n +2 "${ctd}/${sample}_model1.txt" | sort -k1 --parallel=48) > "${ctd}/${sample}_model1_m6A_sorted.txt" || die "sort failed"

  # launch model 2
  printf "$(date) .... starting model 2 ${sample}\n\n\n"

  python3 "${cheui}/scripts/CHEUI_predict_model2.py" -i "${ctd}/${sample}_model1_m6A_sorted.txt" \
    -m "${dm_m6a_2}" \
    -o "${ctd}/${sample}_model2_m6A.txt" || die "model 2 failed\n"

  printf "$(date) ..... finished model 2 ${sample}\n\n\n"

elif [ ${mode} == "C" ]
then

  printf "$(date) .... Using m5C models\n\n\n"

  # preprocess the data
  printf "$(date) .... starting to preprocess ${sample}\n\n\n"

  python3 "${cheui}/scripts/CHEUI_preprocess_m5C.py" -i "${events}" -m ${k_model} -s "${sample}_m5C" -n 48 -o "${ctd}" || die "$(date) ..... preprocssing failed\n"

  printf "$(date) ..... done preprocessing ${sample}\n\n\n"

  # run model I
  printf "$(date) .... starting model 1 ${sample}\n\n\n"

  python3 "${cheui}/scripts/CHEUI_predict_model1.py" -i "${ctd}/${sample}_signals+IDS.p" -m "${dm_m5C_1}" -l "${condition}" -o "${ctd}/${sample}_model1_m5C.txt" || die "$(date) ..... model 1 failed\n"

  printf "$(date) ..... done with model 1${sample}\n\n\n"

  # sort the output
  (head -1 "${ctd}/${sample}_model1.txt" && tail -n +2 "${ctd}/${sample}_model1.txt" | sort -k1 --parallel=48) > "${ctd}/${sample}_model1_m5C_sorted.txt" || die "sort failed"

  # launch model 2
  printf "$(date) .... starting model 2 ${sample}\n\n\n"

  python3 "${cheui}/scripts/CHEUI_predict_model2.py" -i "${ctd}/${sample}_model1_m5C_sorted.txt" \
    -m "${dm_m5c_2}" \
    -o "${ctd}/${sample}_model2_m5C.txt" || die "model 2 failed\n"

  printf "$(date) ..... finished model 2 ${sample}\n\n\n"

else

  die "user supplied invalid mode ${mode}"
fi
