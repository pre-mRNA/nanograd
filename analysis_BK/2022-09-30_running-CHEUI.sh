#!/bin/bash
# written by BK on 2022-09-30
# aim: script for running CHEUI

#############################################################
#############################################################

# define nanopolish data which is the output file of the nanopolish
export events="/g/data/xc17/as7425/sharing/2022-10-10_nanopolish-first6"

# define library conditions
export condition = "degraded/undegraded"

# defining output directory for CHEUI
export out_dir="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/2022-10-10_run-CHEUI"

# defining CHEUI paths
export run_CHEUI="/g/data/xc17/bk9031/apps/CHEUI"
export k_mer_models="${run_CHEUI}/kmer_models/model_kmer.csv"
export trained_model1_m6A="${run_CHEUI}/CHEUI_trained_models/CHEUI_m6A_model1.h5"
export trained_model2_m6A="${run_CHEUI}/CHEUI_trained_models/CHEUI_m6A_model2.h5"

# define library names
export sample=$(basename${events}.txt)

# for gadi install
export PYTHONPATH="${PYTHONPATH}:/g/data/xc17/pm1122/Epinano_IVT/scripts/"
module load tensorflow/2.3.0
module load cuda/10.1

##########################################################################################

# run CHEUI

echo "$(date)...starting CHEUI preprocessing ${sample}\n\n\n"
time python3 "${run_CHEUI}/scripts/CHEUI_preprocess_m6A.py" -i "${events}" -m ${k_mer_models} -s "${sample}" -n 48 -o "${out_dir}" || echo "$(date)...preprocessing failed\n"
echo "$(date)...done CHEUI preprocessing ${sample}\n\n\n"

# run model 1
echo "$(date)...CHEUI Model 1 started ${sample}\n\n\n"
time python3 "${run_CHEUI}/scripts/CHEUI_predict_model1.py" -i "${out_dir}/${sample}_signal+IDS.p" -m "${trained_model1_m6A}" -l "${condition}" -o "${out_dir}/${sample}_model1.txt" || echo "$(date)...failed model 1\n"
echo "$(date)...done CHEUI model 1${sample}\n\n\n"

# sorting the predictions to group all predictions of the same site
sort -k1 --parallel=48 "${out_dir}/${sample}_model1.txt" > "${out_dir}/${sample}_model1_sorted.txt" || echo "$(date)...sorting failed\n"

# running model 2
echo "$(date)...CHEUI model 2 started ${sample}\n\n\n"
time python3 "${run_CHEUI}/scripts/CHEUI_predict_model2.py" -i "${out_dir}/${sample}_model1_sorted.txt" -m "${trained_model2_m6A}" -o "${out_dir}/${sample}_model2.txt" || echo "$(date)...failed model 2\n"
echo "$(date)...done CHEUI model 2 ${sample}\n\n\n"
