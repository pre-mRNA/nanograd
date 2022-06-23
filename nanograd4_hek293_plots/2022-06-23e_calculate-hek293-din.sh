#!/usr/bin/env Rscript

# written by AJ Sethi on 2022-06-23
# calculate DIN scores for all the HEK293 samples
# based on nanograd4 final output at /g/data/lf10/as7425/nanograd/analysis/2022-06-23_hek293-nanograd4-analysis/nanograd4_out

######################################
######################################

module load R

# link to script
export script="/home/150/as7425/nanograd/nanograd4/scripts/calculate_din.R"

# link to nanograd data output directory
export data="/g/data/lf10/as7425/nanograd/analysis/2022-06-23_hek293-nanograd4-analysis/nanograd4_out"

for i in ${data}/*txt; do printf "$(basename $i "_nanograd4_out.txt")\t$(Rscript ${script} ${i})\n" || echo "failed for $i" & done
wait && echo "$(date) ... done for all"
