#!/bin/bash

# link script
export pipeline="/home/150/as7425/nanograd/cheui/2022-10-12_cheui-runnable-script.sh"
export data="/g/data/xc17/degradation_project/Mg_degraded/data"
export output="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/2022-10-10_run-CHEUI"
export annotation="/g/data/xc17/bk9031/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"

# test for sample wt_1
bash ${pipeline} "${annotation}" "${output}/m6A" "wt" "A"

bash ${pipeline} "${annotation}" "${output}/m5C" "wt" "C"

# test for sample wt_2 
bash ${pipeline} "${annotation}" "${output}/m6A" "wt" "A"

bash ${pipeline} "${annotation}" "${output}/m5C" "wt" "C"

# test for mild_deg_1
bash ${pipeline} "${annotation}" "${output}/m6A" "mild_deg" "A"

bash ${pipeline} "${annotation}" "${output}/m5C" "mild_deg" "C"

# test for mild_deg_2
bash ${pipeline} "${annotation}" "${output}/m6A" "mild_deg" "A"

bash ${pipeline} "${annotation}" "${output}/m5C" "mild_deg" "C"

# test for deg_1
bash ${pipeline} "${annotation}" "${output}/m6A" "deg" "A"

bash ${pipeline} "${annotation}" "${output}/m5C" "deg" "C"

# test for deg_2
bash ${pipeline} "${annotation}" "${output}/m6A" "deg" "A"

bash ${pipeline} "${annotation}" "${output}/m5C" "deg" "C"
