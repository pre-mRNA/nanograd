#!/bin/bash
# written by BK on 2022-10-14
# aim: script for R2Dtools

#############################################################
#############################################################

# download R2Dtool from github
git clone git@github.com:comprna/R2Dtool.git
cd R2Dtool

export r2dtool=""
export m6A_bed=""
export m6C_bed=""
export annotate="/g/data/xc17/bk9031/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"
export annotated_output_m6A=""
export annotated_output_m5C=""
export liftover_output_m6A=""
export liftover_output_m5C=""

# converting CHEUI m6A model2 to bed like input file
bash "${r2dtool}/scripts/cheui_to_bed.sh" ["${output}/m6A"] ["m6A_bed"]

bash "${r2dtool}/scripts/cheui_to_bed.sh" ["${output}/m5C"] ["m5C_bed"]

# annotating transcriptomic sites
Rscript "${r2dtool}/scripts/R2_annotate.R" ["${m6A_bed}"] ["${annotation}"] ["${annotated_output_m6A}"]

Rscript "${r2dtool}/scripts/R2_annotate.R" ["${m5C_bed}"] ["${annotation}"] ["${annotated_output_m5C}"]

# liftover transcriptomic sites
Rscript "${r2dtool}/scripts/R2_lift.R" ["${m6A_bed}"] ["${annotation}"] ["${liftover_output_m6A}"]

Rscript "${r2dtool}/scripts/R2_lift.R" ["${m5C_bed}"] ["${annotation}"] ["${liftover_output_m5C}"] 
