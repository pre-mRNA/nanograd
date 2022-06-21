#!/bin/bash

# Aim: Extract minknow read end classification calls
# Written by AJ Sethi on 2022-06-06

##################################################
##################################################

export ss1="/g/data/lf10/as7425/2020-11-05_feloniusGru/data/2021-10-07_srsf3-DIRCLIP-rep1/20211005_0420_MN33901_FAQ68968_22bee822/sequencing_summary_FAQ68968_c2ab58a6.txt"
# path to sequencing_summary 1

export ss2="/g/data/lf10/as7425/2020-11-05_feloniusGru/data/2022-02-11_kw-hmga1-jcsmr-round1/HMGA1A_non_UV_crosslinked_S1/20220130_1406_MN18034_FAR38530_889f2d78/sequencing_summary_FAR38530_232bfccf.txt"
# path to sequencing_summary 1

# path to working directory
export lwd="/scratch/lf10/as7425/tmp5"

### extract data
cat $ss1 | cut -f22 | tail -n +2 | sort | uniq -c > ${lwd}/ss1.txt
cat $ss2 | cut -f22 | tail -n +2 | sort | uniq -c > ${lwd}/ss2.txt
