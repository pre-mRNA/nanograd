#!/bin/bash

export in="/Users/asethi/localGadiData/2022-05-14_cheui-degradation/"
export out="/Users/asethi/localGadiData/2022-05-14_cheui-degradation/"; mkdir -p ${out}
export anno="/Users/asethi/localGadiData/2022-05-11_nanograd-metaplots-R/Homo_sapiens.GRCh38.105.chr.gtf"

cd ${out}
for i in ${in}/*_model2.txt; do

  sample=$(basename $i _model2.txt)
  echo $sample
  # for rat
  bash /Users/asethi/Documents/asethi-jcsmr/femtoSplice/scripts/2022-01-31_cheui-annotate2.sh "II" ${anno} "${i}" "${out}/${sample}_liftover.txt"

done
