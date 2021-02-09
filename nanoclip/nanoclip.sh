#!/bin/bash

export myFasta="/g/data/lf10/as7425/genomes/mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa"
export wtBam="/g/data/lf10/as7425/2020-11-05_feloniusGru/analysis/2021-02-02_nikolay-HL-1-crosslink/Nick_LK1-SQK-RNA002-NXL/Nick_LK1-SQK-RNA002-NXL.primary.sorted.bam"
export modBam="/g/data/lf10/as7425/2020-11-05_feloniusGru/analysis/2021-02-02_nikolay-HL-1-crosslink/Nick_LKL-SQK-RNA002-XL/Nick_LKL-SQK-RNA002-XL.primary.sorted.bam"
export window="5:134900000-149900000"
export outDir="/g/data/lf10/as7425/2020-11-05_feloniusGru/analysis/2021-02-08_build-nanoclip"
export myGenome="/g/data/lf10/as7425/genomes/mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly..genome"

mkdir -p ${output}

# do mpileup
# time samtools mpileup -d 0 -f ${myFasta} -B ${wtBam} ${modBam} -C 20 -Q 7 | awk '$7>0{print}'> ${scratch}/shrek.txt

# samtools extract window around chr5:134,900,000-149,900,000
for i in $wtBam $modBam
do echo $i; samtools view -b "${i}" "$window" > ${outDir}/window_${i##*/} || die "Cannot extract window for $i" &
done; wait && echo "almost done for all"
echo "done for all"

# filter for reads in a region
wtFilt="/g/data/lf10/as7425/2020-11-05_feloniusGru/analysis/2021-02-08_build-nanoclip/window_Nick_LK1-SQK-RNA002-NXL.primary.sorted.bam"
modFilt="/g/data/lf10/as7425/2020-11-05_feloniusGru/analysis/2021-02-08_build-nanoclip/window_Nick_LKL-SQK-RNA002-XL.primary.sorted.bam"

# desplit the data
for i in $wtBam $modBam
  do bedtools bamtobed -split -i ${outDir}/window_${i##*/} | bedtools bedtobam -ubam -g ${myGenome} -i - | samtools sort > ${outDir}/desplice_${i##*/}
  samtools index ${outDir}/desplice_${i##*/}
done

# export desplice names for pileup
wtDesplice="/g/data/lf10/as7425/2020-11-05_feloniusGru/analysis/2021-02-08_build-nanoclip/desplice_Nick_LK1-SQK-RNA002-NXL.primary.sorted.bam"
modDesplice="/g/data/lf10/as7425/2020-11-05_feloniusGru/analysis/2021-02-08_build-nanoclip/desplice_Nick_LKL-SQK-RNA002-XL.primary.sorted.bam"

# do mpileup
time samtools mpileup -d 0 -f ${myFasta} -B ${wtDesplice} ${modDesplice} -C 0 -Q 0 -q 0 | awk '$7>0{print}'> ${outDir}/2021-02-09_desplice_filtLib_reducedPileup.txt
gzip ${outDir}/2021-02-09_desplice_filtLib_reducedPileup.txt


##################################################
##################################################
##################################################
################################################## testing region
################################################## testing region
################################################## testing region

bedtools bamtobed -split -i /g/data/lf10/as7425/2020-11-05_feloniusGru/analysis/2021-02-08_build-nanoclip/window_Nick_LKL-SQK-RNA002-XL.primary.sorted.bam | bedtools bedtobam -ubam -g ${myGenome} -i - | samtools sort -l 0 > ${outDir}/desplice_${i##*/}





# file is 9gb
# take approx 5 min

time gzip ${scratch}/shrek.txt
# file is 855mb
# takes approx 2 min

exit
cd ~; scp as7425@gadi.nci.org.au:/scratch/lf10/as7425/shrek.txt.gz $(pwd)
