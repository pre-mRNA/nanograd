#!/bin/bash

export myFasta="/g/data/lf10/as7425/genomes/mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa"
export wtBam="/g/data/lf10/as7425/2020-11-05_feloniusGru/analysis/2021-02-02_nikolay-HL-1-crosslink/Nick_LK1-SQK-RNA002-NXL/Nick_LK1-SQK-RNA002-NXL.primary.sorted.bam"
export modBam="/g/data/lf10/as7425/2020-11-05_feloniusGru/analysis/2021-02-02_nikolay-HL-1-crosslink/Nick_LKL-SQK-RNA002-XL/Nick_LKL-SQK-RNA002-XL.primary.sorted.bam"
export window="5:134900000-149900000"
export outDir="/g/data/lf10/as7425/2020-11-05_feloniusGru/analysis/2021-02-08_build-nanoclip"

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

# do mpileup
time samtools mpileup -d 0 -f ${myFasta} -B ${wtFilt} ${modFilt} -C 20 -Q 7 | awk '$7>0{print}'> ${outDir}/2021-02-08_filtLib_reducedPileup.txt
gzip ${outDir}/2021-02-08_filtLib_reducedPileup.txt



# file is 9gb
# take approx 5 min

time gzip ${scratch}/shrek.txt
# file is 855mb
# takes approx 2 min

exit
cd ~; scp as7425@gadi.nci.org.au:/scratch/lf10/as7425/shrek.txt.gz $(pwd)
