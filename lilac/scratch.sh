# make mouse testing data on sf3b1 locus
{
# original bam
/g/data/lf10/as7425/2020-11_mouseBrain/analysis/2021-07-12_all-pass-read-splicing/sorted.all_pass_reads.bam

# print header
samtools view -H /g/data/lf10/as7425/2020-11_mouseBrain/analysis/2021-07-12_all-pass-read-splicing/sorted.all_pass_reads.bam > /scratch/lf10/as7425/lilac/sf3b1.sam


# print sam entries
samtools view -b /g/data/lf10/as7425/2020-11_mouseBrain/analysis/2021-07-12_all-pass-read-splicing/sorted.all_pass_reads.bam "1:54,928,488-55,090,062" | samtools sort | samtools view >> /scratch/lf10/as7425/lilac/sf3b1.sam

# convert to bam
samtools view -b /scratch/lf10/as7425/lilac/sf3b1.sam >  /scratch/lf10/as7425/lilac/sf3b1.bam
samtools index /scratch/lf10/as7425/lilac/sf3b1.bam
}

# test intron extraction on sf3b1 locus
{

  samtools view -h /scratch/lf10/as7425/lilac/sf3b1.bam | awk '$0 ~ /^@/ || $6 ~ /N/' | samtools view -bS - | bamToBed -bed12 > /scratch/lf10/as7425/lilac/sf3b1.bed12 || die "cannot make spliced reads"

  bed12ToBed6 -i /scratch/lf10/as7425/lilac/sf3b1.bed12 > /scratch/lf10/as7425/lilac/sf3b1.bed6 || die "cannot make bed6 spliced reads"
  subtractBed -a /scratch/lf10/as7425/lilac/sf3b1.bed12 -b /scratch/lf10/as7425/lilac/sf3b1.bed6 -s | cut -f 1-6 > /scratch/lf10/as7425/lilac/sf3b1.introns.bed || die "cannot make intronic reads"




}

# test regtools
{
  /g/data/lf10/as7425/apps/regtools/build/regtools
  # extract junctions
  regtools junctions extract -a 3 -m 40 -s 0 /scratch/lf10/as7425/lilac/sf3b1.bam > /scratch/lf10/as7425/lilac/sf3b1_introns.bed
  # regtools returns exons, let's see if we can extract introns


  regtools junctions extract -a 3 -m 50 -M 500000 -s 0 /scratch/lf10/as7425/lilac/sf3b1.bam | awk -F'\t' -v OFS='\t' '{split($11, blockSizes, ","); print $1,$2+blockSizes[1],$3-blockSizes[2],$4,$5,$6}' > /scratch/lf10/as7425/lilac/sf3b1_introns.bed


  # count the number of introns
  head /scratch/lf10/as7425/lilac/sf3b1_introns.bed
  totalIntronCount=$(awk '{ s+=$5 } END { print s }' /scratch/lf10/as7425/lilac/sf3b1_introns.bed)
  echo $totalIntronCount
}

# compare junctions to reference introns
{

# reference introns
refInt="/scratch/lf10/as7425/lilac/mouse-introns-repaired.bed"

# observed introns
obsInt="/scratch/lf10/as7425/lilac/sf3b1_introns.bed"

bedtools intersect -f 1 -F 1 -a $obsInt -b $refInt > /scratch/lf10/as7425/lilac/match_perfect.bed

# count the number of gtf-supported introns
supportedIntronCount=$(awk '{ s+=$5 } END { print s }' /scratch/lf10/as7425/lilac/match_perfect.bed)
echo $supportedIntronCount

}
