

export dataDir="/g/data/xc17/bk9031/2022_nanograd_bk/analysis/2021_HEK293-degradation-first4-AR/primaryAlignments"
export wd="/g/data/lf10/as7425/nanograd/analysis/2022-07-27_uLTRA-genomic-bigwigs-degradation"
module load samtools

for i in ${dataDir}/*bam; do
  echo $i
  nameLong=${i##*/}
  name=${nameLong%.*}
  echo $name
  bamCoverage -b ${i} -o "${wd}/${name}.bw" || echo "failed to make bw for ${name}" &
done

wait && echo "done for all"
