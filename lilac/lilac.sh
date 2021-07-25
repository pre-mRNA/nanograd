#!/bin/bash

# written by A.J. Sethi on 2021-07-25
# aim: this script calculates base alignment rate and splice-in-reference rate for nanopore direct RNA alignments
# intron analysis only occurs for introns with more than 3 reads, can be configured at the regtools step

# dependencies: recent versions of: samtools, bedtools, regtools,  R, GNU parallel...

# testing;
# rm -rf "/scratch/lf10/as7425/lilac/lilac_test_data_sf3b1/output"; module load clairo; module load R; module load java; time bash /home/150/as7425/nanograd/lilac/lilac.sh -a "/g/data/lf10/as7425/genomes/mouse_genome" -b "/g/data/lf10/as7425/2020-11_mouseBrain/analysis/2021-07-12_all-pass-read-splicing/" -o "/scratch/lf10/as7425/lilac/lilac_test_data_sf3b1/output" -t "8"

# submit with

# -a /path/to/folder/with/annotations (i.e. fasta + gtf)
# -b /path/to/query/alignment.bam (musted be sorted and indexed with index located at /path/to/query/alignment.bam.bai)
# -o /path/to/output/output/directory
# -t /preferred/threadcount

####################################

# establish admin functions prior to proceeding

# die function
die() { printf "$(date +%F)\t$(date +%T)\t[Error] lilac died because it $*\n"; exit 1; }; export -f die

# new section notification
nsec() { printf "$(date +%F)\t$(date +%T)\t$*\n"; }; export -f nsec

# subsection notification
ssec() { printf "$(date +%F)\t$(date +%T)\t\t> $*\n"; }; export -f ssec

# verbose-subsection
vsec() { if [ "${VERBOSE}" == "TRUE" ]; then printf "$(date +%F)\t$(date +%T)\t\t> $*\n"; fi; }; export -f vsec


####################################

# check for packages and initialize some variables prior to parsing inputs
nsec "lilac initiated"

# test that all the requisite packages are avialable in path

# define a function that returns status 1 if a package isn't available
function checkAvailable {
  builtin type -P "$1" &> /dev/null
}

# iterate over the required packages
for targetMod in samtools bedtools parallel R java regtools; do
  checkAvailable $targetMod || die "cannot locate package $targetMod in user PATH"
done

# specify the default verbosity and threadcount
export threadCount="8"
export VERBOSE="false"

# initilalize runmode
export MODE="die"

# check the path to lilac
[ -n "${SCRIPTPATH+set}" ] || SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )" || die "cannot get script path"; export SCRIPTPATH="${SCRIPTPATH}" # get the script path, adapted from https://stackoverflow.com/questions/4774054/reliable-way-for-a-bash-script-to-get-the-full-path-to-itself
nsec "lilac is running from ${SCRIPTPATH} with parameters\n\t###\t$@\n\t###\tcommit\t$(git rev-parse --short HEAD 2>/dev/null | sed "s/\(.*\)/@\1/")"

####################################

# process the input argument

# flagged arguments
ARGS=""
while [ $# -gt 0 ]
do
  unset OPTIND
  unset OPTARG

  while getopts a:b:o:t: options
  do

    case $options in

      a) # annotation
      [ -d ${OPTARG} ] && export input="${OPTARG}" && inputFiles=`cd ${OPTARG} && ls *.*` || die "cannot access files in annotation directory" # get a list of files in the primary input
      gunzipCount=$(echo $inputFiles | grep -c -e ".gz"); if [ ${gunzipCount} -gt "0" ]; then ssec "unzipping ${gunzipCount} files in ${1}"; for i in ${1}/*; do gunzip ${i} || die "cannot unzip ${i} " & done; wait; fi # unzip any .gz files if they are present

      # ref fasta
      fastaCount=$(echo $inputFiles |  tr " " "\n" | grep -v ".fai" | grep -c -e ".fasta" -e ".fa" -e "fna"); [ ${fastaCount} -eq "1" ] || die "did not identify a single input fasta in ${1}, identified ${fastaCount}" # check that a fasta of some sort is present
      export myFasta="${input}/$(ls ${OPTARG} |  tr " " "\n" | grep -v ".fai" | grep -e ".fasta" -e ".fa" -e "fna")"

      # fasta index
      faiCount=$(echo $inputFiles |  tr " " "\n" | grep -c -e ".fai"); [ ${faiCount} -gt "1" ] && die "did not identify a single input fasta index in input, identified ${faiCount} -- try using an annotation directory with a single fasta index" # check that a fasta of some sort is present
      if [ ${faiCount} -gt "0" ]; then samtools faidx ${myFasta} > ${myFasta}.fai || die "cannot generate fasta index"; fi
      export myFai="${input}/$(ls ${OPTARG} |  tr " " "\n" | grep ".fai")"

      # genome file
      echo $inputFiles
      genomeCount=$(echo $inputFiles |  tr " " "\n" | grep -c -e ".genome"); [ ${genomeCount} -gt "1" ] && die "did not identify a single input genome index in input, identified ${genomeCount} -- try using an annotation directory with a single bedtools genome file" # check that a fasta of some sort is present
      if [ ${genomeCount} -eq "0" ]; then cat ${myFai} | awk -v OFS='\t' {'print $1,$2'} > ${myFai%??????}.genome  || die "cannot generate fasta index"; fi
      export myGenome="${input}/$(ls ${OPTARG} |  tr " " "\n" | grep ".genome")"

      # gtf
      GTFCount=$(echo $inputFiles | grep -c ".gtf"); [ "${GTFCount}" -eq "1" ] || die "did not identify a single input gtf in ${1}" # check that a single gtf is present
      export myAnnotation="${input}/$(ls ${OPTARG} | grep -e ".gtf")"
      ;;


      b) # path to alignment
      [ -d ${OPTARG} ] && export cbamDir="${OPTARG}" && cbamList=`ls -d ${OPTARG}/*.*` || die "ccannot find any binary alignments in mature bam directory" # get a list of files in the primary input
      export cbamCount=$(echo $cbamList | grep -c -e ".bam")
      export firstBam=$(echo $cbamList | tr " " "\n" | grep -e ".bam" | cut -f1 | head -n 1); # echo "your bam is ${firstBam}"
      export firstBamCount=$(samtools view -F 4 ${firstBam} | wc -l 2>/dev/null) #&& echo ${firstBamCount}
      [ ${cbamCount} -gt "0" ] || die "no binary alignments found in ${cbamDir}"
      ;;

      o) # output directory
      export outputDirectory="${OPTARG}"
      mkdir -p ${outputDirectory} 2>/dev/null
      [ -d ${outputDirectory} ] || die "cannot access the user-specified output directory, '${outputDirectory}'"
      ;;

      t) # threadCount
      if [ ! -z "${OPTARG}" ];
      then
        [[ "${OPTARG}" =~ ^-?[0-9]+$ ]] && export threadCount="${OPTARG}" || die "detected that user supplied invalid threadcount";
      else
        export threadCount="4" && ssec "setting threadCount to 4 by default"
      fi
      ;;

    esac;
  done
  shift $((OPTIND-1))
  ARGS="${ARGS} $1 "
  shift
done

####################################

# validate the input arguments

# validate the annotation
[ -f ${myAnnotation} ] && [ -f ${myFasta} ] || die "user supplied invalid annotation"
vsec "your annotation is ${myAnnotation}"
vsec "your reference sequence is ${myFasta}"

# validate output directory
[ -d "${outputDirectory}" ] || die "did not detect a valid output directory"
vsec "Your output directory is ${outputDirectory}"

# define some variables (doing this here because functions are being modularised and these variables might be required if an antecedent function is skipped)
export od=${outputDirectory}
export gi="${od}/generateIntrons"
export mp="${od}/mpileup"

# validate threadCount
ssec "proceeding with ${threadCount} threads"

# also define the quarter threadcount for downstream processing
export quarter_tc=$((${threadCount}/4)) # divid threadcout by 4
export qtc=$(echo $quarter_tc | awk '{print int($1+0.5)}') # round down

####################################

# we use the "main" function to call other functions as reqired
main() {
  analyseSplicing || die "cannot analyse splicing"
  echo "lilac exiting" && exit 0
}

# define a function to generate bed intervals of introns from a reference gtf
function analyseSplicing() {

  # make a directory for intronic bed intevals
  mkdir -p ${gi}

  # read the gtf and convert exons to bed intervals
  perl ${SCRIPTPATH}/source/gtf2bed.pl ${myAnnotation} > ${gi}/exons.bed || die "cannot extract exons"

  # read the bed file of exons and from this extract ranges which correspond to introns
  perl ${SCRIPTPATH}/source/get-intron-bed-from-exon-bed.pl ${gi}/exons.bed ${gi}/intron.bed || die "cannot extract introns"

  # get unique intronic interval range
  sort -k1,1 -k2,2n -k3,3n -u ${gi}/intron.bed | tail -n +2 > ${gi}/unique_intron_unrepaired.bed || die "cannot get unique introns"

  # repair intron file by shifting all starts by -1 units
  awk 'BEGIN {OFS="\t"}; {$2 = $2 - 1; print}' ${gi}/unique_intron_unrepaired.bed > ${gi}/unique_intron_repaired.bed || die "cannot repair intron file"

  regtools junctions extract -a 3 -m 50 -M 500000 -s 0 ${firstBam} | awk -F'\t' -v OFS='\t' '{split($11, blockSizes, ","); print $1,$2+blockSizes[1],$3-blockSizes[2],$4,$5,$6}' > ${gi}/observed-junctions.bed

  # count the number of introns
  totalIntronCount=$(awk '{ s+=$5 } END { print s }' ${gi}/observed-junctions.bed)
  #echo $totalIntronCount


  # reference introns
  refInt="${gi}/unique_intron_repaired.bed"

  # observed introns
  obsInt="${gi}/observed-junctions.bed"

  bedtools intersect -f 1 -F 1 -a $obsInt -b $refInt > ${gi}/observed-supported-junctions.bed || die "cannot identify supported junctions"

  # count the number of gtf-supported introns
  supportedIntronCount=$(awk '{ s+=$5 } END { print s }' ${gi}/observed-supported-junctions.bed)
  #echo $supportedIntronCount

  printf "total introns\t${totalIntronCount}\nsupported introns\t${supportedIntronCount}\n" > ${od}/splicing_metrics.txt

}; export -f analyseSplicing


main || die "cannot do main"

####################################
