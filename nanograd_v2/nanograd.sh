#!/bin/bash

# written by A.J. Sethi on 2020-11-18
# aim: this script computes the coverage gradient for each individual transcript, based on the alignemnt of nanopore direct RNA alignments to a reference CDS

# dependencies: recent versions of: samtools, bedtools, R, GNU parallel
# R modules; tidyverse, runner, zoo, data.table, scales, ggpubr

# testing;
# rm -rf "/scratch/lf10/as7425/bamTest"; module load clairo; module load R; module load java; time bash ~/nanograd/nanograd.sh cluster -a "/scratch/lf10/as7425/sequinData" -b "/scratch/lf10/as7425/sequinData/" -o "/scratch/lf10/as7425/bamTest" -t "8"

# submit with

# -a /path/to/folder/with/annotations (i.e. fasta + gtf)
# -b /path/to/query/alignment.bam (musted be sorted and indexed with index located at /path/to/query/alignment.bam.bai)
# -o /path/to/output/output/directory
# -t /preferred/threadcount
# -m /preferred/runmode (partially implemented)
# -v /preferred/verbosity (partially implemented)

####################################

# establish admin functions prior to proceeding

# die function
die() { printf "$(date +%F)\t$(date +%T)\t[Error] nanograd died because it $*\n"; exit 1; }; export -f die

# new section notification
nsec() { printf "$(date +%F)\t$(date +%T)\t$*\n"; }; export -f nsec

# subsection notification
ssec() { printf "$(date +%F)\t$(date +%T)\t\t> $*\n"; }; export -f ssec

# verbose-subsection
vsec() { if [ "${VERBOSE}" == "TRUE" ]; then printf "$(date +%F)\t$(date +%T)\t\t> $*\n"; fi; }; export -f vsec


####################################

# check for packages and initialize some variables prior to parsing inputs
nsec "nanograd initiated"

# test that all the requisite packages are avialable in path

# define a function that returns status 1 if a package isn't available
function checkAvailable {
  builtin type -P "$1" &> /dev/null
}

# iterate over the required packages
for targetMod in samtools bedtools parallel R java; do
  checkAvailable $targetMod || die "cannot locate package $targetMod in user PATH"
done

# specify the default verbosity and threadcount
export threadCount="8"
export VERBOSE="false"

# initilalize runmode
export MODE="die"

# check the path to nanograd
[ -n "${SCRIPTPATH+set}" ] || SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )" || die "cannot get script path"; export SCRIPTPATH="${SCRIPTPATH}" # get the script path, adapted from https://stackoverflow.com/questions/4774054/reliable-way-for-a-bash-script-to-get-the-full-path-to-itself
nsec "Nanograd is running from ${SCRIPTPATH} with parameters\n\t###\t$@\n\t###\tcommit\t$(git rev-parse --short HEAD 2>/dev/null | sed "s/\(.*\)/@\1/")"

####################################

# process the input arguments
# process runmode and flagged arguments separately

# runmode
export urm=$1 # user run mode
case $urm in
  (decay | Decay) export runmode="d" ;;
  (cluster | Cluster) export runmode="c" ;;
  (*) nsec "user provided invalid runmode -- try 'decay' or 'cluster' next time" && die "detected invalid user runmode" ;;
esac

# flagged arguments
ARGS=""
while [ $# -gt 0 ]
do
  unset OPTIND
  unset OPTARG

  while getopts a:b:o:t:m:v:p: options
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

      v) # verbosity
      if [ ${OPTARG} == "TRUE" ]; then export VERBOSE="TRUE"; else export VERBOSE="false"; fi
      ssec "Proceeding in verbose mode"
      ;;



      # runmode currently deprecated
      #m) # runmode; test if we want to run the entire index or just a specific region

      # note:
      # each mode has a numerical start (nm) and end integer (nmf)
      # the start integer tells nanograd what process to start with
      # nanograd will proceed until the end integer is reached
      # or if no end integer is defined
      # then until the script ends

    #   # complete mode means we do every module (i.e. start at 1 and end at n)
    #   if [[ ${OPTARG} == "complete" ]]; then
    #     export MODE="complete"
    #     export nm="1"
    #     export nmf="999"

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

# report confidence level
export cl="10"
vsec "your cluster confidence level is 25 alignments"

# define some variables (doing this here because functions are being modularised and these variables might be required if an antecedent function is skipped)
export od=${outputDirectory}
export mp="${od}/mpileup"
export sb="${od}/splitBamThreadCount/"
export ec="${od}/multiSamToEndCoordinate/" # end coordinate
export ac="${od}/analyseClusters/"
export cr="${od}/assignClusterToRead"
export ao="${od}/assembleOutput"

# validate threadCount
ssec "proceeding with ${threadCount} threads"

# also define the quarter threadcount for downstream processing
export quarter_tc=$((${threadCount}/4)) # divid threadcout by 4
export qtc=$(echo $quarter_tc | awk '{print int($1+0.5)}') # round down

# validate the runmode
#[ ${MODE} == "die" ] && die "user did not define runmode"
# [[ ${MODE} == "die" ]] && nsec "exiting due to mode die" && exit 0
# ssec "your runmode is ${MODE}"

echo "testing done" && exit 0
####################################

# we use the "main" function to call oter functions as reqired
main() {

  splitBamThreadCount || die "cannot split"
  multiSamToEndCoordinate || die "cannot fetch end coordinates"
  clusterStatistics || die "cannot analyse 3' end clusters"
  assignClusterToRead || die "cannot assign clusters to reads"
  assembleOutput || die "cannot assembly output"
  nsec "nanograd completed succesfully - thanks for testing it :)"
  exit 0


  # if [ "$nm" -lt "2" ]
  # then
  #   nsec "Calculating per base read coverage using samtools mpileup"
  #   perBaseCoverage || die "unable to produce whippet index"
  # fi
  #
  # # generate whippet index
  # if [ "$nm" -lt "3" ]
  # then
  #   echo "Generating whippet index"
  #   makeWhippetIndex || die "unable to prepare whippet index"
  # fi
  #
  # # generate fastq
  # if [ "$nm" -lt "4" ]
  # then
  #   echo "Preparing to analyse cytoplasmic alignments"
  #   generateSequence || die "unable to prepare whippet index"
  # fi
  #
  # # call whippet for mature transcripts
  # if [ "$nm" -lt "5" ]
  # then
  #   echo "Preparing to analyse cytoplasmic alignments"
  #   callWhippet || die "unable to prepare whippet index"
  # fi

  die "no further instructions provided"
}
