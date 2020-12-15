#!/bin/bash

# written by A.J. Sethi on 2020-11-18
# aim: this script computes the coverage gradient for annotated transcripts using nanopore direct RNA alignments

# dependencies: recent versions of: samtools, bedtools, R, GNU parallel

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
      [ -d ${OPTARG} ] && export input="${OPTARG}" && inputFiles=`ls -d ${OPTARG}/*.*` || die "cannot access files in annotation directory" # get a list of files in the primary input
      gunzipCount=$(echo $inputFiles | grep -c -e ".gz"); if [ ${gunzipCount} -gt "0" ]; then ssec "unzipping ${gunzipCount} files in ${1}"; for i in ${1}/*; do gunzip ${i} || die "cannot unzip ${i} " & done; wait; fi # unzip any .gz files if they are present

      # ref fasta
      fastaCount=$(echo $inputFiles |  tr " " "\n" | grep -v ".fai" | grep -c -e ".fasta" -e ".fa" -e "fna"); [ ${fastaCount} -eq "1" ] || die "did not identify a single input fasta in ${1}, identified ${fastaCount}" # check that a fasta of some sort is present
      export myFasta="${input}/$(ls ${OPTARG} |  tr " " "\n" | grep -v ".fai" | grep -e ".fasta" -e ".fa" -e "fna")"

      # fasta index
      faiCount=$(echo $inputFiles |  tr " " "\n" | grep -c -e ".fai"); [ ${faiCount} -gt "1" ] && die "did not identify a single input fasta index in input, identified ${faiCount} -- try using an annotation directory with a single fasta index" # check that a fasta of some sort is present
      if [ ${faiCount} -gt "0" ]; then samtools faidx ${myFasta} > ${myFasta}.fai || die "cannot generate fasta index"; fi
      export myFai="${input}/$(ls ${OPTARG} |  tr " " "\n" | grep ".fai")"

      # genome file
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
export cl="25"
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

# define a function to split the original long-read alignemnts into ${threadCount} subsequent alignments to enable downstream parallelisation
function splitBamThreadCount() {

  # test that this function hasn't already been done -- START
  export splitBamThreadCountComplete="0"
  if [ -f ${od}/splitBamThreadCount/.complete ]; then

    # check when the completion veficication was last modified
    splitFileModDate=$(stat ${od}/splitBamThreadCount/.complete | grep "Modify" | tr " " "\t" | cut -f2,3 | cut -f1 -d".")

    # tell the user and return from the module
    ssec "Input bam has already been split and was last modified at ${splitFileModDate} --- Returning"
    export splitBamThreadCountComplete="1"
  fi
  [ ${splitBamThreadCountComplete} -eq "1" ] && return 0
  # test that this function hasn't already been done -- FINISH

  # make the output directory
  mkdir -p ${sb}/splitRaw 2>/dev/null
  [ -d ${sb}/splitRaw ] || die "cannot access sb"

  # write the bam header
  samtools view -H ${firstBam} > ${sb}/header.txt 2>/dev/null || die "cannot write header"
  headCount=$(cat ${sb}/header.txt | wc -l)
  [ ${headCount} -eq "0" ] && die "header empty"
  export header="${sb}/header.txt"

  # calculate how many subfiles to split the original bam into
  export splitCountTarget=$(python3 -c "print (round((${firstBamCount}.0 + ${threadCount}.0) / ${threadCount}.0))")

  # initially, view the BAM as a SAM and split line by line
  ssec "Splitting original alignment for parallelisation"
  cd ${sb}/splitRaw && samtools view -F 4 ${firstBam} | split -l ${splitCountTarget} -d - ${sb}/splitRaw/split || die "cannot split into n lines"

  # define a function to append the header to each split subsection
  function repairHeader() {
    cat ${header} ${1} | samtools view -b -u - > ${1%.*}.bam || die "cannot attach header for $1"
    rm ${1} || die "cannot remove $1"
  }; export -f repairHeader

  ssec "repairing BAM header"
  cd ${sb}/splitRaw && ls *split*  | parallel -j "${threadCount}" repairHeader {} || die "cannot repair header for all split sams"
  ssec "succesfully repaired BAM header"

  # mark that the function is complete
  touch ${od}/splitBamThreadCount/.complete || die "cannot finish splitbamthreadcount"
}; export -f splitBamThreadCount

# define a function to convert SAM intervals into BED intervals of last (i.e. 3') mapped cooridnate
function multiSamToEndCoordinate() {

  nsec "extracting alignment 3' end coordinates"

  # test that this function hasn't already been done -- START
  export multiSamToEndCoordinateComplete="0"
  if [ -f ${ec}/.complete ]; then

    # check when the completion veficication was last modified
    ecModDate=$(stat ${ec}/.complete | grep "Modify" | tr " " "\t" | cut -f2,3 | cut -f1 -d".")

    # tell the user and return from the module
    ssec "End coordinates already extracted and were last modified at ${ecModDate} --- Returning"
    export multiSamToEndCoordinateComplete="1"
  fi

  [ ${multiSamToEndCoordinateComplete} -eq "1" ] && return 0
  # test that this function hasn't already been done -- FINISH

  # make a subdirectory for the raw beds
  mkdir -p ${ec}/rawBed 2>/dev/null
  [ -d ${ec}/rawBed ] || die "cannot access raw bam"

  # define a subfunction to convert SAM To BED
  function multiBamToBed() {
    local bamName=${1##*/}
    local bedName=${bamName%.*}
    #echo "converting ${1}"
    bedtools bamtobed -i ${1} > ${ec}/rawBed/${bedName}.bed || die "cannot generate BAM for -"
  }; export -f multiBamToBed

  # call subfunction in parallel and convert our bam alignents to bed
  ssec "converting alignments to bed format"
  cd ${sb}/splitRaw && ls *bam  | parallel -j "${threadCount}" multiBamToBed {} || die "cannot repair header for all split sams"
  #ssec "done converting to bed"

  # collect the bedfiles and check how many records we have
  cat ${ec}/rawBed/split*bed > ${ec}/rawBed/totalAligned.bed || die "cannot collect bed files"
  totalBedCount=$(cat ${ec}/rawBed/totalAligned.bed | wc -l)
  [ ${totalBedCount} -gt "0" ] || die "no bed records found"

  # split the bedfile by strand
  mkdir ${ec}/strandBed/ 2>/dev/null
  [ -d "${ec}/strandBed/" ] || die "cannot access strandBed"
  awk -v ecp="${ec}" '{print>ecp"/strandBed/"$6}' ${ec}/rawBed/totalAligned.bed || die "cannot split bed file by strands"
  for i in ${ec}/strandBed/+ ${ec}/strandBed/-; do
    preName=${i##*/}
    name=${preName%.*}
    #echo "strand_${name}.bed"
    mv $i ${ec}/strandBed/strand_${name}.bed
  done

  # select the terminal 3' base for bed records on each strand
  ssec "selecting teminal base"
  cat ${ec}/strandBed/strand_+.bed | awk -F "\t" '{print$1,$3-1,$3,$4,$5,$6}' > ${ec}/strandBed/totalTranscriptionalEnds.bed
  cat ${ec}/strandBed/strand_-.bed | awk -F "\t" '{print$1,$2,$2+1,$4,$5,$6}' >> ${ec}/strandBed/totalTranscriptionalEnds.bed
  bedtools sort -i <(cat ${ec}/strandBed/totalTranscriptionalEnds.bed | tr " " "\t") > ${ec}/strandBed/sortedTotalTranscriptionalEnds.bed

  # count the number of 3' end coordinates
  endCount=$(cat ${ec}/strandBed/sortedTotalTranscriptionalEnds.bed | wc -l)
  [ ${endCount} -gt 0 ] || die "no reads found"
  ssec "identified ${endCount} polyA events"

  # cluster polyA sites using bedtools cluster
  # given the noise in nanopore data, polyA clusters will tend to be slightly 5' of their actual position
  ssec "clustering polyA sites"
  bedtools cluster -s -d 15 -i ${ec}/strandBed/sortedTotalTranscriptionalEnds.bed > ${ec}/clusters.txt
  clusterCount=$(cat ${ec}/clusters.txt | cut -f7 | sort -u | wc -l)
  ssec "${clusterCount} clusters were observed"

  # mark this function as complete
  touch ${ec}/.complete

}; export -f multiSamToEndCoordinate

# define a function to perform statistical analyis on the observed 3' end clusters
function clusterStatistics() {


  nsec "analysing 3' end coordinates"

  # test that this function hasn't already been done -- START
  export clusterStatisticsComplete="0"
  if [ -f ${ac}/.complete ]; then

    # check when the completion veficication was last modified
    acModDate=$(stat ${ac}/.complete | grep "Modify" | tr " " "\t" | cut -f2,3 | cut -f1 -d".")

    # tell the user and return from the module
    ssec "End coordinates already extracted and were last modified at ${acModDate} --- Returning"
    export clusterStatisticsComplete="1"

  fi

  [ ${clusterStatisticsComplete} -eq "1" ] && return 0
  # test that this function hasn't already been done -- FINISH

  # make the output directory
  mkdir -p ${ac} 2>/dev/null
  [ -d ${ac} ] || die "cannot access AC"

  # call the R script which analyses clusters
  # arguments are
    # 1: path to output (a file of high-confidence clusters)
  touch ${ac}/highConfidenceClusters.bed || die "test A"
  Rscript ${SCRIPTPATH}/scripts/cluster.R --args "${ac}/highConfidenceClusters.bed" "${ec}/clusters.txt" &>/dev/null || die "cannot run cluster.R succesfully"

  # mark this function as complete
  touch ${ac}/.complete || die "test B"

}; export -f clusterStatistics

# assign each read to a high-confidence cluster for subsequent analysis
function assignClusterToRead() {

  nsec "analysing read clusters"

  # test that this function hasn't already been done -- START
  export assignClusterToReadComplete="0"
  if [ -f ${cr}/.complete ]; then

    # check when the completion veficication was last modified
    crModDate=$(stat ${cr}/.complete | grep "Modify" | tr " " "\t" | cut -f2,3 | cut -f1 -d".")

    # tell the user and return from the module
    ssec "End coordinates already extracted and were last modified at ${crModDate} --- Returning"
    export assignClusterToReadComplete="1"
  fi

  [ "${assignClusterToReadComplete}" -eq "1" ] && return 0
  # test that this function hasn't already been done -- FINISH

  # first, read the bedtools cluster data and make a separate file containing the read names for transcipts within each cluster
  mkdir -p ${cr}/clusterTargets || die "cannot make cluster target directory"
  cd ${cr}/clusterTargets || die "cannot access cluster target directory"

  #echo "awk" && time
  awk -F "\t" '{print$4>$7}' ${ec}/clusters.txt || die "cannot split cluster file"

  # delete clusters with supporting reads less than confidence level cl;
  cd ${cr}/clusterTargets || die "cannot get to cluster targets"
  # find ${cr}/clusterTargets/ -type f -exec awk -v x=${cl} 'NR==x{exit 1}' {} \; -exec rm -f {} \;
  # ^ migrated to function fitlerCluster


  # filter for clusters with reads >= cl
  function filterCluster() {
  #echo "loop"
  local clust=$1
  local clusterCount=$(cat $clust | wc -l)
  if [ ${clusterCount} -lt "${cl}" ]; then rm $clust || die "cannot remove cluster $clust"; fi
  }; export -f filterCluster

  ssec "filtering for clusters with at least ${cl} supporting reads"
  cd ${cr}/clusterTargets/ && ls * | parallel -j ${threadCount} filterCluster {} || die "cannot filter for clusters with more than ${cl} read"

  # remove empty clusters
  cd ${cr}/clusterTargets && for i in *; do len=$(cat $i | wc -l); [ ${len} -eq "0" ] && rm ${i}; done

  # make the clustered bam directory
  export clusteredBam=${cr}/clusteredBam
  mkdir -p ${clusteredBam} 2>/dev/null
  [ -d ${clusteredBam} ] || die "cannot find ${clusteredBam}"

  # define a subfunction to iterate over each cluster and recover the reads into a new bam file
  # each iteration of the function operates over one cluster
  function recoverClusterBam() {
    clusterLong=${1##*/}
    cluster=${clusterLong%.*}

    grep -Fwf ${1} ${localSam} > ${clusteredBam}/${cluster}.sam
    samtools view -u -b <(cat ${localHeader} ${clusteredBam}/${cluster}.sam) | bedtools bamtobed -splitD -i - | bedtools bedtobam -ubam -g ${myGenome} -i - | samtools sort -l 0 > ${clusteredBam}/${cluster}.bam && rm ${clusteredBam}/${cluster}.sam || die "cannot generate ${clusteredBam}/${cluster}.bam"

    # java -jar ${picardPath}/picard.jar FilterSamReads I=${firstBam}  O=${clusteredBam}/${cluster}.bam READ_LIST_FILE=${1} SORT_ORDER=coordinate FILTER=includeReadList || die "cannot recover mates for ${1}"
  }; export -f recoverClusterBam

  # for each high-confidence cluster, recover supporting reads into a separate bam file
  ssec "aggregating alignment clusters"
  samtools view -H ${firstBam} > ${clusteredBam}/../header.txt && export localHeader="${clusteredBam}/../header.txt" || die "cannot generate sam header"
  #samtools view ${firstBam} | awk -F'\t' -v OFS='\t' '{gsub(/\N/, "D", $6)} 1' > ${clusteredBam}/../full.sam && export localSam="${clusteredBam}/../full.sam" || die "cannot generate full sam"
  #^mask gsub as we try to overcome the splicing error using bedtools

  samtools view ${firstBam} > ${clusteredBam}/../full.sam && export localSam="${clusteredBam}/../full.sam" || die "cannot generate full sam"





  cd ${cr}/clusterTargets && ls *  | parallel -j "${threadCount}" recoverClusterBam {} || die "cannot repair header for all split sams"
  #cd ${cr}/clusterTargets && ls *  | parallel -j 4 recoverClusterBam {} || die "cannot repair header for all split sams"

  export pileDir="${cr}/pileup"
  mkdir -p ${pileDir} 2>/dev/null
  [ -d ${pileDir} ] || die "cannot find ${pileDir}"

  # define a function which performs samtools mpileup for each new cluster
  function pileCluster() {
    local longBase=${1##*/}
    local clusterName=${longBase%.*}
    samtools mpileup -d 0 -f ${myFasta} -B ${1} -C 0 -Q 0 2>/dev/null | cut -f4 | awk '($1 > 8)' > ${pileDir}/${clusterName} || die "cannot perform pileup for ${clusterName}"
  }; export -f pileCluster

  # call pileCluster for each cluster
  ssec "piling reads for each cluster"
  cd ${clusteredBam} && ls *.bam | parallel -j "${threadCount}" pileCluster {} || die "cannot repair header for all split sams"

  # fetch the number of bases per cluster
  cd ${pileDir} && wc -l * | sed '$d' | awk '{print $2, $1}' | sort -k1 -n |  sed 's/.txt//' > ../clusterLength.txt

  # identify clusters with zero length
  cat ../clusterLength.txt | tr " " "\t" | awk '($2 == 0) {print $1}' > ../zeroTargets.txt
  for i in $(cat ../zeroTargets.txt); do rm ${pileDir}/${i} 2>/dev/null || die "cannot remove ${i}"; done

  # remove .DS_Store if it exists (required for mac)
  rm "${pileDir}/.DS_Store" 2>/dev/null

  # mark this function as complete
  touch ${cr}/.complete

}; export -f assignClusterToRead

# compile all the results into a bed-like file
function assembleOutput() {

  nsec "collating output"

  # test that this function hasn't already been done -- START
  export assembleOutputComplete="0"
  if [ -f ${ao}/.complete ]; then

    # check when the completion veficication was last modified
    aoModDate=$(stat ${ao}/.complete | grep "Modify" | tr " " "\t" | cut -f2,3 | cut -f1 -d".")

    # tell the user and return from the module
    ssec "End coordinates already extracted and were last modified at ${aoModDate} --- Returning"
    export assembleOutputComplete="1"
  fi

  [ "${assembleOutputComplete}" -eq "1" ] && return 0
  # test that this function hasn't already been done -- FINISH

  mkdir -p ${ao} 2>/dev/null
  [ -d "${ao}" ] || die "cannot access ${ao}"

  ssec "calculating decay gradients"
  python3 ${SCRIPTPATH}/scripts/gradient.py "${pileDir}" &>/dev/null || die "cannot calculate gradients"

  ssec "collecting output"
  Rscript ${SCRIPTPATH}/scripts/merge.R --args "${ac}/highConfidenceClusters.bed" "${cr}/output.txt" "${cr}/clusterLength.txt" "${ao}/nanograd_out.txt" && nsec "wrote output to ${ao}/nanograd_out.txt" &>/dev/null || die "cannot call merge.R"

  touch ${ao}/.complete || die "cannot make ${ao}/.complete"

}; export -f assembleOutput


main || die "cannot do main"

####################################
