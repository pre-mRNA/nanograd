#!/bin/bash

# nanoGrad
# Last modified on 2022-08-02
# Written by A.J. Sethi
# This script computes the coverage gradient for each individual transcript, based on the alignemnt of nanopore direct RNA alignments to a reference CDS

# arg1 bam (transcriptome alignment)
# arg2 gtf
# arg3 output

# dependencies: recent versions of: samtools, R,
# R modules; tidyverse, runner, zoo, data.table, scales, ggpubr

# testing;
# melbourne iMac
# time bash nanograd4.sh /Users/asethi/localGadiData/2022-03-15_develop-nanograd4/allReadsTranscriptome_EnsemblnoPseudo.bam /Users/asethi/localGadiData/2022-03-15_develop-nanograd4/Mus_musculus.GRCm39.104.chr.gtf /Users/asethi/localGadiData/2022-03-15_develop-nanograd4/test_out.txt

# additional flags (to be implemented)

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
for targetMod in samtools R; do
  checkAvailable $targetMod || die "cannot locate package $targetMod in user PATH"
done

# specify the default verbosity and threadcount
export threadCount="48" # default is set to 48 for my personal usage Gadi but subject to change
export VERBOSE="false"

# check the path to nanograd
[ -n "${SCRIPTPATH+set}" ] || SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )" || die "cannot get script path"; export SCRIPTPATH="${SCRIPTPATH}" # get the script path, adapted from https://stackoverflow.com/questions/4774054/reliable-way-for-a-bash-script-to-get-the-full-path-to-itself
nsec "Nanograd is running from ${SCRIPTPATH} with parameters\n\t###\t$@\n\t###\tcommit\t$(git rev-parse --short HEAD 2>/dev/null | sed "s/\(.*\)/@\1/")"

####################################
export bam=$1
export anno=$2
export out=$3

ssec "Processing bam"
samtools view -F 256 ${bam} | cut -f1,3,6,10 > ${out}.temp || die "cannot process bam"

ssec "Running RScript"
Rscript ${SCRIPTPATH}/scripts/process_bam.R "${anno}" "${out}.temp" "${out}" || die "cannot process reads"
# rm ${out}.temp
