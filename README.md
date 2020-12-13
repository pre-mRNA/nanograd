# nanograd

- nanograd performs annotation-free identification of polyA sites from nanopore direct-RNA sequencing alignments
- nanograd computes coverage decay (average loss in read coverage per base pair) upon each polyA cluster

# dependancies
recent version of:
- bedtools
- samtools
- gnu parallel
- python3, numpy
- R, tidyverse
- picard tools


# usage
sh nanograd.sh [mode] [options] -a "/path/to/fasta/directory" -b "/path/to/bam/" -o "/path/to/output/directory"

# run modes
 - cluster: perform annotation-free identification of poly(A) sites (see output)
 - decay: (cluster), along with information on the decay for each cluster (see output)


# options
-t (threadcount) # currently working               
-v (verbose mode) # partially implemented                 
-c (cluster confidece level, 25 by default) # partially implemented        

# output
- a bed-like register of cluster position, cluster support level, high-confidence cluster length (to be clarufied), and decay constant      
    col 1: chromosome       
    col 2: poly(A) cluster 3' position start
    col 3: poly(A) cluster 3' position end
    col 4: strand       
    col 5: cluster score (i.e. the number of supporting reads)
    col 6: cluster maximum spliced length (the number of base-pairs in endogenous reads with more than 5 reads)
    col 7: [decay mode only] the average decay for transcripts in the cluster, measured in (total reads at position)(base-pairs)^1 along the length of the cluster
