# About


- Nanopore direct RNA sequencing suffers from a decay in coverage towards the 5' ends of transcripts
- Nanograd performs annotation-free identification of poly(A) clusters and their associated 3’→5’ “sequencing decay” from nanopore direct-RNA sequencing alignments.

# Dependencies
Nanograd relies on recent version of:
- bedtools
- samtools
- gnu parallel
- python3 + numpy
- R + tidyverse

# usage
sh nanograd.sh [run mode] [options] -a "/path/to/fasta/directory" -b "/path/to/bam/directory/" -o "/path/to/output/directory"

# run modes:
 - *cluster*: perform annotation-free identification of poly(A) sites (see output)
 - decay: (cluster), with added information on the decay coefficient for each cluster (see output)

# options
-t [int] threadcount (default: 8)                   
-c [int] cluster confidence level (default: 25) # partially implemented        
-v nanograd verbose mode # partially implemented

# output
- a bed-like register of cluster position, cluster support level, high-confidence cluster length (to be clarufied), and decay constant      
    col 1: chromosome           
    col 2: poly(A) cluster 3' position start        
    col 3: poly(A) cluster 3' position end      
    col 4: strand               
    col 5: cluster score (i.e. the number of supporting reads)      
    col 6: cluster maximum spliced length (the maximum length in base-pairs for bases in the cluster with more than 5 supporting reads)     
    col 7: [decay mode only] the average decay for transcripts in the cluster, measured in (total reads at position)(base-pairs)^1 along the length of the cluster      

https://twitter.com/kylemorgenstein/status/1338551698332258304?s=21

# notes 

- The transcriptome annotation (GTF) must contain both 'transcript' and 'exon' entries (field 3), otherwise nanograd will calculate transcript length as including introns
