# nanograd

- nanograd performs transcriptome-free identification of polyA sites from nanopore direct-RNA sequencing alignments 
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
sh nanograd.sh [options] -a "/path/to/annotation (fasta and optionally gtf)" -b "/path/to/bam/directory" -o "/path/to/output/directory" 

# options 
-t (threadcount)          
-m (runmode, to be fully implemented)       
-v (verbose mode, to be fully implemented)          
-c (cluster confidece level, 25 by default, to be fully implemented)        

# output 
- a bed-like register of cluster position, cluster support level, high-confidence cluster length (to be clarufied), and decay constant      
    col 1: chromosome       
    col 2: poly(A) cluster 3' position      
    col 3: strand       
    col 4: cluster support level        
    col 5: cluster maximum spliced length       
    col 6: cluster decay        

    