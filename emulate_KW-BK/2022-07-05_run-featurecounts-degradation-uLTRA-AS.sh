export annotation="/Users/AJlocal/localGadiData/0_a_tmp/Homo_sapiens.GRCh38.104.chr.gtf"
# brew install brewsci/bio/subread
cd /Users/AJlocal/localGadiData/2022-06-22_HEK293-degradation-first4-AR_liqa-genome-alignments_BK
featureCounts --primary -L -s 1 -a ${annotation} `ls *bam` -o fc.txt
