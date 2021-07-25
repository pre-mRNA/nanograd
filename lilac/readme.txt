# written by AJ Sethi on 2021-07-25
# the purpose of lilac is to QC alignments of nanopore direct RNA sequencing

# lilac returns two values for a given combination of reference sequence (.fasta), alignment (.bam), and reference transcriptome (.gtf)
# base alignment rate, which is defined as the ratio of sequenced bases which produce primary alignments with mapq (minimum) against the reference genome; and
# splice in reference rate (SIRR), defined as the number of splice events which perfectly correspond to events in the reference gtf
