#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <sam.h>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <bam_file>\n", argv[0]);
    exit(1);
  }

  char *bam_file = argv[1];

  // Open the BAM file for reading
  samFile *in = sam_open(bam_file, "r");
  if (in == NULL) {
    fprintf(stderr, "Error opening BAM file: %s\n", bam_file);
    exit(1);
  }

  // Read the header from the BAM file
  bam_hdr_t *header = sam_hdr_read(in);
  if (header == NULL) {
    fprintf(stderr, "Error reading header from BAM file: %s\n", bam_file);
    exit(1);
  }

  // Allocate memory for the alignment record
  bam1_t *aln = bam_init1();

  // Read each alignment record from the BAM file
  while (sam_read1(in, header, aln) >= 0) {
    // Check if this is a primary alignment
    if (!(aln->core.flag & BAM_FSECONDARY)) {
      // Get the cigar string for this alignment
      uint32_t *cigar = bam_get_cigar(aln);
      int n_cigar = aln->core.n_cigar;

      // Loop through the cigar operations for this alignment
      int n_aligned_bases = 0;
      for (int i = 0; i < n_cigar; i++) {
        // Get the length and operation for this cigar element
        int op_len = cigar[i] >> BAM_CIGAR_SHIFT;
        int op = cigar[i] & BAM_CIGAR_MASK;

        // Check if this is a match / alignment operation
        if (op == BAM_CMATCH) {
          // Increment the number of aligned bases
          n_aligned_bases += op_len;
        }
      }

      // Print the read id, reference sequence name, and number of aligned bases
      printf("%s\t%s\t%d\n", bam_get_qname(aln), header->target_name[aln->core.tid], n_aligned_bases);
    }
  }

  // Clean up
  bam_destroy1(aln);
  bam_hdr_destroy(header);
  sam_close(in);

  return 0;
}
