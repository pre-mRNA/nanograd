#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>

// Record type to hold alignment information
typedef struct {
  char* read_name;
  int alignment_length;
  char* reference_sequence_name;
} alignment_t;

// Function to free the memory used by an alignment_t record
void free_alignment(alignment_t* alignment) {
  free(alignment->read_name);
  free(alignment->reference_sequence_name);
  free(alignment);
}

int main(int argc, char** argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <input.bam> <output.txt>\n", argv[0]);
    return 1;
  }

  // Open the BAM file for reading
  samFile* in = sam_open(argv[1], "r");
  if (in == NULL) {
    fprintf(stderr, "Error: failed to open BAM file '%s'\n", argv[1]);
    return 1;
  }

  // Open the output file for writing
  FILE* out = fopen(argv[2], "w");
  if (out == NULL) {
    fprintf(stderr, "Error: failed to open output file '%s'\n", argv[2]);
    return 1;
  }

  // Read the BAM header
  bam_hdr_t* header = sam_hdr_read(in);
  if (header == NULL) {
    fprintf(stderr, "Error: failed to read BAM header\n");
    return 1;
  }

  // Create an empty alignment record
  alignment_t* alignment = malloc(sizeof(alignment_t));
  if (alignment == NULL) {
    fprintf(stderr, "Error: failed to allocate memory for alignment\n");
    return 1;
  }

  // Read each alignment from the BAM file
  bam1_t* record = bam_init1();
  while (sam_read1(in, header, record) >= 0) {
    // Check if this is a primary alignment
    if (!(record->core.flag & BAM_FPRIMARY)) continue;

    // Extract the read name
    alignment->read_name = strdup(bam_get_qname(record));
    if (alignment->read_name == NULL) {
      fprintf(stderr, "Error: failed to allocate memory for read name\n");
      return 1;
    }

    // Extract the alignment length
    alignment->alignment_length = bam_endpos(record);

    // Extract the reference sequence name
    alignment->reference_sequence_name = strdup(header->target_name[record->core.tid]);
    if (alignment->reference_sequence_name == NULL) {
      fprintf(stderr, "Error: failed to allocate memory for reference sequence name\n");
      return 1;
    }
