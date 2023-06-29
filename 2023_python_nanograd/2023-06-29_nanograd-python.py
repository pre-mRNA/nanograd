#!/usr/bin/env python3

"""\
Basic python3 implementation of Nanograd4
Written by AS on 2023-06-29
aditya.sethi@anu.edu.au
"""

import argparse
import gffutils
import logging
import os
import pandas as pd
import pysam
import tempfile
import time

# set logging parameters 
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# function to parse bam :: for each primary alignment, returns read_id, sequence_length, and then length of all cigar operations 
def parse_bam(bam_file):
    start = time.time()

    # use pysam to read the bam
    with pysam.AlignmentFile(bam_file, "rb") as bam_file:
        
        # create a dict for bam data
        bam_data = {}

        # loop over records in the bam file
        for read in bam_file.fetch():
            # discard secondary/supplementary
            if read.is_secondary or read.is_supplementary:
                continue

            # basecalled sequence length
            seq_len = len(read.query_sequence)

            # get cigar operations
            cigar_ops = read.cigartuples

            # calculate cumulative length for each cigar operator
            cigar_lengths = {}
            soft_clip_left = 0
            soft_clip_right = 0
            internal_soft_clip = 0
            for i, cigartuple in enumerate(cigar_ops):
                cigar_op, cigar_len = cigartuple
                if cigar_op in cigar_lengths:
                    cigar_lengths[cigar_op] += cigar_len
                else:
                    cigar_lengths[cigar_op] = cigar_len

                # handle soft clips
                # I've quickly tested the logic here, but it would be good still to double check in individual cases 
                # the internal soft clip logic seems unnecssary, but I've left it here in case...
                if cigar_op == 4 or cigar_op == 5:  # soft clips
                    if i == 0:  # first position
                        if read.is_reverse:
                            soft_clip_right += cigar_len
                        else:
                            soft_clip_left += cigar_len
                    elif i == len(cigar_ops) - 1:  # last position
                        if read.is_reverse:
                            soft_clip_left += cigar_len
                        else:
                            soft_clip_right += cigar_len
                    else:  # internal position
                        internal_soft_clip += cigar_len


            # strip version number from transcript id
            reference_transcript = read.reference_name.split('.')[0]

            # write this data into our bam dict
            bam_data[read.query_name] = {"sequence_length": seq_len, "5' softclip": soft_clip_left,
                                         "3' softclip": soft_clip_right, "Internal softclip": internal_soft_clip,
                                         "cigar_lengths": cigar_lengths, "reference_transcript": reference_transcript}

        elapsed_time = round(time.time() - start, 2)
        logging.info(f'Parsed {len(bam_data)} lines from BAM file in {elapsed_time} seconds')

        # print a preview of bam_data if verbosity == 2
        if args.verbosity == 2:
            header = "\n".join([f"{key}: {value}" for key, value in list(bam_data.items())[:5]])
            logging.info(f"BAM Data Header:\n{header}")

        # write the bam dict to a temporary file, if the keep flag is enabled 
        temp_file_path = None
        if args.keep_temp:
            temp_file_path = tempfile.mktemp(suffix="_parsed_bam.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
            df = pd.DataFrame.from_dict(bam_data, orient="index")
            df.to_csv(temp_file_path, sep="\t", index=True)

    return bam_data, temp_file_path 

# function to filter gtf for transcripts seen in the bam (this speeds up the gtf parsing later on)
def filter_gtf(bam_data, gtf_file_path):
    start = time.time()

    # get transcript_ids of interest from the BAM we just parsed 
    transcript_ids = set([info["reference_transcript"] for info in bam_data.values()])
    temp_file_path = tempfile.mktemp(dir=os.path.dirname(os.path.abspath(gtf_file_path)))

    lines = []
    with open(gtf_file_path, 'r') as gtf_file:

        # keep cases where the transcript_id in the gtf file coincides with a transcript_id in the bam file 
        lines = [line for line in gtf_file if not line.startswith('#') and line.split('\t')[8].split(';')[1].split('"')[1].split('.')[0] in transcript_ids]

    # record how many entries we keep 
    original_line_count = sum(1 for line in open(gtf_file_path, 'r') if not line.startswith('#'))
    filtered_line_count = len(lines)

    # write the filtered lines to the temp file
    with open(temp_file_path, 'w') as temp_file:
        temp_file.writelines(lines)

    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Filtered {filtered_line_count} GTF features from {original_line_count} GTF features in {elapsed_time} seconds')

    # print header of parsed gtf if verbosity == 2 
    if args.verbosity == 2:
        with open(temp_file_path, 'r') as temp_file:
            lines = [next(temp_file) for x in range(10)]
            logging.info("First 10 lines of the filtered GTF file:\n" + "\n".join(lines))

    return temp_file_path

# function to parse gtf :: returns, transcript_id (w/o version), transcript_biotype, transcript_length 
def parse_gtf(gtf_file_path):
    start = time.time()

    # use gffutils to parse the gtf file 
    # note, there are small differences between ensembl, gencode, other gtf files 
    # this will need to be tested and refined depending on what annotation is being used 
    gtf_db = gffutils.create_db(gtf_file_path, dbfn='gtf.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True, disable_infer_genes = True, disable_infer_transcripts = True)
    
    # create a dict for transcripts 
    gtf_data = {}
    
    for feature in gtf_db.all_features(featuretype='transcript'):
        
        # remove transcript_version from transcript_id (otherwise, there can be discrecpency with the bam)
        transcript_id = feature.id.split(".")[0]  # remove version

        # fetch biotype or return NA
        # important: GENCODE uses "transcript_type" whereas other annotations use different field names to store the transcript biotype information 
        # I have seen before, transcript_biotype, gene_type etc. 
        biotype = feature.attributes['transcript_type'][0] if 'transcript_type' in feature.attributes else 'NA'
        
        # get transcript length 
        # note, we calculate trancsript length from the sum of the cumulative exons 
        # otherwise, we get the unspliced length
        length = sum(len(i) for i in gtf_db.children(feature, featuretype='exon'))  

        # store the data in a dict 
        gtf_data[transcript_id] = {"biotype": biotype, "length": length}

    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Parsed {len(gtf_data)} lines from GTF file in {elapsed_time} seconds')
    
    # print 10 lines if verbosity == 2
    if args.verbosity == 2:
        preview_data = dict(list(gtf_data.items())[:10])
        logging.info(f"First 10 lines of the GTF data:\n{preview_data}")

    # write the bam dict to a temporary file, if the keep flag is enabled 
    temp_file_path = None
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix="_parsed_gtf.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
        df = pd.DataFrame.from_dict(gtf_data, orient="index")
        df.to_csv(temp_file_path, sep="\t", index=True)
 
    return gtf_data, temp_file_path

# function to parse sequencing_summary file :: returns, read_id, read_end_reason 
def parse_summary(summary_file_path):
    start = time.time()

    # read in the data and select the read_id and read_end_reason 
    summary_df = pd.read_csv(summary_file_path, sep='\t')
    summary_data = {row["read_id"]: row["end_reason"] for _, row in summary_df.iterrows()}
    
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Parsed {len(summary_data)} lines from summary file in {elapsed_time} seconds')
    
    # print 10 lines if verbosity == 2
    if args.verbosity == 2:
        preview_data = dict(list(summary_data.items())[:10])
        logging.info(f"First 10 lines of the summary data:\n{preview_data}")

    # write the bam dict to a temporary file, if the keep flag is enabled 
    temp_file_path = None
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix="_parsed_summary.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
        df = pd.DataFrame.from_dict(summary_data, orient="index")
        df.to_csv(temp_file_path, sep="\t", index=True)

    return summary_data, temp_file_path

# function to merge data from bam, gtf, sequencing_summary, using transcript_id as key 
def merge_data(gtf_data, bam_data, summary_data):
    start = time.time()

    # header for merged data 
    merged_data = "read_id\tref_transcript\tsequence_length\t5_softclip\t3_softclip\tcigar_match\tcigar_ins\tcigar_del\tcigar_splice\tcigar_softclip_total\tcigar_hardclip_total\tend_reason\ttranscript_annotated_length\ttranscript_biotype\n"
    
    # make list to store any transcripts which have merging error with gtf data 
    gtf_error_transcripts = []
    
    # make list for transcripts which have a merging error with sequencing_summary 
    ss_error_transcripts = []

    # create list to store read names which have discrepancy between cigar lengths and read length
    cigar_discrepancy_reads = []
    cigar_discrepancy_dict = {}

    # loop over read_id in bam_data 
    for read_id, bam_info in bam_data.items():

        # select reference transcript 
        ref_transcript = bam_info["reference_transcript"]

        # as a sanity check, get total cigar length and compare to read length. They should be the same. 
        read_length = bam_info.get("sequence_length", 0) 
        
        # get cigar lengths - note, 0 = match, 1 = ins, 2 = del, 3 = splice, 4 = softclip, 5 = hardclip 
        cigar_lengths = bam_info.get("cigar_lengths", {})
        cigar_values = "\t".join(str(cigar_lengths.get(i, 0)) for i in range(6))

        # check if aligned length matches read length 
        #aligned_length = sum(cigar_lengths[0] + cigar_lengths[1] + cigar_lengths[3] + cigar_lengths[4] + cigar_lengths[5] + 1)
        aligned_length = cigar_lengths.get(0, 0) + cigar_lengths.get(1, 0) + cigar_lengths.get(3, 0) + cigar_lengths.get(4, 0) + cigar_lengths.get(5, 0) 

        # compare the CIGAR length and sequenced read length
        if aligned_length != read_length:
            cigar_discrepancy_reads.append(read_id) # if there's an error, write the read_id to our discrepency list 
            cigar_discrepancy_dict[read_id] = aligned_length - read_length # note the size of the error as well 

        # get the softclip info from the bam file 
        soft_clip_5 = bam_info.get("5' softclip", 0)
        soft_clip_3 = bam_info.get("3' softclip", 0)

        # now, check if the reference transcript in the bam is represented in gtf_data, else print to gtf_error_transcript 
        if ref_transcript in gtf_data:
            gtf_info = gtf_data[ref_transcript] # fetch the annotated transcript length and the annotated biotype 
            
            # now, check if the reference transcript is in the sequencing_summary file 
            if read_id in summary_data:
                end_reason = summary_data[read_id] # fetch the read_end_reason 
            
                # finally, combine everything into merged_data 
                merged_data += f'{read_id}\t{ref_transcript}\t{read_length}\t{soft_clip_5}\t{soft_clip_3}\t{cigar_values}\t{end_reason}\t{gtf_info["length"]}\t{gtf_info["biotype"]}\n'
            else:
                ss_error_transcripts.append(ref_transcript)
        else:
            gtf_error_transcripts.append(ref_transcript)

    # log read/transcripts pairs that have errors merging with the gtf 
    # note, we get a few in the test data, because the gtf is chromosome only whereas the FASTA that was used has scaffolds, alt contigs etc. 
    if gtf_error_transcripts:
        unique_gtf_error_transcripts = set(gtf_error_transcripts)
        with open(os.path.join(os.path.dirname(args.output_file), 'gtf_merge_errors.txt'), 'w') as f:
            f.write("\n".join(unique_gtf_error_transcripts))
        logging.error(f'Error merging data for {len(gtf_error_transcripts)} reads with gtf data. Error transcripts written to gtf_merge_errors.txt in the output directory.')
    
    # similarly, log read/transcript pairs with errors merging against the sequencing summary 
    if ss_error_transcripts:
        unique_ss_error_transcripts = set(ss_error_transcripts)
        with open(os.path.join(os.path.dirname(args.output_file), 'ss_merge_errors.txt'), 'w') as f:
            f.write("\n".join(unique_ss_error_transcripts))
        logging.error(f'Error merging data for {len(ss_error_transcripts)} reads with sequencing summary. Error transcripts written to ss_merge_errors.txt in the output directory.')

    # record cases where we have any discrepency in the cigar
    if cigar_discrepancy_reads:
        with open(os.path.join(os.path.dirname(args.output_file), '_cigar_discrepancy_errors.txt'), 'w') as f:
            for read_id, discrepancy in cigar_discrepancy_dict.items():
                f.write(f'{read_id}\t{discrepancy}\n')
        logging.error(f'Error in cigar lengths for {len(cigar_discrepancy_reads)} reads. Discrepancies written to cigar_discrepancy_errors.txt in the output directory.')
    
    # relay stats about the merging -- currently output by default, but could be moved to verbosity == 1
    logging.info(f'Total Reads Merged: {len(bam_data)}')
    logging.info(f'Total Transcripts Merged: {len(gtf_data)}')
    logging.info(f'Reads in BAM which clash with with GTF: {len(gtf_error_transcripts)}, transcripts these reads map to: {len(set(gtf_error_transcripts))}')
    logging.info(f'Reads in BAM which do not match sequencing summary: {len(ss_error_transcripts)}, transcripts these reads map to: {len(set(ss_error_transcripts))}')

    # print 10 lines of merged data if verbosity == 2
    if args.verbosity == 2:
        preview_data = merged_data.splitlines()[:10]
        preview_data_str = '\n'.join(preview_data)
        logging.info("First 10 lines of the merged data:\n%s", preview_data_str)
    
    # write merged data to tmpfile
    temp_file_path = None
    if args.keep_temp:
        temp_file_path = tempfile.mktemp(suffix="_merged_data.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
        with open(temp_file_path, 'w') as file:
            file.write(merged_data)
    logging.info(f'Merged data written to temporary file: {temp_file_path}')

    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Merged data in {elapsed_time} seconds')

    return merged_data, temp_file_path

# calculate DIN using median coverage 
def calculate_din_using_median(merged_data_file):

    start = time.time()

    # read in the data 
    merged_data = pd.read_csv(merged_data_file, sep='\t')

    # calculate read_coverage_ratio
    merged_data["read_coverage_ratio"] = merged_data["sequence_length"] / merged_data["transcript_annotated_length"]
    
    # group by ref_transcript
    grouped_by_transcript = merged_data.groupby("ref_transcript")

    # calculate the median read_coverage_ratio and total number of reads for each ref_transcript
    median_coverage_df = grouped_by_transcript["read_coverage_ratio"].median().reset_index(name="median_read_coverage_ratio")
    count_reads_df = grouped_by_transcript.size().reset_index(name="number_of_reads")

    # merge the two dataframes
    result_df = pd.merge(median_coverage_df, count_reads_df, on="ref_transcript")

    # write result_df to a temp file called TIN
    temp_file_path = tempfile.mktemp(suffix="_TIN.tsv", dir=os.path.dirname(os.path.abspath(args.output_file)))
    result_df.to_csv(temp_file_path, sep="\t", index=False)
    logging.info(f'transcript coverage data written to temporary file: {temp_file_path}')

    # calculate the overall DIN by taking the sum of TIN and dividing by the total number of transcripts
    overall_din = result_df["median_read_coverage_ratio"].sum() / len(result_df)
    logging.info(f'Overall DIN: {overall_din}')

    # update user 
    elapsed_time = round(time.time() - start, 2)
    logging.info(f'Calculated DIN in {elapsed_time} seconds')
    
    return overall_din

# calculate cuts per transcript 
def calculate_cuts_per_transcript(merged_data): 
    # TODO: Implement this function
    pass

# main 
def main(args):
    
    # logging verbosity 
    if args.verbosity == 1:
        logging.basicConfig(level=logging.INFO)
    elif args.verbosity >= 2:
        logging.basicConfig(level=logging.DEBUG)

    # parse BAM 
    logging.info('Processing BAM file')
    bam_data, bam_temp_file = parse_bam(args.bam_file)

    # filter GTF (speeds up the parsing step significantly)
    logging.info('Filtering GTF annotation')
    filtered_gtf_file = filter_gtf(bam_data, args.gtf_file)

    # parse GTF 
    logging.info('Processing GTF annotation')
    gtf_data, gtf_temp_file = parse_gtf(filtered_gtf_file)

    # parse sequencing summary file 
    logging.info('Processing basecalling sequencing summary file')
    summary_data, summary_temp_file = parse_summary(args.summary_file)

    # merge the GTF and SS data with the reads parsed from the bam file 
    logging.info('Merging data')
    merged_data, merge_temp_file = merge_data(gtf_data, bam_data, summary_data)

    # run the din calculation if din mode is selected
    if args.mode == 'din':
        logging.info('calculating sample DIN')
        calculate_din_using_median(merge_temp_file)

    # run the cuts per transcript function if mode is selected s
    elif args.mode == 'cuts':
        logging.info('calculating cuts per transcript')
        calculate_cuts_per_transcript(merge_temp_file) # to be implemented 

    # delete temp files, if keep == false (default)
    if not args.keep_temp:
        for temp_file in [bam_temp_file, gtf_temp_file, summary_temp_file, merge_temp_file]:
            if temp_file is not None and os.path.exists(temp_file):
                os.remove(temp_file)

if __name__ == "__main__":

    # parse arguments
    parser = argparse.ArgumentParser(description='Process GTF, BAM, and sequencing summary files.')
    parser.add_argument('bam_file', help='Input sorted + indexed bam file.')
    parser.add_argument('gtf_file', help='Input gtf/gff2 annnotation.')
    parser.add_argument('summary_file', help='Input sequencing summary file from guppy.')
    parser.add_argument('output_file', help='Output table.')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use (not implemented).') ## not implemented 
    parser.add_argument('-v', '--verbosity', type=int, default=0, help='Verbosity: 0 = Minimum, 1 = Information, 2 = Debugging.')
    parser.add_argument('-k', '--keep_temp', action='store_true', default=False, help='Keep temporary files in output directory.')
    parser.add_argument('-m', '--mode', type=str, default='preprocess', choices=['preprocess', 'din', 'cuts'], help='Mode of operation: preprocess (default), din, or cuts.')
    args = parser.parse_args()

    # Set temporary directory to the output file's directory
    os.environ['TMPDIR'] = os.path.dirname(os.path.abspath(args.output_file))

    main(args)