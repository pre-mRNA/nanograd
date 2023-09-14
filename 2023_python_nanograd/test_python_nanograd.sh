# undegraded pass 1 bam
export bam="/g/data/xc17/as7425/sharing/2023-06-29_nanograd-python/all.undegraded_hek293_pass1.fastq.gz.sorted.bam"

# annotation 
export annotation="/g/data/xc17/as7425/sharing/2023-06-29_nanograd-python/gencode.v38.annotation.gtf"

# sequencing summary 
export ss="/g/data/xc17/as7425/sharing/2023-06-29_nanograd-python/sequencing_summary_FAQ86281_15c37cc7.txt"

# output dir 
export out="/g/data/xc17/as7425/sharing/2023-06-29_nanograd-python/test"
mkdir ${out} 2>/dev/null 

# run
time python3 ~/nanograd/2023_python_nanograd/2023-06-29_nanograd-python.py ${bam} ${annotation} ${ss} "${out}/test.txt" -t 48 -v 2 -k -m din 

# preprocess data (-k to keep temporary files, -v 2 for high verbosity)
python3 ~/nanograd/2023_python_nanograd/2023-06-29_nanograd-python.py ${bam} ${annotation} ${ss} "${out}/test.txt" -k -v 2 -t 48 

# calculate DIN with maximum verbosity 
python3 ~/nanograd/2023_python_nanograd/2023-06-29_nanograd-python.py ${bam} ${annotation} ${ss} "${out}/test.txt" -k -v 2 -m din 

# calculate cuts per nucleotide (not implemented with medium verbosity)
python3 ~/nanograd/2023_python_nanograd/2023-06-29_nanograd-python.py ${bam} ${annotation} ${ss} "${out}/test.txt" -k -v 1 -m cuts 

####################

# test new masking code 
export bam="/g/data/xc17/as7425/sharing/2023-06-29_nanograd-python/all.undegraded_hek293_pass1.fastq.gz.sorted.bam" # undegraded pass 1 bam
export annotation="/g/data/xc17/as7425/sharing/2023-06-29_nanograd-python/gencode.v38.annotation.gtf" # annotation 
export ss="/g/data/xc17/as7425/sharing/2023-06-29_nanograd-python/sequencing_summary_FAQ86281_15c37cc7.txt" # sequencing summary 
export out="/g/data/xc17/as7425/sharing/2023-09-12_nanograd-python/test_masking/"; mkdir -p ${out} 2>/dev/null # output dir 

# go to /g/data/ to prevent temp files from accumulating in /home/
# gtf index is produced in working directory rather than in the parent folder of the 4th argument...
cd $out 

# preprocess data (-k to keep temporary files, -v 2 for high verbosity)
# note, the 4th argument doesn't do anything in this mode, since we're not returning any specific output from nanograd 
# the preprocess table is saved as a temp file in the output directory due to use of the -k (keep temp files) flag 
python3 ~/nanograd/2023_python_nanograd/2023-09-12_nanograd-python-censored.py ${bam} ${annotation} ${ss} "${out}/generate_preprocess.txt" -k -t 48 

# mask the temporary file 
# note, the temporary prefix is random and needs to be accomodated in this step for --input_file
python3 ~/nanograd/2023_python_nanograd/2023-09-12_nanograd-python-censored.py -m mask --input_file ${out}/*_merged_data.tsv --output_file_mask "${out}/merged_data_censored.txt"

# test output 
cd $out 
awk -F '\t' 'NR==1 || $16 != $15' *censored* | wc -l
awk -F '\t' 'NR==1 || $16 != $15' *censored* | grep "unblock" | wc -l



####################

# 2023-09-14
# implement Alice's saturatino code into the 'cuts' mode of nanograd 

# in R, function per transcript-group 
# Saturation<-function(y)
#     return( which.max(diff(cumsum(sapply(1:max(y),function(x) {sum(y==x)})),lag=window))+window ))

# in R, to apply this function: 
# Saturation_df = grouped_by_transcript["mdi_read_length"].apply(lambda x:saturation(x)).reset_index(name="Saturation_3p")

# in python, we can run this in one step on the `merge_temp_file` object, which looks like this: 

# read_id ref_transcript  sequence_length 5_softclip      3_softclip      cigar_match     cigar_ins       cigar_del       cigar_splice    cigar_softclip_total    cigar_hardclip_total    end_reason      transcript_annotated_length     transcript_biotype      mdi_read_length
# 8726b371-ee82-486d-b346-2f11aab03213    ENST00000361390 1146    258     79      796     13      160     0       337     0       signal_positive 956     protein_coding
# 8726b371-ee82-486d-b346-2f11aab03213    ENST00000361390 1146    258     79      796     13      160     0       337     0       signal_positive 956     protein_coding  969
# a0885835-c178-4087-bcd3-5e9604a97c15    ENST00000361390 1052    111     39      890     12      66      0       150     0       signal_positive 956     protein_coding

# test new masking code 
export bam="/g/data/xc17/as7425/sharing/2023-06-29_nanograd-python/all.undegraded_hek293_pass1.fastq.gz.sorted.bam" # undegraded pass 1 bam
export annotation="/g/data/xc17/as7425/sharing/2023-06-29_nanograd-python/gencode.v38.annotation.gtf" # annotation 
export ss="/g/data/xc17/as7425/sharing/2023-06-29_nanograd-python/sequencing_summary_FAQ86281_15c37cc7.txt" # sequencing summary 
export out="/g/data/xc17/as7425/sharing/2023-09-12_nanograd-python/test_cuts/"; mkdir -p ${out} 2>/dev/null # output dir 

# go to /g/data/ to prevent temp files from accumulating in /home/
# gtf index is produced in working directory rather than in the parent folder of the 4th argument...
cd $out 

# run the censor mode 
python3 ~/nanograd/2023_python_nanograd/2023-09-14_nanograd-python-censored-cuts.py ${bam} ${annotation} ${ss} "${out}/generate_preprocess.txt" -m cuts -k -t 48 

awk -F '\t' 'NR==1 || $16 != $15' *censored* | wc -l
awk -F '\t' 'NR==1 || $16 != $15' *censored* | grep "unblock" | wc -l