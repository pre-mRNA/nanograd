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


# preprocess data (-k to keep temporary files, -v 2 for high verbosity)
python3 ~/nanograd/2023_python_nanograd/2023-06-29_nanograd-python.py ${bam} ${annotation} ${ss} "${out}/test.txt" -k -v 2 -t 48 

# mask the temporary file with 48 cores 
python3 ~/nanograd/2023_python_nanograd/2023-06-29_nanograd-python.py -m mask --input-file "${out}/tmpi1wzzfny_merged_data.tsv" --output-file "${out}/tmpi1wzzfny_merged_data_censored.tsv" -t 48 
