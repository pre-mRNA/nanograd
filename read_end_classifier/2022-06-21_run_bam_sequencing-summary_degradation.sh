# make a working directory
wd="/g/data/lf10/as7425/nanograd/analysis/2022-06-21_merge-SS-temp/"

# link the nanograd temp files
export deg_1="/g/data/lf10/as7425/nanograd/analysis/2022-04-05_test-nanograd4-decay/all.5mM_MgCl_degrdation_pass1.fastq.gz.sorted_nanograd4_out.txt.temp"
export deg_2="/g/data/lf10/as7425/nanograd/analysis/2022-04-05_test-nanograd4-decay/all.5mM_MgCl_degrdation_pass2.fastq.gz.sorted_nanograd4_out.txt.temp"
export wt_1="/g/data/lf10/as7425/nanograd/analysis/2022-04-05_test-nanograd4-decay/all.undegraded_hek293_pass1.fastq.gz.sorted_nanograd4_out.txt.temp"
export wt_2="/g/data/lf10/as7425/nanograd/analysis/2022-04-05_test-nanograd4-decay/all.undegraded_hek293_pass2.fastq.gz.sorted_nanograd4_out.txt.temp"

# link the sequencing summary files
export deg_1_ss="/g/data/xc17/degradation_project/Mg_degraded/data/5mM_MgCl_degrdation_pass1/sequencing_summary_FAR08510_cd74bc1a.txt"
export wt_1_ss="/g/data/xc17/degradation_project/Mg_degraded/data/undegraded_hek293_pass1/sequencing_summary_FAQ86281_15c37cc7.txt"
export deg_2_ss="/g/data/xc17/degradation_project/Mg_degraded/data/5mM_MgCl_degrdation_pass2/sequencing_summary_FAP73846_c4e33c9d.txt"
export wt_2_ss="/g/data/xc17/degradation_project/Mg_degraded/data/undegraded_hek293_pass2/sequencing_summary_FAP73818_93b719e9.txt"

# link the script of interest
export script="/home/150/as7425/nanograd/read_end_classifier/2022-06-16_process_bam_sequencing-summary.R"

# load R
export R_LIBS=/g/data/lf10/as7425/apps/Rlibs
module load R

# link the annotation
export anno="/g/data/lf10/as7425/genomes/human_genome/Homo_sapiens.GRCh38.104.chr.gtf"
# updated to use the full annotation rather than the transcript columns only ### 2022-06-22

# call the Rscript for each of our data
Rscript "${script}" "${anno}" "${deg_1}" "${deg_1_ss}" "deg_rep1" ${wd}/deg_rep1.txt || echo "failed for deg1" &
Rscript "${script}" "${anno}" "${deg_2}" "${deg_2_ss}" "deg_rep2" ${wd}/deg_rep2.txt || echo "failed for deg2" &
Rscript "${script}" "${anno}" "${wt_1}" "${wt_1_ss}" "wt_rep1" ${wd}/wt_rep1.txt || echo "failed for wt 1" &
Rscript "${script}" "${anno}" "${wt_2}" "${wt_2_ss}" "wt_rep2" ${wd}/wt_rep2.txt || echo "failed for wt2" &
wait && echo "done for all"


# merge the data keeping only one copy of the header
cat ${wd}/deg_rep1.txt <(cat ${wd}/deg_rep2.txt | tail -n +2) <(cat ${wd}/wt_rep1.txt | tail -n +2) <(cat ${wd}/wt_rep2.txt | tail -n +2) | gzip -c > ${wd}/all_degradation_combined.txt.gz
