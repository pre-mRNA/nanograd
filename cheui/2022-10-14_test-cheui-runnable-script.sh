#!/bin/bash

# link script
export pipeline="/home/150/as7425/nanograd/cheui/2022-10-12_cheui-runnable-script.sh"
export data="/home/150/as7425/nanograd/cheui/CHEUI/test/nanopolish_output_test.txt"
export output="/scratch/lf10/as7425/2022-10-14_test-cheui-pipeline"
export annotation=""

# test for m6A
bash ${pipeline} "${data}" "${output}/m6A" "test_data" "A"

# test for m5C
bash ${pipeline} "${data}" "${output}/m5C" "test_data" "C"
