# written by A.J. Sethi on 2020-12-10
# splits sam records based on lookup table

# arg1 = path to sam
# arg2 = path to lookup
# arg3 = output directory

import os
import sys
from pathlib import Path

# check that inputs are provided
if len(sys.argv) < 3:
    print("Sufficient arguments not provided to gradient.py")
    exit(1)

with gzip.open("values.gz", "r") as eqtl:
    for line in eqtl.readlines():
        row = line.split()
        if row[1] in lookup:
           # assign unique file name here, then
           with open ("UniqueFileNameForThisGene","w") as outfile:
              outfile.write(line)


# time grep -Fwf "/scratch/lf10/as7425/2020-12-08_mousebrain-nanograd/assignClusterToRead/clusterTargets/10014" "/scratch/lf10/as7425/mousetest/shrek.sam"
