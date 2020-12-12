# written by A.J. Sethi on 2020-12-09
# caculates the gradient for each pileup text file
# provide /path/to/pileups as arg1

import os
import sys
from pathlib import Path
import numpy as np

# check that arguents are provided
if len(sys.argv) < 2:
    print("Argument was not provided to gradient.py")
    exit(1)



# assign $1 to variable $gradpath
gradPath = sys.argv[1]

# check that gradpath is a valid directory that contains files
os.path.isdir(gradPath) or print("cannot find gradpath") and exit(1)

# set working directory to gradpath
os.chdir(gradPath)

# open the output file
f = open('../output.txt','w')

# redirec the stdout to the otput file
sys.stdout = f

# iterate over files in the directory
for entry in os.scandir(gradPath):
    data = np.loadtxt(entry, delimiter=",")
    gradient = np.gradient(data) # calculate the rate of change in coverage
    #print(data_gradient_2)
    average_gradient=sum(gradient)/len(gradient) # take the average rate of change
    #print(data_gradient_2)
    #print(type(entry))
    print(os.path.splitext(os.path.basename(entry.path))[0], abs(round(average_gradient, 5))) # write the cluster name and the average gradient to the output file

f.close() # close the output file

# reset the stdout (previously going to the output file)
sys.stdout = sys.__stdout__

# print("done printing") and exit(0)
exit(0)
