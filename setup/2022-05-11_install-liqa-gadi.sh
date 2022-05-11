#!/bin/bash

# written by AJ Sethi on 2022-05-11
# aim: install LIQA on Gadi into /g/data/ and set up a conda environment for LIQA

############################################################
############################################################

# install LIQA from https://github.com/WGLab/LIQA/blob/master/doc/Install.md
# requires python=3.7 and R=3.5.2, and modules
# Python: Pysam, Numpy, Scipy, Lifelines
# R:gcmr, betareg

# create a conda environmetn for liqa with python=3.7
conda create -n liqa python==3.7 && conda activate liqa

# install R/3.5.2 in conda
# conda config --add channels r
# conda install r:r-essentials=3.5.2

# for now, use gadi's module load R/3.6.1
module load R/3.6.1

# check which pip we are using (it should be from our new enironment)
which pip

# install LIQA
pip install liqa

# test liqa
liqa

# runs succesfully
