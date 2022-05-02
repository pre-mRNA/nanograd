#!/bin/bash

# written by AJ Sethi on 2022-05-02
# installs squiggleplot, and use a separate conda environment to set it up because it's good practice

# go to apps
cd /g/data/lf10/as7425/apps

# download squigglekit
git clone https://github.com/Psy-Fer/SquiggleKit.git

##################################################
##################################################

# install python dependencies using conda
# follow the tutorial at https://www.activestate.com/resources/quick-reads/how-to-manage-python-dependencies-with-conda/

conda update conda --all
conda update anaconda

conda create --name squigglekit python=3.8
source activate squigglekit

# check pip is running in the conda environment
which pip
# /g/data/lf10/as7425/apps/miniconda3/envs/squigglekit/bin/pip
# it is!

# now, install python dependencies
pip install numpy h5py sklearn matplotlib
# do this after to get around an annoying version check bug
pip install pyslow5

# required for ont_fast5_api
pip install ont_fast5_api
