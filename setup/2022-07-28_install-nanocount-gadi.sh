#!/bin/bash

# written by AJ Sethi on 2022-07-28
# installs nanocount on gadi

# go to apps
cd /g/data/lf10/as7425/apps

##################################################
##################################################

# install python dependencies using conda
# follow the tutorial at https://www.activestate.com/resources/quick-reads/how-to-manage-python-dependencies-with-conda/

conda update conda --all
conda update anaconda

conda create -n NanoCount python=3.6
conda activate NanoCount
