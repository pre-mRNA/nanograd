#!/bin/bash

# written by Bhavika Kumar on 2022-05-03
# aim: install uLTRA on Gadi into /g/data/xc17/bk9031/apps/ultra and set up a conda environment for uLTRA

############################################################
############################################################

# install uLTRA from https://github.com/ksahlin/ultra#INSTALLATION
# requires python=3
# requires MEM finder slaMEM and aligner minimap2

# create a conda environment for ultra
conda create -n ultra python==3.7

# activating uLTRA
conda activate ultra

# check which pip we are using (it should be from our new environment)
which pip

# install uLTRA
pip install ultra-bioinformatics

# make a directory to install slaMEM for uLTRA
mkdir -p /g/data/lf10/as7425/apps/ultra
cd /g/data/lf10/as7425/apps/ultra

# install MEM finder slaMEM
git clone https://github.com/fjdf/slaMEM.git
cd slaMEM
make
export PATH=$PATH:$(pwd)

# install aligner minimap2
conda install -c bioconda minimap2

# install strobemap
cd /g/data/lf10/as7425/apps/ultra
wget https://github.com/ksahlin/strobemers/raw/main/strobemers_cpp/binaries/Linux/StrobeMap-0.0.2
mv StrobeMap-0.0.2 StrobeMap
chmod +x StrobeMap
export PATH=$PATH:$(pwd)

# install mummer
conda install --yes -c bioconda mummer

# install R
conda config --add channels conda-forge r
conda install -c r r-essentials

conda install -c conda-forge r-base
conda upgrade -c conda-forge readline

# check if uLTRA is installed
uLTRA --help

# appears to work
