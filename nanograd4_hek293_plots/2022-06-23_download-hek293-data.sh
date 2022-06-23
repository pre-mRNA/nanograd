#!/bin/bash

# download published direct RNA sequencing data from HEK293
# Written by A.J. Sethi on 2022-05-23

##################################################
##################################################

# make a data directory
export hek293_public_data="/g/data/lf10/as7425/nanograd/data/2022-06-23_hek293-nanograd4-analysis"
# get data from Soneson 2019
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/006/ERR3218376/ERR3218376.fastq.gz || echo "cannot get soneson_1" &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/007/ERR3218377/ERR3218377.fastq.gz || echo "cannot get soneson_2" &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/008/ERR3218378/ERR3218378.fastq.gz || echo "cannot get soneson_3" &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/009/ERR3218379/ERR3218379.fastq.gz || echo "cannot get soneson_4" &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/000/ERR3218380/ERR3218380.fastq.gz || echo "cannot get soneson_5" &


# get data from Soneson 2019
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/006/ERR3218376/ERR3218376.fastq.gz || echo "cannot get soneson_1" &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/007/ERR3218377/ERR3218377.fastq.gz || echo "cannot get soneson_2" &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/008/ERR3218378/ERR3218378.fastq.gz || echo "cannot get soneson_3" &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/009/ERR3218379/ERR3218379.fastq.gz || echo "cannot get soneson_4" &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/000/ERR3218380/ERR3218380.fastq.gz || echo "cannot get soneson_5" &
