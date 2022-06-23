#!/bin/bash

# download published direct RNA sequencing data from HEK293
# Written by A.J. Sethi on 2022-05-23

# refer to 'HEK293 Studies BK AS on Nanograd ANU OneDrive account'
# 2022-06-23 at https://anu365-my.sharepoint.com/:x:/r/personal/u6081208_anu_edu_au/_layouts/15/Doc.aspx?sourcedoc=%7BCFA6BAB3-538E-4999-869F-523EB1C6BDFD%7D&file=HEK293%20studies%20BK%20AS%20.xlsx&action=default&mobileredirect=true

##################################################
##################################################

# make a data directory
export hek293_public_data="/g/data/lf10/as7425/nanograd/data/2022-06-23_hek293-nanograd4-analysis"
cd ${hek293_public_data}

# get data from Soneson 2019
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/006/ERR3218376/ERR3218376.fastq.gz && mv ERR3218376.fastq.gz soneson_ERR3218376.fastq.gz || echo "cannot get soneson_1" &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/007/ERR3218377/ERR3218377.fastq.gz && mv ERR3218377.fastq.gz soneson_ERR3218377.fastq.gz || echo "cannot get soneson_2" &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/008/ERR3218378/ERR3218378.fastq.gz && mv ERR3218378.fastq.gz soneson_ERR3218378.fastq.gz || echo "cannot get soneson_3" &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/009/ERR3218379/ERR3218379.fastq.gz && mv ERR3218379.fastq.gz soneson_ERR3218379.fastq.gz || echo "cannot get soneson_4" &
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR321/000/ERR3218380/ERR3218380.fastq.gz && mv ERR3218380.fastq.gz soneson_ERR3218380.fastq.gz || echo "cannot get soneson_5" &

# get data from Hassan 2022
module load aspera
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR132/094/SRR13261194/SRR13261194_1.fastq.gz . && mv SRR13261194_1.fastq.gz hassan_SRR13261194.fastq.gz || echo "cannot get Hassan 1" &
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR132/095/SRR13261195/SRR13261195.fastq.gz . && mv SRR13261195.fastq.gz hassan_SRR13261195.fastq.gz || echo "cannot get Hassan 2" &
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR132/096/SRR13261196/SRR13261196_1.fastq.gz . && mv SRR13261196_1.fastq.gz hassan_SRR13261196.fastq.gz || echo "cannot get Hassan 3" &

# get data from Lorenz 2022
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR964/001/SRR9646141/SRR9646141_1.fastq.gz . && mv SRR9646141_1.fastq.gz lorenz_SRR9646141.fastq.gz || echo "cannot get Lorenz 1" &
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR964/002/SRR9646142/SRR9646142_1.fastq.gz . && mv SRR9646142_1.fastq.gz lorenz_SRR9646142.fastq.gz || echo "cannot get Lorenz 2" &

# copy over our data from Kat and Agin
cp /g/data/xc17/degradation_project/Mg_degraded/basecalled/undegraded_hek293_pass1/all.undegraded_hek293_pass1.fastq.gz ${hek293_public_data}/wt_rep1.fastq.gz
cp /g/data/xc17/degradation_project/Mg_degraded/basecalled/undegraded_hek293_pass2/all.undegraded_hek293_pass2.fastq.gz ${hek293_public_data}/wt_rep2.fastq.gz
cp /g/data/xc17/degradation_project/Mg_degraded/basecalled/5mM_MgCl_degrdation_pass1/all.5mM_MgCl_degrdation_pass1.fastq.gz ${hek293_public_data}/deg_rep1.fastq.gz
cp /g/data/xc17/degradation_project/Mg_degraded/basecalled/5mM_MgCl_degrdation_pass2/all.5mM_MgCl_degrdation_pass2.fastq.gz ${hek293_public_data}/deg_rep2.fastq.gz
