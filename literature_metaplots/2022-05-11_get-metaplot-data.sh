#!/bin/bash

# written by AJ Sethi on 2022-05-11
# Download metaplot data for nanograd folder

# Study	SRA number
# Nanopore native RNA sequencing of a human poly(A) transcriptome	https://github.com/nanopore-wgs-consortium/NA12878
# Accurate expression quantification from nanopore direct RNA sequencing with NanoCount	PRJEB39347.
# Identification of differential RNA modifications from nanopore direct RNA sequencing with xPore	PRJEB40872
# Nanopore direct RNA sequencing detects DUX4-activated repeats and isoforms in human muscle cells	PRJDB8318
# Direct RNA sequencing enables m6A detection in endogenous transcript isoforms at base-specific resolution	PRJNA549597

############################################################
############################################################

# make a data directory
export data="/g/data/lf10/as7425/nanograd/data/2022-05-11_metaplot_data"; mkdir -p ${data} 2>/dev/null

# get NA12878 direct RNA fastq data
wget https://s3.amazonaws.com/nanopore-human-wgs/rna/fastq/NA12878-DirectRNA_All_Guppy_4.2.2.fastq.gz

# get data from SRA
module load aspera

ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/DRR178/DRR178490/DRR178490_1.fastq.gz . && mv DRR178490_1.fastq.gz DRR178490_PromethION_sequencing_of_SAMD00171761_1.fastq.gz &
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR435/004/ERR4352444/ERR4352444.fastq.gz . && mv ERR4352444.fastq.gz ERR4352444_MinION_sequencing.fastq.gz &
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR435/003/ERR4352443/ERR4352443.fastq.gz . && mv ERR4352443.fastq.gz ERR4352443_MinION_sequencing.fastq.gz &
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/DRR178/DRR178487/DRR178487_1.fastq.gz . && mv DRR178487_1.fastq.gz DRR178487_PromethION_sequencing_of_SAMD00171758_1.fastq.gz &
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR964/003/SRR9646143/SRR9646143_1.fastq.gz . && mv SRR9646143_1.fastq.gz SRR9646143_GSM3897647_HMEC_WT_Homo_sapiens_RNA-Seq_1.fastq.gz &
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR964/001/SRR9646141/SRR9646141_1.fastq.gz . && mv SRR9646141_1.fastq.gz SRR9646141_GSM3897645_HEK293_WT_Homo_sapiens_RNA-Seq_1.fastq.gz &
