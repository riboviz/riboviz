#!/usr/bin/env bash

wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda2.sh
bash miniconda2.sh -b -p $HOME/miniconda2

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
bash miniconda3.sh -b -p $HOME/miniconda3

source $HOME/miniconda2/bin/activate
conda install -y -c bioconda cutadapt 
conda install -y -c bioconda pysam
conda list
source $HOME/miniconda2/bin/deactivate

source $HOME/miniconda3/bin/activate
conda install -y -c bioconda cutadapt 
conda install -y -c bioconda pysam
conda list
source $HOME/miniconda3/bin/deactivate
