#!/usr/bin/env bash

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
bash miniconda3.sh -b -p $HOME/miniconda3

source $HOME/miniconda3/bin/activate
conda install -y pyyaml
conda install -y pytest
conda install -y pytest-cov
conda install -y pylint
conda install -y pycodestyle
conda install -y gitpython
conda install -y pandas
conda install -y -c bioconda cutadapt
conda install -y -c bioconda pysam
conda install -y -c anaconda biopython
pip install gffutils
conda install -y -c anaconda h5py
conda install -y -c bioconda umi_tools
conda list
source $HOME/miniconda3/bin/deactivate
