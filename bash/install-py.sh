#!/usr/bin/env bash

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
bash miniconda3.sh -b -p $HOME/miniconda3

source $HOME/miniconda3/bin/activate
conda create -y --name riboviz python=3.7
conda activate riboviz
conda install -y pyyaml
conda install -y gitpython
conda install -y pytest
conda install -y pandas
conda install -y -c bioconda cutadapt
cutadapt --version
conda install -y -c bioconda pysam
conda install -y -c bioconda samtools=1.9
samtools --version
conda install -y -c anaconda biopython
pip install gffutils
conda install -y -c anaconda h5py
conda install -y -c bioconda umi_tools
umi_tools -v
conda install -y -c bioconda nextflow=20
javac -version
java -version
nextflow -v
nextflow -version
# Developer dependencies
conda install -y pytest-cov
conda install -y pylint
conda install -y pycodestyle
pip install sphinx
