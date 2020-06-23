#!/usr/bin/env bash

echo "# lsb_release -a"
lsb_release -a
echo "# python --version"
python --version
echo "# R --version"
R --version
echo "# hisat2 --version"
hisat2 --version
echo "# bowtie --version"
bowtie --version
echo "# umi_tools -v"
umi_tools -v
echo "# pip list"
pip list
echo "# conda list"
conda list
echo "# Rscript rscripts/list-r-packages.R"
Rscript rscripts/list-r-packages.R
