#!/usr/bin/env bash
#
# Print version information for the current operating system,
# command-line tools, Python packages and R packages required
# by RiboViz.
#
# Run this script in the "riboviz" home directory, as it uses
# rscripts/list-r-packages.R. For example:
#
#     $ source bash/environment.sh

echo "# RiboViz environment."
echo "# lsb_release -a"
lsb_release -a
echo "# git --version"
git --version
echo "# curl --version"
curl --version
echo "# bedtools --version"
bedtools --version
echo "# h5diff --version # Tool from hdf5tools"
h5diff --version
echo "# pigz --version"
pigz --version
echo "# pandoc --version"
pandoc --version
echo "# dot -v # Tool from GraphViz"
dot -V
echo "# zip -h"
zip -h | head -n 2
echo "# unzip -h"
unzip -h | head -n 1
echo "# python --version"
python --version
echo "# conda list"
conda list
echo "# pip list"
pip list
echo "# cutadapt --version"
cutadapt --version
echo "# samtools --version"
samtools --version
echo "# umi_tools -v"
umi_tools -v
echo "# javac -version"
javac -version
echo "# java -version"
java -version
echo "# nextflow -v"
nextflow -v
echo "# hisat2-build --version"
hisat2-build --version
echo "# hisat2 --version"
hisat2 --version
echo "# bowtie --version"
bowtie --version
echo "# R --version"
R --version
echo "# Rscript rscripts/list-r-packages.R"
Rscript rscripts/list-r-packages.R
