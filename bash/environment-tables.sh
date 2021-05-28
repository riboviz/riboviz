#!/usr/bin/env bash
#
# Print version information for the command-line tools, Python
# packages and R packages required by RiboViz as 3 tables, each in
# Markdown format.
#
# Run this script in the "riboviz" home directory, as it uses
# rscripts/list-r-packages.R. For example:
#
#     $ source bash/environment-tables.sh
#
# Examples of the version strings output by each command are given to
# help show how these are parsed to extract version numbers.
echo "| Command-line tool | Version |"
echo "| ----------------- | ------- |"
GIT_VERSION=$(git --version | cut -d" " -f3)
# git version 2.17.1
echo "| Git | $GIT_VERSION |"
CURL_VERSION=$(curl --version | head -n 1 | cut -d" " -f2)
# curl 7.68.0 (x86_64-conda_cos6-linux-gnu) libcurl/7.68.0 OpenSSL/1.1.1 zlib/1.2.11 libssh2/1.8.2
echo "| cURL | $CURL_VERSION |"
BEDTOOLS_VERSION=$(bedtools --version | cut -d" " -f2 | tr -d 'v')
# bedtools v2.26.0
# Use tr to remove leading "v".
echo "| bedtools | $BEDTOOLS_VERSION |"
HD5TOOLS_VERSION=$(h5diff --version | cut -d" " -f3)
# h5diff: Version 1.10.4
echo "| hdf5tools (h5diff) | $HD5TOOLS_VERSION |"
PIGZ_VERSION=$(pigz --version 2>&1 | cut -d" " -f2)
# pigz 2.4
# pigz outputs version information to standard error stream.
echo "| pigz | $PIGZ_VERSION |"
PANDOC_VERSION=$(pandoc --version | head -n 1 | cut -d" " -f2)
# pandoc 1.19.2.4
echo "| pandoc | $PANDOC_VERSION |"
GRAPHVIZ_VERSION=$(dot -V 2>&1 | cut -d" " -f5)
# dot - graphviz version 2.40.1 (20161225.0304)
# dot outputs version information to standard error stream.
echo "| GraphViz (dot) | $GRAPHVIZ_VERSION |"
ZIP_VERSION=$(zip -h | head -n 2 | tail -n 1 | cut -d" " -f2)
# Copyright (c) 1990-2008 Info-ZIP - Type 'zip "-L"' for software license.
# Zip 3.0 (July 5th 2008). Usage:
echo "| zip | $ZIP_VERSION |"
UNZIP_VERSION=$(unzip -h | head -n 1 | cut -d" " -f2)
# UnZip 6.00 of 20 April 2009, by Debian. Original by Info-ZIP.
echo "| unzip | $UNZIP_VERSION |"
PYTHON_VERSION=$(python --version | cut -d" " -f2)
# Python 3.7.3
echo "| Python | $PYTHON_VERSION |"
CUTADAPT_VERSION=$(cutadapt --version)
# 1.18
echo "| Cutadapt | $CUTADAPT_VERSION |"
SAMTOOLS_VERSION=$(samtools --version | head -n 1 | cut -d" " -f2)
# samtools 1.9
echo "| samtools | $SAMTOOLS_VERSION |"
UMITOOLS_VERSION=$(umi_tools -v | cut -d" " -f3)
# UMI-tools version: 1.0.1
echo "| UMI-tools | $UMITOOLS_VERSION |"
JAVAC_VERSION=$(javac -version 2>&1 | head -n 1 | cut -d" " -f2)
# javac 1.8.0_152-release
echo "| Java (javac) | $JAVAC_VERSION |"
JAVA_VERSION=$(java -version 2>&1 | head -n 1 | cut -d" " -f3 | tr -d '"'x)
# openjdk version "1.8.0_152-release"
# java outputs version information to standard error stream.
# Use tr to remove enclosing quotes.
echo "| Java (java) | $JAVA_VERSION |"
NEXTFLOW_VERSION=$(nextflow -v | cut -d" " -f3)
# nextflow version 20.01.0.5264
echo "| Nextflow | $NEXTFLOW_VERSION |"
HISAT2_VERSION=$(hisat2 --version)
HISAT2_VERSION=$(echo "$HISAT2_VERSION" | head -n 1  | cut -d" " -f3)
# /home/ubuntu/hisat2-2.1.0/hisat2-align-s version 2.1.0
# Use two steps as combining into a single line gives:
# (ERR): hisat2-align exited with value 141
echo "| Hisat2 | $HISAT2_VERSION |"
BOWTIE_VERSION=$(bowtie --version | head -n 1 | cut -d" " -f3)
echo "| Bowtie | $BOWTIE_VERSION |"
# /home/ubuntu/bowtie-1.2.2-linux-x86_64/bowtie-align-s version 1.2.2
R_VERSION=$(R --version | head -n 1 | cut -d" " -f3)
# R version 3.4.4 (2018-03-15) -- "Someone to Lean On"
echo "| R | $R_VERSION |"
echo " "

echo "| Python Package | Version | Package Manager |"
echo "| -------------- | ------- | --------------- |"
CONDA_PKGS=$(conda list)
CONDA_LIST=(biopython cutadapt gitpython h5py nextflow pandas pycodestyle pylint pysam pytest pytest-cov pyyaml samtools sphinx umi_tools)
for pkg in ${CONDA_LIST[@]}; do
    PKG_VERSION=$(echo "$CONDA_PKGS" | grep -iw "$pkg " | tr -s " " | cut -d" " -f2)
    # pkg     M.N    ...
    # Use 'echo "$CONDA_PKGS"' not 'echo $CONDA_PKGS' to preserves
    # newlines so grep can be used.
    # grep for "$pkg " to ensure exact matches. "-iw" is not enough to
    # stop "pytest" matching "pytest-cov", for example.
    # Delimiter between columns can be 1 or more spaces, so use tr
    # to remove multiple spaces.
    echo "| $pkg | $PKG_VERSION | conda | |"
done
PIP_PKGS=$(pip list)
PIP_LIST=(gffutils sphinx)
for pkg in ${PIP_LIST[@]}; do
    PKG_VERSION=$(echo "$PIP_PKGS" | grep -iw "$pkg " | tr -s " " | cut -d" " -f2)
    echo "| $pkg | $PKG_VERSION | pip |"
done
echo " "
echo "| R Package | Version |"
echo "| --------- | ------- |"
R_PKGS=$(Rscript rscripts/list-r-packages.R)
R_LIST=(Biostrings devtools git2r glue here knitr lintr optparse plotly RcppRoll rhdf5 rmarkdown roxygen2 Rsamtools rtracklayer shiny ShortRead styler testthat tidyverse withr)
for pkg in ${R_LIST[@]}; do
    PKG_VERSION=$(echo "$R_PKGS" | grep -iw "$pkg " | tr -s " ")
    PKG_VERSION=$(echo $PKG_VERSION | cut -d" " -f2)
    #   pkg     M.N
    # Delimiter between columns can be 1 or more spaces, so tr is used
    # to remove multiple spaces.
    # To remove leading space, use echo $PKG_VERSION outwith a string.
    echo "| $pkg | $PKG_VERSION |"
done
