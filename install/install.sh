#!/usr/bin/env bash

sudo apt-get install -y git
sudo apt-get install -y emacs
sudo apt-get install -y bedtools
bedtools -version
sudo apt-get install -y samtools
samtools --version
sudo apt-get install -y hdf5-tools

wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip -O hisat2-2.1.0.zip
unzip hisat2-2.1.0.zip -d ~

wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/bowtie-1.2.2-linux-x86_64.zip/download -O bowtie-1.2.2.zip
unzip bowtie-1.2.2.zip -d ~

echo "#"'!'"/usr/bin/env bash" > setenv.sh
echo 'export PATH=~/hisat2-2.1.0:$PATH' >> setenv.sh
echo 'export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH' >> setenv.sh

source setenv.sh

hisat2 --version
bowtie --version

echo $PATH
