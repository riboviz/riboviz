#!/usr/bin/env bash

sudo yum update -y
sudo yum install -y git
git --version
sudo yum install -y curl
curl --version
sudo yum install -y epel-release
sudo yum install -y BEDTools
bedtools -version
sudo yum install -y hdf5-devel
h5diff --version
sudo yum install -y pigz
pigz --version
sudo yum install -y pandoc
pandoc --version
sudo yum install -y graphviz
dot -V
zip -v
unzip -v
sudo yum install -y libxml2-devel
sudo yum install -y openssl-devel
sudo yum install -y libcurl-devel
sudo yum install -y R
sudo yum install -y R-devel
R --version
