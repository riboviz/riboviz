#!/usr/bin/env bash

sudo apt update -y
sudo apt install -y git
git --version
sudo apt install -y curl
curl --version
sudo apt install -y bedtools
bedtools -version
sudo apt install -y hdf5-tools
h5diff --version
sudo apt install -y pigz
pigz --version
sudo apt install -y pandoc
pandoc --version
sudo apt install -y graphviz
dot -V
sudo apt install -y zip
zip -v
sudo apt install -y unzip
unzip -v
sudo apt install -y libxml2-dev
sudo apt install -y libssl-dev
sudo apt install -y libcurl4-openssl-dev
sudo apt install -y libgit2-dev
sudo apt install -y r-base
sudo apt install -y r-base-dev
R --version
