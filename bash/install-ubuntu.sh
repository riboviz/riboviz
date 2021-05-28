#!/usr/bin/env bash

sudo apt update -y
sudo apt-get install -y git
sudo apt-get install -y curl
sudo apt-get install -y emacs
sudo apt-get install -y bedtools
bedtools -version
sudo apt-get install -y hdf5-tools
sudo apt-get install -y pigz
pigz --version
sudo apt-get install -y graphviz
dot -V
sudo apt-get install -y zip
zip -v
sudo apt-get install -y unzip
unzip -v
