#!/usr/bin/env bash

sudo yum install -y git
sudo yum install -y curl
sudo yum install -y emacs
sudo yum install -y epel-release
sudo yum install -y BEDTools
bedtools -version
sudo yum install -y hdf5-devel
sudo yum install -y pigz
pigz --version
