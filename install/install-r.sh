#!/usr/bin/env bash

sudo echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" | sudo tee -a /etc/apt/sources.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update
sudo apt-get install -y r-base
sudo apt-get install -y r-base-dev
sudo apt-get install -y libxml2-dev
sudo apt-get install -y libcurl4-openssl-dev

R --version
