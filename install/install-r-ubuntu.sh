#!/usr/bin/env bash

sudo apt-get update -y
sudo apt-get install -y r-base
sudo apt-get install -y r-base-dev
sudo apt-get install -y libxml2-dev
sudo apt-get install -y libssl-dev
sudo apt-get install -y libcurl4-openssl-dev

R --version
