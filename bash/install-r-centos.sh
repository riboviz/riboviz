#!/usr/bin/env bash

sudo yum update -y
sudo yum install -y R
sudo yum install -y R-devel
sudo yum install -y libxml2-devel
sudo yum install -y openssl-devel
sudo yum install -y libcurl-devel

R --version
