#!/usr/bin/env bash

mkdir sample-data
wget https://github.com/shahpr/RiboViz/blob/master/scripts/yeast_CDS_w_250utrs.fa -P sample-data
wget http://gdurl.com/KGnn/download -O sample-data/rrna.fa
wget https://github.com/shahpr/RiboViz/blob/master/scripts/yeast_CDS_w_250utrs.gff3 -P sample-data
