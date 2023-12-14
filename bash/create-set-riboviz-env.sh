#!/usr/bin/env bash
echo "#"'!'"/usr/bin/env bash" > set-riboviz-env.sh
echo 'export PATH=$HOME/hisat2-2.1.0:$PATH' >> set-riboviz-env.sh
echo 'export PATH=$HOME/bowtie-1.2.2-linux-x86_64/:$PATH' >> set-riboviz-env.sh
echo 'source $HOME/miniconda3/bin/activate' >> set-riboviz-env.sh
echo 'conda activate riboviz' >> set-riboviz-env.sh
source set-riboviz-env.sh
