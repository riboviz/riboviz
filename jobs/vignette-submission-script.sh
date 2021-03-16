#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N nxfl-vig-symlk
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 16
#$ -o $JOB_NAME-$JOB_ID-$HOSTNAME.o
#$ -e $JOB_NAME-$JOB_ID-$HOSTNAME.e

#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 48 hours: -l h_rt
#  ask for 8 Gbyte RAM: -l h_vmem
#  use shared memory parallel environment, request 16 CPUs
#  redirect output with format jobname-jobID-hostname (jobname -N)
#  redirect error with same format as output
# Initialise the environment modules
. /etc/profile.d/modules.sh

#!/usr/bin/env bash
export R_LIBS=/exports/csce/eddie/biology/groups/wallace_rna/Rlibrary
module load openmpi
module load igmm/apps/BEDTools
module load igmm/apps/bowtie
module load igmm/apps/hdf5
module load igmm/apps/HISAT2
module load igmm/apps/pigz
module load igmm/apps/R/3.6.3
module load anaconda
source activate riboviz


## NOTE: this script is designed to work best with the
## following file structure:
##   $HOME/riboviz/riboviz  # (branch: develop)
##   $HOME/riboviz/example-datasets  # (branch: master)
## To recreate this:

## at $HOME:
##   mkdir riboviz
##   cd riboviz
##   git clone https://github.com/riboviz/example-datasets.git
##   cd example-datasets
##   git checkout master
##   cd ..   # back to $HOME/riboviz

## at $HOME/riboviz:
##   git clone https://github.com/riboviz/riboviz.git
##   cd riboviz
##   git checkout develop

## then check you have this file structure:
## $HOME/riboviz/riboviz  # (branch: develop)
## $HOME/riboviz/example-datasets  # (branch: master)


DATAFOLDER="vignette"

echo "Running riboviz on dataset: ${DATAFOLDER}"

# move to scratch space
cd /exports/eddie/scratch/$USER

# make folder there
mkdir remote-vignette
mkdir remote-vignette/input

# copy vignette input files there
cp -r $HOME/riboviz/riboviz/vignette/input/ remote-vignette/

# move to riboviz folder
cd $HOME/riboviz/riboviz

# make system link at riboviz folder to folder on scratch
ln -s /exports/eddie/scratch/$USER/remote-vignette

# copy the relevant yaml into the riboviz/scratch folder, rename it
cp $HOME/riboviz/riboviz/vignette/remote_vignette_config.yaml remote-vignette/local_vignette_config.yaml
# NOTE: this example uses a separate yaml configured to point to the scratch locations 
# and is a modified version of vignette/vignette_config.yaml

# This remote_vignette_config.yaml has been edited to point to remote-vignette/input and remote-vignette/index etc
# The yaml can point to files in these ways:
#  orf_fasta_file: /home/$USER/riboviz/example-datasets/fungi/saccharomyces/annotation/Saccharomyces_cer$
#  orf_gff_file: ../../riboviz/example-datasets/fungi/saccharomyces/annotation/Saccharomyces_cerevisiae_yea$

# presuming you're in correct branch on riboviz

echo "nextflow validation"

#run nextflow validation: 
nextflow run prep_riboviz.nf -params-file remote-vignette/local_vignette_config.yaml -work-dir remote-vignette/work -ansi-log false --validate_only

echo "nextflow validation complete"

echo "running nextflow"

# run nextflow: 
nextflow run prep_riboviz.nf -params-file remote-vignette/local_vignette_config.yaml -work-dir remote-vignette/work -ansi-log false

echo "nextflow run complete"

# hopefully success
echo "run_nextflow-vignette.sh complete"
