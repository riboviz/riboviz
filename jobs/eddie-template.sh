#!/bin/sh
# Grid Engine options (lines prefixed with #$)
# Job name:
#$ -N %%job_name%%
# Use the current working directory:
#$ -cwd
# Runtime limit:
#$ -l h_rt=%%runtime%%
# RAM
#$ -l h_vmem=%%memory%%
# Use shared memory parallel environment and request number of CPUs:
#$ -pe sharedmem %%num_cpus%%
# Redirected output file name format:
#$ -o $JOB_NAME-$JOB_ID-$HOSTNAME.o
# Redirected error file name format:
#$ -e $JOB_NAME-$JOB_ID-$HOSTNAME.e

# Initialise the environment modules.
. /etc/profile.d/modules.sh

#!/usr/bin/env bash
export R_LIBS=%%r_libs%%
module load openmpi
module load igmm/apps/BEDTools 
module load igmm/apps/bowtie
module load igmm/apps/hdf5
module load igmm/apps/HISAT2
module load igmm/apps/pigz
module load igmm/apps/R/3.6.3
module load anaconda
source activate riboviz

DATAFOLDER="%%data_folder%%"

echo "Running riboviz on dataset: ${DATAFOLDER}"

# Move to scratch space, create and move to $DATAFOLDER.
cd /exports/eddie/scratch/$USER
mkdir -p $DATAFOLDER/input
cd $DATAFOLDER/input
echo "${PWD}"

## Get the dataset read files.
module load igmm/apps/sratoolkit/2.10.8
## Pre-fetch with Aspera client.
%%pre_fetch%%
# Use pigz to zip .fastq files into .fastq.gz files.
pigz *.fastq

echo "(Hopefully) downloaded and pigz'd files into /exports/eddie/scratch/$USER/${DATAFOLDER}/input"

# Following presumes SRA files were downloaded and they are in
# /exports/eddie/scratch/$USER/$DATAFOLDER/input/.

# Move to riboviz folder.
cd $HOME/riboviz/riboviz
echo "Moved to ${PWD}"
# Make system link in riboviz/ scratch folder.
ln -s /exports/eddie/scratch/$USER/$DATAFOLDER
# Copy YAML configuration file into riboviz/$DATAFOLDER/.
cp %%config_path%% $DATAFOLDER/

# Following presumes we are in correct branch on riboviz.
# Move back up to riboviz folder (or Nextflow won't run).
cd $HOME/riboviz/riboviz
echo "Moved to ${PWD}. Ready to run!"

# Following presumes YAML configuration file exists in
# $HOME/riboviz/riboviz/${DATAFOLDER} AND it points to
# ${DATAFOLDER}/input. for files as required.

# Run Nextflow validation.
nextflow run prep_riboviz.nf -params-file ${DATAFOLDER}/%%config_file%% -work-dir ${DATAFOLDER}/work -ansi-log false --validate_only

# Run Nextflow.
nextflow run prep_riboviz.nf -params-file ${DATAFOLDER}/%%config_file%% -work-dir ${DATAFOLDER}/work -ansi-log false

# (Hopefully) success. 
echo "nextflow riboviz ${DATAFOLDER} data run complete!"
