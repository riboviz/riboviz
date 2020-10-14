#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N %%JOB_NAME%%
#$ -cwd
#$ -l h_rt=%%RUNTIME%
#$ -l h_vmem=%%MEMORY%%
#$ -pe sharedmem %%NUM_CPUS%%
#$ -o $JOB_NAME-$JOB_ID-$HOSTNAME.o
#$ -e $JOB_NAME-$JOB_ID-$HOSTNAME.e
#  These options are:
#  -N: job name: -N
#  -cwd: use the current working directory
#  -l h_rt: runtime limit
#  -l h_vmwm: RAM
#  -pe: sharedmem use shared memory parallel environment, request CPUs
#  -o: redirect output with format jobname-jobID-hostname (jobname -N)
#  -e: redirect error with same format as output
# Initialise the environment modules
. /etc/profile.d/modules.sh

#!/usr/bin/env bash
export R_LIBS=%%R_LIBS%%
module load openmpi
module load igmm/apps/BEDTools 
module load igmm/apps/bowtie
module load igmm/apps/hdf5
module load igmm/apps/HISAT2
module load igmm/apps/pigz
module load igmm/apps/R/3.6.3
module load anaconda
source activate riboviz

DATAFOLDER="%%DATA_FOLDER%%"

echo "Running riboviz on dataset: ${DATAFOLDER}"

# move to scratch space
cd /exports/eddie/scratch/$USER

# make folder there
mkdir $DATAFOLDER
mkdir $DATAFOLDER/input

cd $DATAFOLDER/input

echo "${PWD}"

## get the dataset read files
module load igmm/apps/sratoolkit/2.10.8

##prefetch with Aspera client
%%PREFETCH%%

# use pigz to zip .fastq files into .fastq.gz files:
pigz *.fastq

echo "hopefully downloaded and pigz'd the files into /exports/eddie/scratch/$USER/${DATAFOLDER}/input"

# presumes I've already downloaded the SRA files & they're in /exports/eddie/scratch/$USER/$DATAFOLDER/input

# move to riboviz folder: 
cd $HOME/riboviz/riboviz

echo "moved to $HOME/riboviz/riboviz"

# make system link at riboviz folder to folder on scratch
ln -s /exports/eddie/scratch/$USER/$DATAFOLDER

# copy yaml into the riboviz/$DATAFOLDER folder, rename it
cp %%CONFIG_PATH%%/%%CONFIG_FILE%% $DATAFOLDER/

# presuming I'm in correct branch on riboviz

# move back up to riboviz folder (or nextflow won't run)
cd $HOME/riboviz/riboviz
echo "now in folder: ${PWD} ready to run"

# presuming .yaml config exists in $HOME/riboviz/riboviz/${DATAFOLDER} AND it points to ${DATAFOLDER}/input for files as required

# run nextflow validation: 
nextflow run prep_riboviz.nf -params-file ${DATAFOLDER}/%%CONFIG_FILE%% -work-dir ${DATAFOLDER}/work -ansi-log false --validate_only

# run nextflow: 
nextflow run prep_riboviz.nf -params-file ${DATAFOLDER}/%%CONFIG_FILE%% -work-dir ${DATAFOLDER}/work -ansi-log false

# hopefully success. 
echo "nextflow riboviz ${DATAFOLDER} data run complete"
