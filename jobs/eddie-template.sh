#!/bin/sh
# Grid Engine options (lines prefixed with #$)
# Job name:
#$ -N %%job_name%%
# Use the current working directory:
#$ -cwd
# Runtime limit:
#$ -l h_rt=%%job_runtime%%
# RAM:
#$ -l h_vmem=%%job_memory%%
# Use shared memory parallel environment and request number of CPUs:
#$ -pe sharedmem %%job_num_cpus%%
# Redirected output file name format:
#$ -o $JOB_NAME-$JOB_ID-$HOSTNAME.o
# Redirected error file name format:
#$ -e $JOB_NAME-$JOB_ID-$HOSTNAME.e
# E-mail address for batch job events:
#$ -M %%job_email%%
# Events triggering emails about batch job.
# Any combination of b - begin, e - end, a - abort, s - suspend:
#$ -m %%job_email_events%%

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
module load igmm/apps/sratoolkit/2.10.8
module load anaconda
source activate riboviz

echo "Running Nextflow riboviz..."

nextflow run prep_riboviz.nf -params-file %%config_file%% -work-dir %%nextflow_work_dir%% -ansi-log false -with-report %%nextflow_report_file%% %%validate_only%%

echo "Nextflow riboviz run complete!"
