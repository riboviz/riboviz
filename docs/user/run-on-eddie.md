# Running the RiboViz workflow on Eddie

This page describes how you can run **RiboViz** on [Eddie](https://www.ed.ac.uk/information-services/research-support/research-computing/ecdf/high-performance-computing), The University of Edinburgh ECDF Linux Compute Cluster.

**Note:** This information is for University of Edinburgh users only.

The Eddie service documentation is on the University of Edinburgh [wiki](https://www.wiki.ed.ac.uk/display/ResearchServices/Eddie),

These guidelines may be useful for running **RiboViz** in other HPC systems.

All Python and R packages required to run **RiboViz** have been installed in `/exports/csce/eddie/biology/groups/wallace_rna` on Eddie

Contents:

* [Logging in](#logging-in)
* [Configure Anaconda enviroment](#configure-anaconda-enviroment)
* [Get Riboviz and example-datasets](#get-riboviz-and-example-datasets)
* [Interactive sessions](#interactive-sessions)
* [Set up your environment from scratch (optional)](#set-up-your-environment-from-scratch-optional)
* [Create `set-riboviz-env.sh`](#create-set-riboviz-envsh)
* [Run a "vignette" of the RiboViz workflow in an interactive node](#run-a-vignette-of-the-riboviz-workflow)
* [Job submission](#job-submission)
  - [Submitting jobs](#submitting-jobs)
  - [Monitoring jobs](#monitoring-jobs)
  - [Cancelling jobs](#cancelling-jobs)
  - [Job accounting](#job-accounting)
 * [Run a full-size example dataset](#run-a-full-size-example-dataset)
   - [Create directories for input paths](#create-directories-for-input-paths)
   - [Download SRR files](#download-srr-files)
   - [Create `qsub` script](#create-qsub-script)

---

## Logging in

Connect to the cluster using `ssh` from a terminal window (Linux and Mac OS) or use a client such as MobaXterm (Windows):

```console
$ ssh -X <YOUR_UUN>@eddie.ecdf.ed.ac.uk
```

In the rest of this document, we shall abbreviate your universal username as `$USER`.

**Note** that access to the cluster is only available from the University network. External users should first connect to the University network using the [VPN service](https://www.ed.ac.uk/information-services/computing/desktop-personal/vpn).

---

## Configure Anaconda enviroment

Configure your `.condarc file` to point to the `anaconda directory` in `/exports/csce/eddie/biology/groups/wallace_rna/`.

If you do not have a `.condarc` file in your home directory, create it first.

Add:

```
envs_dirs:
  - /exports/csce/eddie/biology/groups/wallace_rna/anaconda/envs
pkgs_dirs:
  - /exports/csce/eddie/biology/groups/wallace_rna/anaconda/pkgs
```

---

## Get **RiboViz** and example-datasets

Get **RiboViz** and example-datasets:

```console
$ mkdir riboviz
$ cd riboviz
$ git clone https://github.com/riboviz/riboviz
$ git clone https://github.com/riboviz/example-datasets
```

**Note:** Your home directory space is enough for running a vignette but is not enough for running a full-size dataset.

We recommend using the cluster filesystem (`/exports/[COLLEGE]/eddie/...`) for storing **RiboViz** and `example-datasets`.

If you do not have a group space, you can use your scratch directory (`/exports/eddie/scratch/<YOUR_UUN>`)

---

## Interactive sessions

There are a limited number of nodes that accept interactive login sessions, to allow you to run interactive jobs or graphical  applications. To start an interactive session run:

```console
$ qlogin -l h_vmem=16G
```

`-l h_vmem` means that you ask for 16GB RAM

If you have access to a priority queue then you can use:

```console
$ qlogin -P <QUEUE_NAME> -l h_vmem=16G
```

See [Interactive sessions](https://www.wiki.ed.ac.uk/display/ResearchServices/Interactive+Sessions) for more information.

**Troubleshooting: fail to enter interactive node**

If you see:

```
Your job 2674903 ("QLOGIN") has been submitted
waiting for interactive job to be scheduled ...timeout (5 s) expired while waiting on socket fd 9

Your "qlogin" request could not be scheduled, try again later.
```

There may be no free nodes at present. Alternatively, Eddie may be under maintenance. You can check Eddie's status on the [Information Systems Alerts](https://alerts.is.ed.ac.uk/).

Either way, you have to wait for a free node to become available or for Eddie to come back up. It usually won't take too long.

---

## Set up your environment from scratch (optional)

### Activate Anaconda

`/exports/csce/eddie/biology/groups/wallace_rna` has Anaconda packages (in `anaconda`) and all the Python packages required by **RiboViz** are there, accessible as a `riboviz` conda environment.

```console
$ source activate riboviz
````

### Configure R packages path

```console
$ export R_LIBS=/exports/csce/eddie/biology/groups/wallace_rna/Rlibrary
```

### Load necessary modules on node

```console
$ module load igmm/apps/BEDTools
$ module load igmm/apps/bowtie
$ module load igmm/apps/hdf5
$ module load igmm/apps/HISAT2
$ module load igmm/apps/pigz
$ module load igmm/apps/R/3.6.3
$ module load anaconda
```

---

## Create `set-riboviz-env.sh`

You can create a script named `set-riboviz-env.sh` for above commands to set up your environment:

```
#!/usr/bin/env bash
export R_LIBS=/exports/csce/eddie/biology/groups/wallace_rna/Rlibrary
module load igmm/apps/BEDTools
module load igmm/apps/bowtie
module load igmm/apps/hdf5
module load igmm/apps/HISAT2
module load igmm/apps/pigz
module load igmm/apps/R/3.6.3
module load anaconda
source activate riboviz
```

In future you need only to run:

```console
$ source set-riboviz-env.sh
````

---

## Run a "vignette" of the **RiboViz** workflow

Change into the **RiboViz** repository:

```console
$ cd riboviz/riboviz
```

To run the Python workflow:

```console
$ python -m riboviz.tools.prep_riboviz -c vignette/vignette_config.yaml
Running under Python: 3.7.6 | packaged by conda-forge | (default, Jun  1 2020, 18:57:50)
[GCC 7.5.0]
Created by: RiboViz Date: 2020-06-06 02:09:31.484844 Command-line tool: /exports/csce/eddie/biology/groups/wallace_rna/riboviz/riboviz/tools/prep_riboviz.py File: /exports/csce/eddie/biology/groups/wallace_rna/riboviz/riboviz/tools/prep_riboviz.py Version: commit 0fe0191f585d55763d00283d129a12e4c3c1e5c1 date 2020-06-04 00:43:45-07:00
Configuration file: vignette/vignette_config.yaml
Command file: run_riboviz_vignette.sh
Number of processes: 1
Build indices for alignment, if necessary/requested
Build indices for alignment (vignette/input/yeast_rRNA_R64-1-1.fa). Log: vignette/logs/20200606-020931/hisat2_build_r_rna.log
Build indices for alignment (vignette/input/yeast_YAL_CDS_w_250utrs.fa). Log: vignette/logs/20200606-020931/hisat2_build_orf.log
Processing samples
Processing sample: WTnone
Processing file: vignette/input/SRR1042855_s1mi.fastq.gz
Cut out sequencing library adapters. Log: vignette/logs/20200606-020931/WTnone/01_cutadapt.log
Remove rRNA or other contaminating reads by alignment to rRNA index files. Log: vignette/logs/20200606-020931/WTnone/02_hisat2_rrna.log
Align remaining reads to ORFs index files using hisat2. Log: vignette/logs/20200606-020931/WTnone/03_hisat2_orf.log
Trim 5' mismatches from reads and remove reads with more than 2 mismatches. Log: vignette/logs/20200606-020931/WTnone/04_trim_5p_mismatch.log
Convert SAM to BAM and sort on genome. Log: vignette/logs/20200606-020931/WTnone/05_samtools_view_sort.log
Index BAM file. Log: vignette/logs/20200606-020931/WTnone/06_samtools_index.log
Calculate transcriptome coverage for + strand and save as a bedgraph. Log: vignette/logs/20200606-020931/WTnone/07_bedtools_genome_cov_plus.log
Calculate transcriptome coverage for - strand and save as a bedgraph. Log: vignette/logs/20200606-020931/WTnone/08_bedtools_genome_cov_minus.log
Make length-sensitive alignments in H5 format. Log: vignette/logs/20200606-020931/WTnone/09_bam_to_h5.log
Create summary statistics, and analyses and QC plots for both RPF and mRNA datasets. Log: vignette/logs/20200606-020931/WTnone/10_generate_stats_figs.log
Finished processing sample: vignette/input/SRR1042855_s1mi.fastq.gz
Processing sample: WT3AT
Processing file: vignette/input/SRR1042864_s1mi.fastq.gz
Cut out sequencing library adapters. Log: vignette/logs/20200606-020931/WT3AT/01_cutadapt.log
Remove rRNA or other contaminating reads by alignment to rRNA index files. Log: vignette/logs/20200606-020931/WT3AT/02_hisat2_rrna.log
Align remaining reads to ORFs index files using hisat2. Log: vignette/logs/20200606-020931/WT3AT/03_hisat2_orf.log
Trim 5' mismatches from reads and remove reads with more than 2 mismatches. Log: vignette/logs/20200606-020931/WT3AT/04_trim_5p_mismatch.log
Convert SAM to BAM and sort on genome. Log: vignette/logs/20200606-020931/WT3AT/05_samtools_view_sort.log
Index BAM file. Log: vignette/logs/20200606-020931/WT3AT/06_samtools_index.log
Calculate transcriptome coverage for + strand and save as a bedgraph. Log: vignette/logs/20200606-020931/WT3AT/07_bedtools_genome_cov_plus.log
Calculate transcriptome coverage for - strand and save as a bedgraph. Log: vignette/logs/20200606-020931/WT3AT/08_bedtools_genome_cov_minus.log
Make length-sensitive alignments in H5 format. Log: vignette/logs/20200606-020931/WT3AT/09_bam_to_h5.log
Create summary statistics, and analyses and QC plots for both RPF and mRNA datasets. Log: vignette/logs/20200606-020931/WT3AT/10_generate_stats_figs.log
Finished processing sample: vignette/input/SRR1042864_s1mi.fastq.gz
File not found: vignette/input/example_missing_file.fastq.gz
Finished processing 3 samples, 1 failed
Collate TPMs across sample results. Log: vignette/logs/20200606-020931/collate_tpms.log
Count reads. Log: vignette/logs/20200606-020931/count_reads.log
Completed
```

To run the Nextflow workflow:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.04.1
Launching `prep_riboviz.nf` [big_shirley] - revision: 281e0a4d55
No such sample file (NotHere): example_missing_file.fastq.gz
[eb/da58ed] Submitted process > cutAdapters (WT3AT)
[40/456217] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[c7/83b32a] Submitted process > cutAdapters (WTnone)
[1c/9c76ba] Submitted process > buildIndicesrRNA (yeast_rRNA)
[75/c7ac3b] Submitted process > hisat2rRNA (WT3AT)
[8a/52ead5] Submitted process > hisat2rRNA (WTnone)
[7c/be343a] Submitted process > hisat2ORF (WT3AT)
[57/d2ac14] Submitted process > hisat2ORF (WTnone)
[db/6b0cd9] Submitted process > trim5pMismatches (WT3AT)
[0e/f163f9] Submitted process > trim5pMismatches (WTnone)
[db/80a57c] Submitted process > samViewSort (WT3AT)
[6b/4e6432] Submitted process > samViewSort (WTnone)
[17/0557d6] Submitted process > outputBams (WT3AT)
[5e/236c46] Submitted process > outputBams (WTnone)
[b0/be630a] Submitted process > makeBedgraphs (WT3AT)
[4a/7dc7de] Submitted process > bamToH5 (WT3AT)
[50/cb98ce] Submitted process > makeBedgraphs (WTnone)
[0d/ab3d58] Submitted process > bamToH5 (WTnone)
[f1/7ba436] Submitted process > generateStatsFigs (WT3AT)
[af/08787f] Submitted process > generateStatsFigs (WTnone)
Finished processing sample: WT3AT
[68/53d480] Submitted process > renameTpms (WT3AT)
Finished processing sample: WTnone
[4d/815c11] Submitted process > renameTpms (WTnone)
[77/8cebac] Submitted process > collateTpms (WT3AT, WTnone)
[d8/a8fcb5] Submitted process > countReads
Workflow finished! (OK)
```

For more information about the vignette, see [Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file](./run-vignette.md).

---

## Job submission

Computational work on Eddie is usually submitted to the cluster as batch jobs initiated from a login node. In order to submit a job you need to write a Grid Engine job submission script containing details of the program to run as well as requests for resources. Then, you submit this job script to the cluster with the `qsub` command.

You can create a job script named `job_riboviz.sh` in your `riboviz` directory to run a **RiboViz** workflow:

```
#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N riboviz_vignette              
#$ -cwd                  
#$ -l h_rt=01:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 16
#$ -o $JOB_NAME-$JOB_ID-$HOSTNAME.o
#$ -e $JOB_NAME-$JOB_ID-$HOSTNAME.e
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 1 hour: -l h_rt
#  ask for 8 Gbyte RAM: -l h_vmem
#  use shared memory parallel environment, request 16 CPUs
#  redirect output with format jobname-jobID-hostname (jobname -N)
#  redirect error with same format as output

# Initialise the environment modules
. /etc/profile.d/modules.sh

export R_LIBS=/exports/csce/eddie/biology/groups/wallace_rna/Rlibrary
module load igmm/apps/BEDTools
module load igmm/apps/bowtie
module load igmm/apps/hdf5
module load igmm/apps/HISAT2
module load igmm/apps/pigz
module load igmm/apps/R/3.6.3
module load anaconda
source activate riboviz

# Uncomment this to run the python workflow:
#python -m riboviz.tools.prep_riboviz -c vignette/vignette_config.yaml

# Run the Nextflow workflow:
nextflow run prep_riboviz.nf -params-file vignette/vignette_config.yaml -ansi-log false
```

### Submitting jobs

Change into the **RiboViz** repository:

```console
$ cd riboviz/riboviz
```

Run:

```console
$ qsub job_riboviz.sh
```

If you have access to a priority queue then you can use:

```console
$ qsub -P <QUEUE_NAME> job_riboviz.sh
```

A job ID will be displayed.

This will output the standard output from `prep_riboviz.py` or `prep_riboviz.nf` (depending on which option you are running) to a file, `riboviz_vignette-$JOB_ID-$HOSTNAME.o`, in the current working directory, and errors to a file, `riboviz_vignette-$JOB_ID-$HOSTNAME.e`.

The contents of `riboviz_vignette-$JOB_ID-$HOSTNAME.o` should be the same as the standard output of [Run a "vignette" of the RiboViz workflow in an interactive node](#run-a-vignette-of-the-RiboViz-workflow) above.

### Monitoring jobs

Active jobs (i.e. pending or running) can be monitored with the `qstat` command

```console
$ qstat
job-ID     prior   name       user         state submit/start at     queue                          jclass                                                     slots ja-task-ID
----------------------------------------------------------------------------------------------------------------------                            --------------------------
   2701173 0.00000 riboviz_vi $USER     qw    06/11/2020 13:22:28                                                                                               1
```

### Cancelling jobs

If you want to kill a job you've submitted to Eddie, use the `qdel` command

```console
$ qdel [Your-Job-ID]
```

or

```console
$ qdel [Your-Job-Name]
```

### Job accounting

Detailed information about completed jobs is available with the `qacct` command:

`$ qacct -j [Your-Job-ID]`

```console
$ qacct -j 2701137
==============================================================
qname        eddie
hostname     node1c17.ecdf.ed.ac.uk
group        eddie_users
owner        $USER
project      uoe_baseline
department   defaultdepartment
jobname      riboviz_vignette
jobnumber    2701173
taskid       undefined
pe_taskid    NONE
account      sge
priority     0
cwd          /home/$USER/riboviz/riboviz
submit_host  login02.ecdf.ed.ac.uk
submit_cmd   qsub /exports/eddie3_homes_local/$USER/job_riboviz.sh
qsub_time    06/11/2020 13:22:28.652
start_time   06/11/2020 13:22:37.680
end_time     06/11/2020 13:27:51.952
granted_pe   NONE
slots        1
failed       0
deleted_by   NONE
exit_status  0
ru_wallclock 314.272
ru_utime     290.658
ru_stime     10.078
ru_maxrss    904316
ru_ixrss     0
ru_ismrss    0
ru_idrss     0
ru_isrss     0
ru_minflt    3064424
ru_majflt    0
ru_nswap     0
ru_inblock   0
ru_oublock   1696
ru_msgsnd    0
ru_msgrcv    0
ru_nsignals  0
ru_nvcsw     198940
ru_nivcsw    117557
wallclock    314.333
cpu          300.736
mem          179.488
io           7.731
iow          0.000
ioops        1566484
maxvmem      1.656G
maxrss       0.000
maxpss       0.000
arid         undefined
jc_name      NONE
bound_cores  0,4

```

See [Job submission](https://www.wiki.ed.ac.uk/display/ResearchServices/Job+Submission) for more information.

---

## Run a Full-size Example Dataset

In this example, we're using the Wallace et al. 2020 *Cryptococcus neoformans* 'JEC21' dataset from the [Example-Datasets repository](https://github.com/riboviz/example-datasets). This example dataset repository contains .yaml config files, annotation files and contaminant files for a range of different publically available datasets across a range of organisms.

To run the `Wallace_2020_JEC21` dataset on Eddie, logout from any interactive node you may be logged into (for example, if you were running the vignette example above) and ensure you are within the `example-datasets` repository at $HOME/riboviz/example-datasets and that you are in the correct git branch for both riboviz and example-datasets repositories.

NOTE: the following sections are here for information as it might be helpful in explaining to new users what these steps do, and how to adjust these steps for your own data.  These steps (except [initial setup of SRA Toolkit](### Download fastq data files from the Short Read Archive (SRA): initial setup)) are included in the sample job submission script for the Wallace_2020_JEC21 dataset, so you don't have to carry out these steps manually if you plan to run the [job submission script](#Create-qsub-script) and you can skip straight there if you want to try running the script.

Please also note the paths in the YAML configuration file we will be using from the `example-datasets` directory are just a reference. You should check and edit the paths according to your directory structure.

```
cd example-datasets
git checkout master
cd ..   # back to $HOME/riboviz
cd riboviz
git checkout develop

## check you have this file structure:
## $HOME/riboviz/riboviz  # (branch: develop)
## $HOME/riboviz/example-datasets  # (branch: master)
```
This will give you access to the correct config.yaml, annotation and contaminants files.

### Create directories for input paths

Create a directory named `Wallace_2020_JEC21` in `/exports/eddie/scratch/$USER/` and a directory within that called `input`
```
$ mkdir Wallace_2020_JEC21
$ mkdir Wallace_2020_JEC21/input
```
### Set up riboviz dataset folder and create system links to scratch folders

Move to the main riboviz folder
```
$ cd $HOME/riboviz/riboviz
```
Create a system link between a new folder for our dataset and the folder on scratch which will hold our inputs and outputs:
```
# make system link at riboviz folder to folder on scratch
$ ln -s /exports/eddie/scratch/$USER/$Wallace_2020_JEC21
```
Now we copy the yaml across from the example-datasets folder, into our main riboviz folder. If you wish to edit the yaml, then it's best to edit this version, in $HOME/riboviz/riboviz/Wallace_2020_JEC21.  

```console
$ mkdir annotation
$ cp /exports/eddie/scratch/s1919303/riboviz/example-datasets/fungi/cryptococcus/annotation/JEC21_10p_up12dwn9_CDS_with_120bputrs.fa annotation

# copy yaml into the riboviz/Wallace_2020_JEC21 folder, rename it
 # cp [example-datasets version yaml] [our 'local' riboviz folder version]

$ cp $HOME/riboviz/example-datasets/fungi/cryptococcus/Wallace_2020_JEC21_2-samples_10p_up12dwn9_CDS_120bpL_120bpR_config.yaml Wallace_2020_JEC21/Wallace_2020_JEC21_2-samples_10p_up12dwn9_CDS_120bpL_120bpR_config.yaml

cd $HOME/riboviz/riboviz
```
We need to make sure we move back into the main riboviz folder, where we will be ready to run the nextflow commands.

### Download fastq data files from the Short Read Archive (SRA): initial setup

Eddie allows us to load the [SRA Toolkit](https://github.com/ncbi/sra-tools) module, including the utility `fasterq-dump` for downloading data files.  This utility has been included in SRA Toolkit since version 2.9.1. We recommend using `fasterq-dump`.

An earlier tool, `fastq-dump`, is also included in SRA Toolkit, however, you may find it is too slow for `fastq-dump` to download a large dataset like `Wallace_2020_JEC21` which is around 50GB uncompressed. Even using the `--gzip` option to directly download the `.gz` file may be too slow.

A faster alternative can be to use `fasterq-dump` and the Aspera client's `prefetch` tool (which is provided as part of the above module), as recommended here.

To get `fasterq-dump`, follow SRA Toolkit's [02. Installing SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) to install the latest version of the SRA Toolkit (follow the instructions for CentOS).

Before your initial use of SRA toolkit, configure download and cache settings, by running:

```
$ module load igmm/apps/sratoolkit/2.10.8
$ vdb-config --interactive
```
then follow the interactive prompts (using tab to navigate through the menus) and edit the `CACHE` > `Set Default Import Path` section to change the workspace location.

This path adjusts where the tool puts your cache directory, which could get very large (100s of GB). We recommend using your scratch space `/exports/eddie/scratch/$USER/ncbi`, where `$USER` is replaced by your username.

You may have to repeat the `vdb-config` step periodically, as data on Eddie scratch space is automatically cleared after one month.
For more information about the configuration utility, see [SRA toolkit installation and configuration Guide](https://ncbi.github.io/sra-tools/install_config.html).

### Download fastq data files from the Short Read Archive (SRA): subsequent uses

To get the Wallace_2020_JEC21 dataset .fastq read files (remember to change to the `Wallace_2020_JEC21/input` directory on scratch first):

```
cd /exports/eddie/scratch/$USER/Wallace_2020_JEC21/input
$ module load igmm/apps/sratoolkit/2.10.8

# prefetch with Aspera client
$ prefetch SRR9620588 SRR9620586

$ fasterq-dump SRR9620588
$ fasterq-dump SRR9620586
```

These download utilities do not have an option to compress (gzip) the files, nor apparently allow you to pipe their output into another program. So we use the `pigz` utility to compress.

```
$ module load igmm/apps/pigz
$ pigz *.fastq
```

It may be helpful to test this functionality with a smaller download file, e.g. [SRR014376](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR014376) from Ingolia 2009.

Alternatively, it is possible to download `.fastq.gz` format files of SRA data from the European Nucleotide Archive, but we have not tested the speed. For example, from Ingolia 2009 again, [ftp link to SRR014376.fastq.gz](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR014/SRR014376/SRR014376.fastq.gz).

### Create `qsub` script

Create the job submission script in `$HOME` or a location of your choosing, and name it something like `run_W-Cn-JEC21_2020.sh`.  

```
#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N W-Cn-JEC21_2020              
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

DATAFOLDER="Wallace_2020_JEC21"

echo "Running riboviz on dataset: ${DATAFOLDER}"

# move to scratch space
cd /exports/eddie/scratch/$USER

# make folder there
mkdir $DATAFOLDER
mkdir $DATAFOLDER/input

cd $DATAFOLDER/input

echo "${PWD}"

# get the dataset read files
module load igmm/apps/sratoolkit/2.10.8

##prefetch with Aspera client
prefetch SRR9620588 SRR9620586
fasterq-dump SRR9620588
fasterq-dump SRR9620586

# use pigz to zip .fastq files into .fastq.gz files:
pigz *.fastq

echo "hopefully downloaded and pigz'd the files into /exports/eddie/scratch/$USER/${DATAFOLDER}/input"

# presumes downloaded the SRA files OK & they're in /exports/eddie/scratch/$USER/$DATAFOLDER/input

# move to riboviz folder:
cd $HOME/riboviz/riboviz

echo "moved to $HOME/riboviz/riboviz"

# make system link at riboviz folder to folder on scratch
ln -s /exports/eddie/scratch/$USER/$DATAFOLDER

# copy yaml into the riboviz/$DATAFOLDER folder, rename it
cp $HOME/riboviz/example-datasets/fungi/cryptococcus/Wallace_2020_JEC21_2-samples_10p_up12dwn9_CDS_120bpL_120bpR_config.yaml $DATAFOLDER/Wallace_2020_JEC21_2-samples_10p_up12dwn9_CDS_120bpL_120bpR_config.yaml

# presuming I'm in correct branch on riboviz

# move back up to riboviz folder (or nextflow won't run)
cd $HOME/riboviz/riboviz
echo "now in folder: ${PWD} ready to run"

# presuming .yaml config exists in $HOME/riboviz/riboviz/${DATAFOLDER} AND it points to ${DATAFOLDER}/input for files as required

# run nextflow validation:
nextflow run prep_riboviz.nf -params-file ${DATAFOLDER}/Wallace_2020_JEC21_2-samples_10p_up12dwn9_CDS_120bpL_120bpR_config.yaml -work-dir ${DATAFOLDER}/work -ansi-log false --validate_only

# run nextflow:
nextflow run prep_riboviz.nf -params-file ${DATAFOLDER}/Wallace_2020_JEC21_2-samples_10p_up12dwn9_CDS_120bpL_120bpR_config.yaml -work-dir ${DATAFOLDER}/work -ansi-log false

# hopefully success.
echo "nextflow riboviz ${DATAFOLDER} data run complete"
```

### Job Submission

Check that you are in the same location as your submission script, or remember to add that path to your `qsub` command.

```console
$ cd /exports/eddie/scratch/s1919303/riboviz/riboviz/
$ qsub [Your-script]
```

If you run the `example-dataset` in your scratch space, remember to move the output data to DataStore or other persistent storage after the jobs have finished.

### Checking outputs

The job submission should create two files: an output file `JOB_NAME-$JOB_ID-$HOSTNAME.o` and an error file `$JOB_NAME-$JOB_ID-$HOSTNAME.e`.  These are the best place to start looking after a job has completed, to check if it has run successfully.

The output file will contain the standard output from the nextflow run, and will give you the locations within `/work` for the folders associated with each process. More information on this can be found in the user documentation for [Debugging and Bash Scripts](../docs/user/prep-riboviz-run-nextflow.md#debugging-and-bash-scripts)

The output files will be in `/exports/eddie/scratch/$USER/Wallace_2020_JEC21/output/`.

Another file worth checking if you are uncertain how a nextflow run performed is the `.nextflow.log` file found in `$HOME/riboviz/riboviz`.

### Moving and Downloading Outputs

If you run the example-dataset in your scratch space as detailed in these instructions, remember to move the output data to DataStore or other persistent storage after the jobs have finished.  

Alternatively, you can use `scp` and `scp -r` to transfer the files and folders via the commandline to a local machine (e.g. `scp -r UUN@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/$USER/Wallace_2020_JEC21/output/* $HOME/Downloads/` to download all of the output files from the scratch space to the user's local Downloads file).

You can check the file sizes using `du -ch` to get an idea of how much space you will need to have available before transferring the files.  The `/input` files for this example dataset total ~18GB, the `/output` files total ~8GB, but the total size of the `Wallace_2020_JEC21` folder on scratch after a successful run currently comes to 551GB!

Files older than one month are removed from this directory automatically.

See [Storage](https://www.wiki.ed.ac.uk/display/ResearchServices/Storage) for more information.
