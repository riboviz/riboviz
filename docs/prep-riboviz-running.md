# Running the RiboViz workflow

This page describes how to run `prep_riboviz.py`. It is assumed you are familiar with [Configuring the RiboViz workflow](./prep-riboviz-config.md).

Contents:

* [Prepare input data](#prepare-input-data)
* [Set up your environment](#set-up-your-environment)
* [Configure number of processes (optional)](#configure-number-of-processes-optional)
* [Dry run `prep_riboviz.py`](#dry-run-prep_ribovizpy)
* [Run `prep_riboviz.py`](#run-prep_ribovizpy)
  - [Troubleshooting: `samtools sort: couldn't allocate memory for bam_mem`](#troubleshooting-samtools-sort-couldnt-allocate-memory-for-bam_mem)
  - [Troubleshooting: `WARNING: dedup_umis was TRUE but extract_umis was FALSE`](#troubleshooting-warning-dedup-umis-was-true-but-extract_umis-was-false)
  - [Troubleshooting: `Configuration parameter error: No sample files or multiplexed files are specified`](#troubleshooting-configuration-parameter-error-no-sample-files-or-multiplexed-files-are-specified)
  - [Troubleshooting: `Configuration parameter error: Both sample files (fq_files) and multiplexed files (multiplex_fq_files) are specified](#troubleshooting-configuration-parameter-error:-both-sample-files-and-multiplexed-files-are-specified)
  - [Troubleshooting: `Configuration parameter error: Multiplexed files ... are specified but no sample sheet`](#troubleshooting-configuration-parameter-error-multiplexed-files-are-specified-but-no-sample-sheet`)
* [Watch your disk space](#watch-your-disk-space)
* [Using-pre-generated hisat2 indices](#using-pre-generated-hisat2-indices)
* [Capturing commands submitted to bash](#capturing-commands-submitted-to-bash)
  - [Capturing bash commands and demultiplexing](#capturing-bash-commands-and-demultiplexing)
* [Customising logging](#customising-logging)
* [Exit codes](#exit-codes)

---

## Prepare input data

For each condition/sample in your experiment, merge all the FASTQ files with the ribosome profiling data for each condition/sample into a single FASTQ file. It is recommended that the files are also gzip-compressed to save space.

Multiple FASTQ files can be merged and compressed as follows:

```console
$ cat <SAMPLE_1_A>.fastq <SAMPLE_1_B>.fastq ... | gzip >> <SAMPLE_1>.fastq.gz
$ cat <SAMPLE_2_Y>.fastq <SAMPLE_2_Z>.fastq ... | gzip >> <SAMPLE_2>.fastq.gz
...etc...
```

Multiple gzipped FASTQ files can be merged as follows:

```console
$ cat <SAMPLE_1_A>.fastq.gz <SAMPLE_1_B>.fastq.gz ... | gzip >> <SAMPLE_1>.fastq.gz
```

Input files files should then be placed into a single input directory (`dir_in`). If this is not possible, then create an input directory and provide symbolic links to the actual files. For example, given:

```
samples/
  sample1/
    sample1.fastq.gz
  sample2/
    sample2.fastq.gz
```

These can be symlinked into a single directory, `data`, a sibling of `samples`, as follows:

```console
$ mkdir data
$ cd data/
$ ln -s ../samples/sample1/sample1.fastq.gz
$ ln -s ../samples/sample2/sample2.fastq.gz
$ ls -l
```
```
total 0
lrwxrwxrwx 1 user user 41 Aug 22 02:55 sample1.fastq.gz -> ../samples/sample1/sample1.fastq.gz
lrwxrwxrwx 1 user user 41 Aug 22 02:55 sample2.fastq.gz -> ../samples/sample2/sample2.fastq.gz
```

---

## Set up your environment

If you have not already done so, activate your Python environment:

* Miniconda Python 3.6+:

```console
$ source $HOME/miniconda3/bin/activate
```

If you have not already done so, set the paths to Hisat2 and Bowtie:

* If you followed [Create `setenv.sh` to configure Hisat2 and Bowtie paths](./install.md#create-setenvsh-to-configure-hisat2-and-bowtie-paths), then run:

```console
$ source $HOME/setenv.sh
```

* Otherwise, run the following (your directory names may be different, depending on the versions of Hisat2 and Bowtie you have):

```console
$ export PATH=~/hisat2-2.1.0:$PATH
$ export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH
```

---

## Configure number of processes (optional)

`prep_riboviz.py` can instruct the tools it invokes to use a specific number of processes. By default this is 1.

To configure `prep_riboviz.py` to use additional processes, in the configuration file, change:

```yaml
nprocesses: 1
```

* to the desired number of processes:

```yaml
nprocesses: 4
```

This parameter is currently used by `hisat2`, `samtools sort`, `bam_to_h5.R` and `generate_stats_figs.R`.

**Note:** for `cutadapt` the number of available processors on the host will be used regardless.

---

## Dry run `prep_riboviz.py`

`prep_riboviz.py` supports a `--dry-run` command-line parameter which can be used to validate the configuration. 

**Tip:** we strongly recommend doing a dry run before doing a live run on data you have not processed before.

Run `prep_riboviz.py` with `--dry-run` enabled:

* Either:

```console
$ python -m riboviz.tools.prep_riboviz --dry-run rscripts/ <CONFIG>
```

* Or:

```console
$ PYTHONPATH=. python riboviz/tools/prep_riboviz.py --dry-run rscripts/ <CONFIG>
```

where:

* `riboviz/tools/prep_riboviz.py`: path to `prep_riboviz.py`, relative to the RiboViz home directory.
* `--dry-run`: flag to enable the dry run.
* `rscripts/`: path to the directory with RiboViz's R scripts, relative to the RiboViz home directory.
* `<CONFIG>`: path to a YAML configuration file.

For more on this feature, see [Capturing commands submitted to bash](#capturing-commands-submitted-to-bash) below.

---

## Run `prep_riboviz.py`

Run `prep_riboviz.py`:

* Either:

```console
$ python -m riboviz.tools.prep_riboviz rscripts/ <CONFIG>
```

* Or:

```console
$ PYTHONPATH=. python riboviz/tools/prep_riboviz.py rscripts/ <CONFIG>
```

where:

* `riboviz/tools/prep_riboviz.py`: path to `prep_riboviz.py`, relative to the RiboViz home directory.
* `rscripts/`: path to the directory with RiboViz's R scripts, relative to the RiboViz home directory.
* `<CONFIG>`: path to a YAML configuration file.

Information on the key steps during processing is displayed. More detailed information, including the causes of any errors, is also added to a timestamped log file, in the current directory, named `riboviz-YYYYMMDD-HHMMSS.log` (for example, `riboviz.20190926-002455.log`).

To check the exit code, run:

```console
$ echo $0
```

See [Exit codes](#exit-codes), below, for a complete list of exit codes.

### Troubleshooting: `samtools sort: couldn't allocate memory for bam_mem`

If using more than one process (`nprocesses` > 1) you might get the error:

```
samtools sort: couldn't allocate memory for bam_mem
```

This leads to a failure to create `.bam` and `.bam.bai` files.

In this case, you may need to explicitly set the amount of memory per thread in calls to `samtools sort`.

Check how much free memory you have e.g.

```console
$ free --mega
              total        used        free      shared  buff/cache   available
Mem:           2017         684        1028           2         303        1181
Swap:           969         619         350
```

Divide the free memory by the number of processes, `nprocesses` e.g. 1024/4 = 256 MB.

Edit `riboviz/workflow.py` and change the lines:

```python
cmd_sort = ["samtools", "sort", "-@", str(run_config.nprocesses),
            "-O", "bam", "-o", bam_file, "-"]
```

to include the `samtools` flag `-m <MEMORY_DIV_PROCESSES>M` e.g.:

```python
cmd_sort = ["samtools", "sort", "-@", str(run_config.nprocesses),
            "-m", "256M",
            "-O", "bam", "-o", bam_file, "-"]
```

### Troubleshooting: `WARNING: dedup_umis was TRUE but extract_umis was FALSE`

This warning means that in your YAML configuration file you have defined:

```yaml
extract_umis: FALSE
dedup_umis: TRUE
```

While `prep_riboviz.py` gives a warning, it will continue execution.

Unless you explicitly want this to occur you should:

* Either, set `extract_umis` to `TRUE`, if you want UMI extraction to occur.
* Or, set `dedup_umis` to `FALSE`, if you do not want deduplication to occur.

### Troubleshooting: `Configuration parameter error: No sample files or multiplexed files are specified`

If you see:

```
Configuration parameter error: No sample files (fq_files) or
multiplexed files (multiplex_fq_files) are specified
```

then this means you provided neither `fq_files` nor `multiplex_fq_files` parameters in your configuration. Add the missing configuration parameter.

### Troubleshooting: `Configuration parameter error: Both sample files and multiplexed files are specified`

If you see:

```
Configuration parameter error: Both sample files (fq_files) and
multiplexed files (multiplex_fq_files) are specified 
```

then this means you provided both `fq_files` and `multiplex_fq_files` parameters in your configuration. Only a group of non-multiplexed files or a single multiplexed file can be provided. Delete one of these configuration parameters.

### Troubleshooting: `Configuration parameter error: Multiplexed files are specified but no sample sheet`

If you see:

```
Configuration parameter error: Multiplexed files (multiplex_fq_files)
are specified but no sample sheet (sample_sheet)
```

then this means you provided a `multiplex_fq_files` parameter in your configuration but not a complementary `sample_sheet` parameter. Add the missing configuration parameter.

---

## Watch your disk space

`prep_riboviz.py` generates many intermediate files and some of these may be unompressed and **large**, i.e. about the same size as the input files. All these files are placed in a temporary directory (`dir_tmp`). The temporary directory's contents can be inspected for troubleshooting, if necessary.

Similarly, `prep_riboviz.py` creates numerous log files.

For example, here is the volume of the outputs from a run of the vignette as documented in [Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file](./run-vignette.md):

| Directory |   MB |
| --------- | ---- |
| `index`   |    9 |
| `tmp`     | 1040 |
| `output`  |    3 |
| `logs`    |    1 |
| Total     | 1053 |

**Tip:** We recommend you delete temporary directories and log directories when you have completed your analyses to your satisfaction.

---

## Using pre-generated hisat2 indices

If you have already generated hisat2 indices for the same organism and annotation, set `build_indices` to `FALSE` and specify the index directory (`dir_index`) which contains the indices.

---

## Capturing commands submitted to bash

When run in `--dry-run` mode or live, `prep_riboviz.py` logs the commands it submits to bash to execute local scripts and third-party tools. These are logged into a bash script (`cmd_file`). This script can be used for debugging. However it can also be run as-is and the effect is the same as running, or rerunning, `prep_riboviz.py`.

For example, if `cmd_file` had value `run_riboviz.sh` then the script can be run as follows:

```console
$ bash run_riboviz.sh
```

### Capturing bash commands and demultiplexing 

When run with demultiplexing enabled, `prep_riboviz.py` determines the demultiplexed files to process based on the samples in the `num_reads.tsv` file, output by `demultipled_fastq.py`, which have one or more reads.

If running `prep_riboviz.py` in dry run mode then this file does not exist. In this case, `prep_riboviz.py` will determine the samples based on the samples in the sample sheet.

This means that the bash script output by `prep_riboviz.py` with `--dry-run` may not match that output during a live run, if, during the live run, some samples had 0 associated reads.

---

## Customising logging

You can customise logging by editing the file `riboviz/logging.yaml`.

If you do not want `riboviz.log` to include a timestamp (i.e. you want the log file to be named `riboviz.log`) then edit `riboviz/logging.yaml` and replace:

```yaml
  handlers: [console, timestamp_file_handler]
```

with:

```yaml  
  handlers: [console, file_handler]
```

---

## Exit codes

`prep_riboviz.py` returns the following exit codes:

* 0: Processing successfully completed.
* 1: A file does not seem to exist.
* 2: Errors occurred loading or accessing configuration e.g. missing configuration parameters, inconsistent configuration parameters.
* 3: Error occurred during processing.
