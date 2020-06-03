# Running the RiboViz workflow

This page describes how to run `riboviz.tools.prep_riboviz`. It is assumed you are familiar with [Configuring the RiboViz workflow](./prep-riboviz-config.md).

Contents:

* [Prepare input data](#prepare-input-data)
* [Set up your environment](#set-up-your-environment)
* [Configure number of processes (optional)](#configure-number-of-processes-optional)
* [Dry run `prep_riboviz`](#dry-run-prep_riboviz)
* [Run `prep_riboviz`](#run-prep_riboviz)
  - [Troubleshooting: `This script needs to be run under Python 3`](#troubleshooting-this-script-needs-to-be-run-under-python-3)
  - [Troubleshooting: `samtools sort: couldn't allocate memory for bam_mem`](#troubleshooting-samtools-sort-couldnt-allocate-memory-for-bam_mem)
  - [Troubleshooting: `WARNING: dedup_umis was TRUE but extract_umis was FALSE`](#troubleshooting-warning-dedup_umis-was-true-but-extract_umis-was-false)
  - [Troubleshooting: `Configuration parameter error: No sample files or multiplexed files are specified`](#troubleshooting-configuration-parameter-error-no-sample-files-or-multiplexed-files-are-specified)
  - [Troubleshooting: `Configuration parameter error: Both sample files and multiplexed files are specified`](#troubleshooting-configuration-parameter-error-both-sample-files-and-multiplexed-files-are-specified)
  - [Troubleshooting: `Configuration parameter error: Multiplexed files are specified but no sample sheet`](#troubleshooting-configuration-parameter-error-multiplexed-files-are-specified-but-no-sample-sheet)
* [Invoking `prep_riboviz` from outwith the RiboViz home directory](#invoking-prep_riboviz-from-outwith-the-riboviz-home-directory)
* [Managing your disk usage](#managing-your-disk-usage)
* [Using-pre-generated hisat2 indices](#using-pre-generated-hisat2-indices)
* [Capturing commands submitted to bash](#capturing-commands-submitted-to-bash)
  - [Capturing bash commands and demultiplexing](#capturing-bash-commands-and-demultiplexing)
* [Customising logging](#customising-logging)
  - [Removing timestamps](#removing-timestamps)
  - [Using custom log configuration files](#using-custom-log-configuration-files)
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

```console
$ source $HOME/miniconda3/bin/activate
$ conda activate riboviz
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

`prep_riboviz` can instruct the tools it invokes to use a specific number of processes. By default this is 1.

To configure `prep_riboviz` to use additional processes, in the configuration file, change:

```yaml
num_processes: 1
```

* to the desired number of processes:

```yaml
num_processes: 4
```

This parameter is currently used by `hisat2`, `samtools sort`, `bam_to_h5.R` and `generate_stats_figs.R`.

**Note:** for `cutadapt` the number of available processors on the host will be used regardless.

---

## Dry run `prep_riboviz`

`prep_riboviz` supports a `-d` (or `--dry-run`) command-line parameter which can "dry run" the workflow. This validates the configuration without executing the workflow steps.

**Tip:** we strongly recommend doing a dry run before doing a live run on data you have not processed before.

Run `prep_riboviz` in dry run mode:

```console
$ python -m riboviz.tools.prep_riboviz -d -c <CONFIG_FILE>
```

where:

* `-d`: flag to enable the dry run.
* `<CONFIG_FILE>`: path to a YAML configuration file.

For more on this feature, see [Capturing commands submitted to bash](#capturing-commands-submitted-to-bash) below.

---

## Run `prep_riboviz`

Run `prep_riboviz`:

```console
$ python -m riboviz.tools.prep_riboviz -c <CONFIG_FILE>
```

where:

* `<CONFIG_FILE>`: path to a YAML configuration file.

Information on the key steps during processing is displayed. More detailed information, including the causes of any errors, is also added to a timestamped log file, in the current directory, named `riboviz-YYYYMMDD-HHMMSS.log` (for example, `riboviz.20190926-002455.log`).

To check the exit code, run:

```console
$ echo $0
```

See [Exit codes](#exit-codes), below, for a complete list of exit codes.

### Troubleshooting: `This script needs to be run under Python 3`

This warning arises if you try and run `prep_riboviz` under Python 2. You can only run `prep_riboviz` with Python 3.

### Troubleshooting: `samtools sort: couldn't allocate memory for bam_mem`

If using more than one process (`num_processes` > 1) you might get the error:

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

Divide the free memory by the number of processes, `num_processes` e.g. 1024/4 = 256 MB.

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

While `prep_riboviz` gives a warning, it will continue execution.

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

## Invoking `prep_riboviz` from outwith the RiboViz home directory

To invoke `prep_riboviz` from outwith the RiboViz home directory:

* Ensure that all the input paths in the YAML configuration file are correctly configured:
  - `dir_in`
  - `rrna_fasta_file`
  - `orf_fasta_file`
  - `orf_gff_file`
* Ensure that all the output paths in the YAML configuration file are correctly configured:
  - `dir_index`
  - `dir_tmp`
  - `dir_out`
  - `dir_logs`
* Ensure that all the data paths in the YAML configuration file are correctly configured:
  - `features_file`
  - `codon_positions_file`
  - `t_rna_file`
  - `asite_disp_length_file`

* Run:

```console
$ PYTHONPATH=<RIBOVIZ> python -m riboviz.tools.prep_riboviz -c <CONFIG_FILE>
```

where:

* `<RIBOVIZ>`: path to RiboViz home directory, which can be relative to the current directory or absolute.
* `<CONFIG>`: path to a YAML configuration file, which can be relative to the current directory or absolute.

For example, given an `analysis/` directory with contents:

```
sample_config.yaml
sample_input/
  sample.fastq
  yeast_rRNA_R64-1-1.fa
  yeast_YAL_CDS_w_250utrs.fa
  yeast_YAL_CDS_w_250utrs.gff3
```

and assuming that `analysis/` is a sibling of `riboviz/`, then file paths in `sample_config.yaml` can be configured as follows:

```yaml
dir_in: analysis/sample_input
dir_out: analysis/sample_output
dir_tmp: analysis/sample_tmp
dir_logs: analysis/sample_logs
rrna_fasta_file: analysis/sample_input/yeast_rRNA_R64-1-1.fa
orf_fasta_file: analysis/sample_input/yeast_YAL_CDS_w_250utrs.fa
orf_gff_file: analysis/sample_input/yeast_YAL_CDS_w_250utrs.gff3
dir_index: analysis/index
features_file: riboviz/data/yeast_features.tsv
t_rna_file: riboviz/data/yeast_tRNAs.tsv
codon_positions_file: riboviz/data/yeast_codon_pos_i200.RData
asite_disp_length_file: riboviz/data/yeast_standard_asite_disp_length.txt
```

`prep_riboviz` can then be run in the parent directory of `analysis/` and `riboviz/` as follows:

```console
$ PYTHONPATH=riboviz/ python -m riboviz.tools.prep_riboviz \
  -c analysis/sample_config.yaml 
```

Absolute paths can also be used. For example, if the file paths were as follows:

```yaml
dir_in: /home/user/analysis/sample_input
dir_out: /home/user/analysis/sample_output
dir_tmp: /home/user/analysis/sample_tmp
dir_logs: /home/user/analysis/sample_logs
rrna_fasta_file: /home/user/analysis/sample_input/yeast_rRNA_R64-1-1.fa
orf_fasta_file: /home/user/analysis/sample_input/yeast_YAL_CDS_w_250utrs.fa
orf_gff_file: /home/user/analysis/sample_input/yeast_YAL_CDS_w_250utrs.gff3
dir_index: /home/user/analysis/index
features_file: /home/user/riboviz/data/yeast_features.tsv
t_rna_file: /home/user/riboviz/data/yeast_tRNAs.tsv
codon_positions_file: /home/user/riboviz/data/yeast_codon_pos_i200.RData
asite_disp_length_file: /home/user/riboviz/data/yeast_standard_asite_disp_length.txt
```

`prep_riboviz` can then be run from any directory as follows:

```console
$ PYTHONPATH=/home/user/riboviz/ python -m riboviz.tools.prep_riboviz \
  -c /home/user/analysis/sample_config.yaml 
```

---

## Managing your disk usage

`prep_riboviz` generates many intermediate files and some of these may be unompressed and **large**, i.e. about the same size as the input files. All these files are placed in a temporary directory (`dir_tmp`). The temporary directory's contents can be inspected for troubleshooting, if necessary.

Similarly, `prep_riboviz` creates numerous log files.

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

When run in `--dry-run` mode or live, `prep_riboviz` logs the commands it submits to bash to execute local scripts and third-party tools. These are logged into a bash script (`cmd_file`). This script can be used for debugging. However it can also be run as-is and the effect is the same as running, or rerunning, `prep_riboviz`.

For example, if `cmd_file` had value `run_riboviz.sh` then the script can be run as follows:

```console
$ bash run_riboviz.sh
```

### Capturing bash commands and demultiplexing 

When run with demultiplexing enabled, `prep_riboviz` determines the demultiplexed files to process based on the samples in the `num_reads.tsv` file, output by `riboviz.tools.demultiplex_fastq`, which have one or more reads.

If running `prep_riboviz` in dry run mode then this file does not exist. In this case, `prep_riboviz` will determine the samples based on the samples in the sample sheet.

This means that the bash script output by `prep_riboviz` with `--dry-run` may not match that output during a live run, if, during the live run, some samples had 0 associated reads.

---

## Customising logging

You can customise logging by editing the file `riboviz/logging.yaml`.

### Removing timestamps

If you do not want `riboviz.log` to include a timestamp (i.e. you want the log file to be named `riboviz.log`) then edit `riboviz/logging.yaml` and replace:

```yaml
  handlers: [console, timestamp_file_handler]
```

with:

```yaml  
  handlers: [console, file_handler]
```

### Using custom log configuration files

A custom log configuration file can be provided by defining a `RIBOVIZ_LOG_CONFIG` environment variable. For example:

```console
$ RIBOVIZ_LOG_CONFIG=custom_logging.yaml
```

---

## Exit codes

`prep_riboviz` returns the following exit codes:

* 0: Processing successfully completed.
* 1: A file does not seem to exist.
* 2: Errors occurred loading or accessing configuration e.g. missing configuration parameters, inconsistent configuration parameters.
* 3: Error occurred during processing.
* 4: User is using Python 2.                          
