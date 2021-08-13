# Running the RiboViz Nextflow workflow

This page describes how to run the Nextflow workflow, `prep_riboviz.nf`. It is assumed you are familiar with [Configuring the RiboViz workflow](./prep-riboviz-config.md).

Contents:

* [Prepare input data](#prepare-input-data)
* [Set up your environment](#set-up-your-environment)
* [Configure number of processes (optional)](#configure-number-of-processes-optional)
* [Defining values for environment variables](#defining-values-for-environment-variables)
* [Validate configuration](#validate-configuration)
  - [Skip checks for ribosome profiling data files parameter](#skip-checks-for-ribosome-profiling-data-files-parameter)
* [Run the Nextflow workflow](#run-the-nextflow-workflow)
  - [Troubleshooting: `Error executing process staticHTML` and `AnalysisOutputs.Rmd`](#troubleshooting-error-executing-process-statichtml-and-analysisoutputsrmd)
  - [Troubleshooting: `samtools sort: couldn't allocate memory for bam_mem`](#troubleshooting-samtools-sort-couldnt-allocate-memory-for-bam_mem)
  - [Troubleshooting: deduplication and memory issues](#troubleshooting-deduplication-and-memory-issues)
  - [Troubleshooting: cutadapt 3.x failure](#troubleshooting-cutadapt-3x-failure)
  - [Troubleshooting: `ln: failed to create symbolic link`](#troubleshooting-ln-failed-to-create-symbolic-link)
* [Help](#help)
* [Incremental build and resuming a workflow](#incremental-build-and-resuming-a-workflow)
  - [Resuming after a failure](#resuming-after-a-failure)
  - [Resuming after adding more samples](#resuming-after-adding-more-samples)
  - [Resuming after changing configuration](#resuming-after-changing-configuration)
* [Multiplexed files](#multiplexed-files)
* [Debugging](#debugging)
* [Generating reports](#generating-reports)
  - [Troubleshooting: WARN To render the execution DAG in the required format it is required to install Graphviz](#troubleshooting-warn-to-render-the-execution-dag-in-the-required-format-it-is-required-to-install-graphviz)
* [Invoking the workflow from outwith the RiboViz home directory](#invoking-the-workflow-from-outwith-the-riboviz-home-directory)

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

* If you followed [Create `set-riboviz-env.sh` to configure paths](./install.md#create-set-riboviz-envsh-to-configure-paths), then run:

```console
$ source $HOME/set-riboviz-env.sh
```

* Otherwise, run the following (your directory names may be different, depending on the versions of Hisat2 and Bowtie you have):

```console
$ export PATH=~/hisat2-2.1.0:$PATH
$ export PATH=~/bowtie-1.2.2-linux-x86_64/:$PATH
```

---

## Configure number of processes (optional)

The workflow can instruct the tools it invokes to use a specific number of processes. By default this is 1.

To configure the workflow to use additional processes, in the configuration file, change:

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

## Defining values for environment variables

To specify values for environment variables cited as tokens in configuration parameters, (see [Environment variables and configuration tokens](./prep-riboviz-config.md#environment-variables-and-configuration-tokens)), you have two options, where:

* `<SAMPLES_DIRECTORY>` is a directory with input files.
* `<ORGANISMS_DIRECTORY>` is a directory with input files.
* `<DATA_DIRECTORY>` is a directory with input files.

The options are:

1. Specify environment variables with the paths to the directories on the same line as your command to run the workflow. The values will be used for this run of the workflow only. For example:

```console
$ RIBOVIZ_SAMPLES=<SAMPLES_DIRECTORY> \
  RIBOVIZ_ORGANISMS=<ORGANISMS_DIRECTORY> \
  RIBOVIZ_DATA=<DATA_DIRECTORY> \
  nextflow run prep_riboviz.nf -params-file <CONFIG_FILE>
```

2. Define values for the environment variables within your bash shell. The values will be available for successive runs of the workflow. The values need to be defined using `export` so they are available to `nextflow` when it runs. For example:

```console
$ export RIBOVIZ_SAMPLES=<SAMPLES_DIRECTORY>
$ export RIBOVIZ_ORGANISMS=<ORGANISMS_DIRECTORY>
$ export RIBOVIZ_DATA=<DATA_DIRECTORY>
$ nextflow run prep_riboviz.nf -params-file <CONFIG_FILE>
```

The above approaches can be combined i.e. you can define variables using `export` (2) but provide other values as part of the command to run the workflow (1). Values provided within the command take precedence over those defined via `export`.

**Note:** If using (1) then `<CONFIG_FILE>` cannot use the value of any variable defined within the same line. For example, if there was a configuration file `/home/$USER/riboviz/example-datasets/simulated/mok/Mok-tinysim_config.yaml` and the user runs:

```console
$ RIBOVIZ_SAMPLES=/home/$USER/tinysim \
  RIBOVIZ_DATA=/home/$USER/riboviz/riboviz/data \
  RIBOVIZ_ORGANISMS=/home/$USER/riboviz/example-datasets/simulated/mok \
  nextflow prep_riboviz.nf -params-file ${RIBOVIZ_ORGANISMS}/Mok-tinysim_config.yaml
```

then this will fail as `-params-file` has value `/Mok-tinysim_config.yaml` and not `/home/$USER/riboviz/example-datasets/simulated/mok/Mok-tinysim_config.yaml`. This is because the bash shell expands all the environment variables in the command (within `${...}`) *before* it runs the command (and, so, before the variables are defined). In such cases, either provide the path to `<CONFIG_FILE>` or define the variables using `export`.

**Note:** If a configuration file contains environment variable tokens then you **must** provide values for these when running the workflow.

---

## Validate configuration

The workflow supports a `--validate_only` command-line parameter which allows for the workflow configuration to be validated without running the workflow.

**Tip:** we strongly recommend validating the configuration before doing a live run on data you have not processed before.

Validate configuration, by running one of the following (depending on whether or not your workflow configuration uses environment variable tokens - see [Defining values for environment variables](#defining-values-for-environment-variables) above):

```console
$ nextflow run prep_riboviz.nf -params-file <CONFIG_FILE> --validate_only
```
```console
$ RIBOVIZ_SAMPLES=<SAMPLES_DIRECTORY> \
  RIBOVIZ_ORGANISMS=<ORGANISMS_DIRECTORY> \
  RIBOVIZ_DATA=<DATA_DIRECTORY> \
  nextflow run prep_riboviz.nf -params-file <CONFIG_FILE> --validate_only
```
```console
$ export RIBOVIZ_SAMPLES=<SAMPLES_DIRECTORY>
$ export RIBOVIZ_ORGANISMS=<ORGANISMS_DIRECTORY>
$ export RIBOVIZ_DATA=<DATA_DIRECTORY>
$ nextflow run prep_riboviz.nf -params-file <CONFIG_FILE> --validate_only
```

where:

* `<CONFIG_FILE>`: is a YAML configuration file.

### Skip checks for ribosome profiling data files parameter

`--validate_only` can be complemented by a `--skip_inputs` command-line parameter. This skips checks for the existence of the ribosome profiling data files (`fq_files`, `multiplexed_fq_files`, `sample_sheet`). An example without `--skip_inputs` might appear as:

```console
$ nextflow run prep_riboviz.nf -params-file vignette/experiment_config.yaml -ansi-log false --validate_only 
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [elated_bartik] - revision: e6bda28069
Validating configuration only
No such sample file (WTone): SRR1234_s1mi.fastq.gz
No such sample file (WTtwo): SRR5678_s1mi.fastq.gz
None of the defined sample files (fq_files) exist
```

And with `--skip_inputs` this might appear as:

```console
$ nextflow run prep_riboviz.nf -params-file vignette/experiment_config.yaml -ansi-log false --validate_only --skip_inputs
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [compassionate_galileo] - revision: e6bda28069
Validating configuration only
Skipping checks of ribosome profiling input files (fq_files|multiplex_fq_files
Validated configuration
```

This can be useful if you want to check that a configuration file you have received or downloaded is valid, before actually downloading, or preparing, the ribosome profiling data files.

---

## Run the Nextflow workflow

Run one of the following (depending on whether or not your workflow configuration uses environment variable tokens - see [Defining values for environment variables](#defining-values-for-environment-variables) above):

```console
$ nextflow run prep_riboviz.nf -ansi-log false -params-file <CONFIG_FILE>
```
```console
$ RIBOVIZ_SAMPLES=<SAMPLES_DIRECTORY> \
  RIBOVIZ_ORGANISMS=<ORGANISMS_DIRECTORY> \
  RIBOVIZ_DATA=<DATA_DIRECTORY> \
  nextflow run prep_riboviz.nf -ansi-log false -params-file <CONFIG_FILE>
```
```console
$ export RIBOVIZ_SAMPLES=<SAMPLES_DIRECTORY>
$ export RIBOVIZ_ORGANISMS=<ORGANISMS_DIRECTORY>
$ export RIBOVIZ_DATA=<DATA_DIRECTORY>
$ nextflow run prep_riboviz.nf -ansi-log false -params-file <CONFIG_FILE>
```

where:

* `-ansi-log false`: requests that each invocation of a Nextflow process is displayed on a separate line.
* `<CONFIG_FILE>`: is a YAML configuration file.
* Configuration parameters can also be provided via the command-line in the form `--<PARAMETER>=<VALUE>` (for example `--make_bedgraph=FALSE`).

The workflow will then execute, displaying information on each step as it is executed:

* Indexing steps are labelled with the index prefix
* Sample-specific steps are labelled with the sample ID.
* Multiplexed file-specific steps are labelled with the file name (minus extension).
* `collateTpms` is labelled with the IDs of all the samples that are collated.

### Troubleshooting: `Error executing process staticHTML` and `AnalysisOutputs.Rmd`

If your version of RiboViz is from a `.zip` or `.tar.gz` file downloaded from GitHub or Figshare (i.e. not a clone of the repository from GitHub) then you may see the following error during execution of a `staticHTML` step of the workflow. For example:

```console
...
[a1/93adf8] Submitted process > staticHTML (WT3AT)
Error executing process > 'staticHTML (WT3AT)'

Caused by:
  Process `staticHTML (WT3AT)` terminated with an error exit status (1)

Command executed:

  Rscript -e "rmarkdown::render('/home/user/riboviz/riboviz/rmarkdown/AnalysisOutputs.Rmd', ...)

Command exit status:
  1

...

Command error:
 
processing file: AnalysisOutputs.Rmd
  Quitting from lines 95-103 (AnalysisOutputs.Rmd) 
  Error in file(filename, "r", encoding = encoding) : 
    cannot open the connection
  Calls: <Anonymous> ... withCallingHandlers -> withVisible -> eval -> eval -> source -> file
  Execution halted
...

Workflow finished! (failed)
...
```

If this arises then create a `.here` file in your RiboViz directory:

```console
$ touch .here
```

and rerun the workflow.

If the problem still arises then it is recommended you use a clone of the `riboviz` repository from GitHub. Please also add information to the issue "Fix bug with AnalysisOutputs.Rmd use of 'here' if `.git` is not present" [#352](https://github.com/riboviz/riboviz/issues/352). It is unclear why this arises but it is hoped, when RiboViz moves its R code into a package that this problem will no longer arise.

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

Either, rerun the workflow and provide the memory via the command line e.g. `--samsort_memory=256M`.

Or, edit your YAML configuration, and add the line:

```yaml
samsort_memory: 256M
```

For information on the allowable values, see [samtools sort](http://www.htslib.org/doc/samtools-sort.html) and its `-m` flag.

### Troubleshooting: deduplication and memory issues

See [Troubleshooting: deduplication and memory issues](memory-storage.md#troubleshooting-deduplication-and-memory-issues) in [Memory and storage](./memory-storage.md).

### Troubleshooting: cutadapt 3.x failure

For some users, the workflow fails during invocation of `cutadapt` 3.2 or 3.4. [Debugging](#debugging) reveals the following error message:

```
Traceback (most recent call last):  
...
-  Traceback (most recent call last):  File "/.../miniconda3/envs/riboviz/bin/cutadapt", line 12, in <module>
    sys.exit(main())TypeError: main() missing 1 required positional argument: 'cmdlineargs'
```

It is unclear why this arises but it seems to arise if the Python `cutadapt` package was installed using `conda`. One solution is to reinstall the Python `cutadapt` package using `pip`. For example:

```console
$ conda remove cutadapt
$ pip install cutadapt
```

If you encounter this issue for another version of cutadapt or if you have other solutions then please add a comment to cutadapt 3.x from conda causes workflow failure [#380](https://github.com/riboviz/riboviz/issues/380).

### Troubleshooting: `ln: failed to create symbolic link`

If the workflow fails and [Debugging](#debugging) reveals an error message of form:

```
ln: failed to create symbolic link '<FILE>': Operation not supported
```

then this can arise if you do not have permission to create symbolic links within your file system. Nextflow requires the ability to be able to create symbolic links.

The solution is to reconfigure your file system to allow you to create symbolic links within it.

---

## Help

Usage and configuration information can be viewed via use of the `--help` flag:

```console
$ nextflow run prep_riboviz.nf --help
```

Note that `--help` displays RiboViz-specific workflow help, whereas `-help` display's the `nextflow run` command's in-built help.

---

## Incremental build and resuming a workflow

Nextflow supports a `-resume` parameter which allows a workflow to be resumed from the point at which it failed or when new samples have been added.

### Resuming after a failure

If the processing of a sample fails then the workflow has been written to ensure that this does not prevent the processing of other samples. For example, if the file for the vignette sample `WTnone` was corrupt in some way, then the workflow might fail as follows:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [mad_heisenberg] - revision: 134fcb8f6f
No such file (NotHere): example_missing_file.fastq.gz
[a4/e3fc5d] Submitted process > cutAdapters (WT3AT)
[03/f1cfc3] Submitted process > cutAdapters (WTnone)
[9d/43b3e0] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[ec/1aafe5] Submitted process > buildIndicesrRNA (yeast_rRNA)
[03/f1cfc3] NOTE: Process `cutAdapters (WTnone)` terminated with an error exit status (1) -- Error is ignored
[35/1c7ac8] Submitted process > hisat2rRNA (WT3AT)
[1a/af8ff1] Submitted process > hisat2ORF (WT3AT)
[c4/2f18a2] Submitted process > trim5pMismatches (WT3AT)
...
```

After fixing the problematic file, the workflow can be resumed, using `-resume` as follows:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false -resume
```

The workflow will rerun but the samples already processed will not be reprocessed. The outputs to date for these samples, cached by Nextflow, will be reused.

### Resuming after adding more samples

`-resume` can also be used to resume a workflow, when more samples have been added to a configuration, without reprocessing the samples already processed. For example, given a `vignette_config.yaml` which specifies only sample `WTnone`, running Nextflow may give:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [serene_noether] - revision: 26bb98ec7b
No such file (NotHere): example_missing_file.fastq.gz
[bf/ddd7de] Submitted process > cutAdapters (WTnone)
[0d/f82601] Submitted process > buildIndicesrRNA (yeast_rRNA)
[61/5a1186] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[02/dac1f8] Submitted process > hisat2rRNA (WTnone)
[15/e5f9ed] Submitted process > hisat2ORF (WTnone)
[ed/fc50b9] Submitted process > trim5pMismatches (WTnone)
...
```

If `WT3AT` is then added to `vignette_config.yaml` and Nextflow is run with the `-resume` option, then only the processing for `WT3AT` is done, the outputs to date for `WTnone`, cached by Nextflow, are reused:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false -resume
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [desperate_coulomb] - revision: 26bb98ec7b
No such file (NotHere): example_missing_file.fastq.gz
[80/79a5f4] Submitted process > cutAdapters (WT3AT)
[bf/ddd7de] Cached process > cutAdapters (WTnone)
[61/5a1186] Cached process > buildIndicesORF (YAL_CDS_w_250)
[0d/f82601] Cached process > buildIndicesrRNA (yeast_rRNA)
[02/dac1f8] Cached process > hisat2rRNA (WTnone)
[15/e5f9ed] Cached process > hisat2ORF (WTnone)
[ed/fc50b9] Cached process > trim5pMismatches (WTnone)
[e3/9c39d7] Submitted process > hisat2rRNA (WT3AT)
[a6/49addd] Submitted process > hisat2ORF (WT3AT)
[46/d5db2c] Submitted process > trim5pMismatches (WT3AT)
...
```

### Resuming after changing configuration

`-resume` can also be used to trigger a rerun of parts of the workflow when you have changed the configuration. For example, if the configuration file has:

```yaml
make_bedgraph: FALSE
```

then bedgraphs will not be created when the workflow is run. If you decide you do want the bedgraphs after all you can use `-resume`, in conjunction with setting the `make_bedgraph` parameter to `TRUE`, on the command-line, to create the bedgraphs. For example:

```console
$ nextflow run prep_riboviz.nf -params-file vignette/vignette_config.yaml -ansi-log false -resume --make_bedgraph=TRUE
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [loving_fermi] - revision: 6df32df2e2
No such file (NotHere): example_missing_file.fastq.gz
[51/2f6008] Cached process > cutAdapters (WTnone)
[fe/9f150a] Cached process > cutAdapters (WT3AT)
[4c/71f6ea] Cached process > buildIndicesORF (YAL_CDS_w_250)
[42/3c9bd2] Cached process > buildIndicesrRNA (yeast_rRNA)
[2e/8a1462] Cached process > hisat2rRNA (WT3AT)
[6d/06df71] Cached process > hisat2rRNA (WTnone)
[77/218316] Cached process > hisat2ORF (WT3AT)
[e0/f45966] Cached process > hisat2ORF (WTnone)
[8f/c92188] Cached process > trim5pMismatches (WT3AT)
[bd/87fae1] Cached process > trim5pMismatches (WTnone)
[bb/598220] Cached process > samViewSort (WTnone)
[3b/3a09d2] Cached process > samViewSort (WT3AT)
[dc/a27468] Cached process > outputBams (WTnone)
[f9/3528ca] Cached process > outputBams (WT3AT)
[5a/bef907] Cached process > bamToH5 (WTnone)
[89/cf27fd] Cached process > bamToH5 (WT3AT)
[74/a2f2aa] Cached process > generateStatsFigs (WTnone)
Finished processing sample: WTnone
[c6/bfef26] Cached process > generateStatsFigs (WT3AT)
Finished processing sample: WT3AT
[b0/b534d4] Submitted process > makeBedgraphs (WTnone)
[80/2638f2] Cached process > renameTpms (WTnone)
[5f/4db643] Cached process > renameTpms (WT3AT)
[ff/5a41cb] Submitted process > makeBedgraphs (WT3AT)
[12/8109a2] Cached process > collateTpms (WTnone, WT3AT)
[0c/9d6aa8] Cached process > countReads
Workflow finished! (OK)
```

This works because the workflow uses `make_bedgraph` to determine whether or not to create bedgraphs. In some cases, changing the configuration will not be enough to trigger a reinvocation of the affected steps. For example, changing the value of `adapters` then rerunning the workflow with `-resume` does not trigger a reinvocation of `cutadapt`. For these cases you need to remove the output files produced by the step you want to rerun. When you resume the workflow, Nextflow will notice that these files are missing and reinvoke the step to recreate them, using your updated configuration.

To remove the output files for a specific step of the workflow and rerun the workflow to recreate these files using an updated configuration:

* Use the `nextflow log` command to find the directory in Nextflow's `work` directory for the step. For more information see [Nextflow `work/` directory](./prep-riboviz-operation.md#nextflow-work-directory) and [Debugging](./prep-riboviz-run-nextflow.md#debugging).
* Remove the output files from the step-specific `work/` subdirectory.
* Update the configuration.
* `-resume` the workflow.

---

## Multiplexed files

When running the workflow with a multiplexed FASTQ file (`multiplex_fq_files`) then a list of the sample IDs for which there were matching reads (and output files) and those for which there were no matching reads (and no output files) is printed.

For example, the sample sheet `data/simdata/multiplex_barcodes.tsv`, used in `vignette/simdata_multiplex_config.yaml` specifies a sample, `Tag3`, for which there will be no matching reads. The output of the workflow is as follows:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/simdata_multiplex_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [elated_jones] - revision: 0e94b1fa62
[e9/53999d] Submitted process > cutAdaptersMultiplex (multiplex_umi_barcode_adaptor)
[c3/6566db] Submitted process > buildIndicesrRNA (yeast_rRNA)
[43/752a9a] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[19/981102] Submitted process > extractUmisMultiplex (multiplex_umi_barcode_adaptor)
[e5/c52b82] Submitted process > demultiplex (multiplex_umi_barcode_adaptor)
Demultiplexed samples: [Tag0, Tag1, Tag2]
Non-demultiplexed samples: [Tag3]
[3f/73d0ab] Submitted process > hisat2rRNA (Tag2)
[18/c710a7] Submitted process > hisat2rRNA (Tag1)
[cc/f425d1] Submitted process > hisat2rRNA (Tag0)
[25/dd22c3] Submitted process > hisat2ORF (Tag2)
[b3/2db101] Submitted process > hisat2ORF (Tag1)
[85/9427ed] Submitted process > hisat2ORF (Tag0)
[f9/362265] Submitted process > trim5pMismatches (Tag0)
[d2/082ba5] Submitted process > trim5pMismatches (Tag2)
[de/e42cb0] Submitted process > trim5pMismatches (Tag1)
[63/1da506] Submitted process > samViewSort (Tag0)
[a5/ab590b] Submitted process > samViewSort (Tag2)
[72/4e9acd] Submitted process > samViewSort (Tag1)
[8b/c0936d] Submitted process > dedupUmis (Tag2)
[46/e9f43f] Submitted process > groupUmisPreDedup (Tag0)
[07/d696b3] Submitted process > groupUmisPreDedup (Tag2)
[52/2bdb04] Submitted process > dedupUmis (Tag0)
[36/ef62b5] Submitted process > dedupUmis (Tag1)
[5f/f51fca] Submitted process > groupUmisPreDedup (Tag1)
[9d/688af6] Submitted process > outputBams (Tag2)
[dc/e5d588] Submitted process > groupUmisPostDedup (Tag2)
[c5/690b4c] Submitted process > groupUmisPostDedup (Tag0)
[85/241fd3] Submitted process > outputBams (Tag0)
[9c/d88da8] Submitted process > bamToH5 (Tag2)
[d0/316e3d] Submitted process > makeBedgraphs (Tag2)
[a4/c2d620] Submitted process > makeBedgraphs (Tag0)
[0e/e327cb] Submitted process > bamToH5 (Tag0)
[62/789bb7] Submitted process > groupUmisPostDedup (Tag1)
[c7/13ea6b] Submitted process > outputBams (Tag1)
[8f/dabdca] Submitted process > bamToH5 (Tag1)
[d4/b4eef0] Submitted process > makeBedgraphs (Tag1)
[f8/6ba638] Submitted process > generateStatsFigs (Tag2)
[0a/06981b] Submitted process > generateStatsFigs (Tag1)
[f6/27751b] Submitted process > generateStatsFigs (Tag0)
Finished processing sample: Tag0
[52/447fbe] Submitted process > renameTpms (Tag0)
Finished processing sample: Tag1
[10/54496e] Submitted process > renameTpms (Tag1)
Finished processing sample: Tag2
[d7/7cc9ea] Submitted process > renameTpms (Tag2)
[6c/28e9e0] Submitted process > collateTpms (Tag0, Tag1, Tag2)
[29/354049] Submitted process > countReads
Workflow finished! (OK)
```

If no samples in the sample sheet have a matching read (and no sample-specific files are output at all) then the output will appear as follows:

```console
...
[53/c42278] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[47/59eaa1] Submitted process > buildIndicesrRNA (yeast_rRNA)
[58/24a15e] Submitted process > cutAdaptersMultiplex (multiplex_umi_barcode_adaptor)
[f8/df225f] Submitted process > extractUmisMultiplex (multiplex_umi_barcode_adaptor)
[53/8df442] Submitted process > demultiplex (multiplex_umi_barcode_adaptor)
Demultiplexed samples: []
Non-demultiplexed samples: [Tag0, Tag1, Tag2, Tag3]
Workflow finished! (OK)
```

---

## Debugging

Nextflow prints information if one or more steps fail. For example:

```console
$ nextflow run -ansi-log false prep_riboviz.nf -params-file vignette/vignette_config.yaml
N E X T F L O W  ~  version 20.04.1
Launching `prep_riboviz.nf` [big_majorana] - revision: 6c6670470d
samples_dir: .
organisms_dir: .
data_dir: .
No such sample file (NotHere): example_missing_file.fastq.gz
[7e/e093a6] Submitted process > cutAdapters (WT3AT)
[42/b1c9a2] Submitted process > buildIndicesrRNA (yeast_rRNA)
[3d/5bcc67] Submitted process > cutAdapters (WTnone)
[0f/0757fc] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[e4/850912] Submitted process > createVizParamsConfigFile
[e9/5f318c] Submitted process > hisat2rRNA (WTnone)
[8d/1ecb5d] Submitted process > hisat2rRNA (WT3AT)
[f5/a34897] Submitted process > hisat2ORF (WTnone)
[11/a2c59f] Submitted process > trim5pMismatches (WTnone)
[7e/2579c7] Submitted process > samViewSort (WTnone)
[d2/68d6aa] Submitted process > outputBams (WTnone)
[53/5bee02] Submitted process > makeBedgraphs (WTnone)
[16/26d72f] Submitted process > bamToH5 (WTnone)
[2d/0330b5] Submitted process > hisat2ORF (WT3AT)
[02/5594ad] Submitted process > trim5pMismatches (WT3AT)
[91/3ad140] Submitted process > samViewSort (WT3AT)
[60/796bd2] Submitted process > outputBams (WT3AT)
[e2/254948] Submitted process > makeBedgraphs (WT3AT)
[14/154ad6] Submitted process > bamToH5 (WT3AT)
[e1/bebe80] Submitted process > generateStatsFigs (WTnone)
[b0/013c22] Submitted process > generateStatsFigs (WT3AT)
Finished processing sample: WTnone
[69/1f6794] Submitted process > renameTpms (WTnone)
[a1/b61f46] Submitted process > staticHTML (WTnone)
Finished processing sample: WT3AT
[a1/6eb17a] Submitted process > renameTpms (WT3AT)
[b0/5a6872] Submitted process > staticHTML (WT3AT)
Finished visualising sample: WTnone
Finished visualising sample: WT3AT
[38/784d89] Submitted process > collateTpms (WTnone, WT3AT)
Error executing process > 'collateTpms (WTnone, WT3AT)'

Caused by:
  Process `collateTpms (WTnone, WT3AT)` terminated with an error exit status (1)

Command executed:

  Rscript --vanilla /home/ubuntu/riboviz/rscripts/collate_tpms.R             --tpms-file=TPMs_all_CDS_all_samples.tsv             WTnone WTnone_tpms.tsv WT3AT WT3AT_tpms.tsv

Command exit status:
  1

Command output:
  ...
  [1] "collate_tpms.R running with parameters:"
  $options
  $options$tpms_file
  [1] "TPMs_all_CDS_all_samples.tsv"
  
  $options$orf_fasta
  [1] NA
  
  $options$sort_orfs
  [1] FALSE
  
  $options$digits
  [1] 1
  
  $options$help
  [1] FALSE
  
  
  $args
  [1] "WTnone"          "WTnone_tpms.tsv" "WT3AT"           "WT3AT_tpms.tsv" 
  
  [1] "Loading ORFs from: WTnone_tpms.tsv"
  [1] "Loading TPMs from: WTnone_tpms.tsv"
  [1] "Loading TPMs from: WT3AT_tpms.tsv"

Command error:
  Parsed with column specification:
  cols(
    ORF = col_character(),
    readcount = col_double(),
    rpb = col_double(),
    tpm = col_double()
  )
  Parsed with column specification:
  cols(
    ORF = col_character(),
    readcount = col_double(),
    rpb = col_double(),
    tpm = col_double()
  )
  Parsed with column specification:
  cols(
    ORF = col_character(),
    readcount = col_double(),
    rpb = col_double(),
    tpm = col_double()
  )
  Error: Can't recycle `ORF` (size 68) to match `WT3AT` (size 4).
  Backtrace:
       ¦
    1. +-global::CollateTpms(...)
    2. ¦ +-`%>%`(...)
    3. ¦ ¦ +-base::eval(lhs, parent, parent)
    4. ¦ ¦   +-base::eval(lhs, parent, parent)
    5. ¦ +-global::MakeTpmTable(orf_fasta, samples, sort_orfs = sort_orfs)
    6. ¦   +-dplyr::bind_cols(ORF = orfs, tpm_list[non_null_elts])
    7. ¦     +-vctrs::vec_cbind(!!!dots, .name_repair = .name_repair)
    8. +-vctrs::stop_incompatible_size(...)
    9.   +-vctrs:::stop_incompatible(...)
   10.     +-vctrs:::stop_vctrs(...)
  Execution halted

Work dir:
  /home/ubuntu/riboviz/work/38/784d89646ff067d5fa9bedcdd4db73

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

If there are many errors, this may result in lots of lines being printed, which can be overwhelming. Nextflow has a `nextflow log` command that allows for the browsing of step-specific information, for both successful and failed steps, in a more manageable way.

Each time Nextflow runs a workflow, it gives the workflow a unique name. For example, in the above, failed, run it was `big_majorana`:

```
Launching `prep_riboviz.nf` [big_majorana] - revision: 6c6670470d
```

`nextflow log` can be used with this unique name to find out information about the workflow's run.

To see every invocation of a process, every step, that was run and the names of these steps, run, for example:

```console
$ nextflow log big_majorana -f process,name
cutAdapters	cutAdapters (WT3AT)
buildIndicesrRNA	buildIndicesrRNA (yeast_rRNA)
cutAdapters	cutAdapters (WTnone)
buildIndicesORF	buildIndicesORF (YAL_CDS_w_250)
createVizParamsConfigFile	createVizParamsConfigFile
hisat2rRNA	hisat2rRNA (WTnone)
hisat2rRNA	hisat2rRNA (WT3AT)
hisat2ORF	hisat2ORF (WTnone)
trim5pMismatches	trim5pMismatches (WTnone)
samViewSort	samViewSort (WTnone)
outputBams	outputBams (WTnone)
makeBedgraphs	makeBedgraphs (WTnone)
bamToH5	bamToH5 (WTnone)
hisat2ORF	hisat2ORF (WT3AT)
trim5pMismatches	trim5pMismatches (WT3AT)
samViewSort	samViewSort (WT3AT)
outputBams	outputBams (WT3AT)
makeBedgraphs	makeBedgraphs (WT3AT)
bamToH5	bamToH5 (WT3AT)
generateStatsFigs	generateStatsFigs (WTnone)
generateStatsFigs	generateStatsFigs (WT3AT)
renameTpms	renameTpms (WTnone)
staticHTML	staticHTML (WTnone)
renameTpms	renameTpms (WT3AT)
staticHTML	staticHTML (WT3AT)
collateTpms	collateTpms (WTnone, WT3AT)
countReads	countReads
```

If the workflow failed, then the names of the failed steps can be displayed as follows, for example:

```console
$ nextflow log big_majorana -f name -filter "status == 'FAILED'"
collateTpms (WT3AT, WTnone)
```

To see only the names of specific steps, invocations of a specific process, run, for example:

```console
$ nextflow log big_majorana -f name -filter "process == 'cutAdapters'"
cutAdapters (WTnone)
cutAdapters (WT3AT)
```
```console
$ nextflow log big_majorana -f name -filter "process == 'collateTpms'"
collateTpms (WT3AT, WTnone)
```

To see information about a specific step including the command that was run by Nextflow, and its exit code, run, for example:

```console
$ nextflow log big_majorana -f script,exit -filter "name == 'cutAdapters (WTnone)'"

        cutadapt --trim-n -O 1 -m 5 -a CTGTAGGCACC             -o trim.fq SRR1042855_s1mi.fastq.gz -j 0
        	0
```
```console
$ nextflow log big_majorana -f script,exit -filter "name == 'collateTpms (WT3AT, WTnone)'"

        Rscript --vanilla /home/ubuntu/riboviz/rscripts/collate_tpms.R             --tpms-file=TPMs_all_CDS_all_samples.tsv             WT3AT WT3AT_tpms.tsv WTnone WTnone_tpms.tsv
        	1
```

Here, the `cutadapt` command run for `cutAdapters (WTnone)` had an exit code of 0 so the step succeeded. The `collate_tpms.R` command run for `collateTpms (WT3AT, WTnone)` had an exit code of 1, indicating an error, so the step failed.

To see any output and error messages printed by the commands run by Nextflow and which are captured by Nextflow, run, for example:

```console
$ nextflow log big_majorana -f stdout,stderr -filter "name == 'cutAdapters (WTnone)'"
This is cutadapt 1.18 with Python 3.7.6Command line parameters: --trim-n -O 1 -m 5 -a CTGTAGGCACC -o trim.fq SRR1042855_s1mi.fastq.gz -j 0Processing reads on 4 cores in single-end mode ...Finished in 7.35 s (8 us/read; 7.87 M reads/minute).-
```
```console
$ nextflow log big_majorana -f stdout,stderr -filter "name == 'collateTpms (WT3AT, WTnone)'"
[1] "Created by: RiboViz"[1] "Date: 2021-06-24 05:43:14"[1] "File: /home/ubuntu/riboviz/rscripts/collate_tpms.R"[1] "Version: commit 144d95ab228a2da71b9a92912a24b26c37f4a64e date 2021-06-21 08:11:03 GMT"[1] "collate_tpms.R running with parameters:"$options$options$tpms_file[1] "TPMs_all_CDS_all_samples.tsv"	Parsed with column specification:cols(  ORF = col_character(),  readcount = col_double(),  rpb = col_double(),  tpm = col_double())Parsed with column specification:cols(  ORF = col_character(),  readcount = col_double(),  rpb = col_double(),  tpm = col_double())Parsed with column specification:cols(  ORF = col_character(),  readcount = col_double(),  rpb = col_double(),  tpm = col_double())Error: Can't recycle `ORF` (size 4) to match `WTnone` (size 68).Backtrace:     ¦  1. +-global::CollateTpms(...)  2. ¦ +-`%>%`(...)  3. ¦ ¦ +-base::eval(lhs, parent, parent)  4. ¦ ¦   +-base::eval(lhs, parent, parent)  5. ¦ +-global::MakeTpmTable(orf_fasta, samples, sort_orfs = sort_orfs)  6. ¦   +-dplyr::bind_cols(ORF = orfs, tpm_list[non_null_elts])  7. ¦     +-vctrs::vec_cbind(!!!dots, .name_repair = .name_repair)  8. +-vctrs::stop_incompatible_size(...)  9.   +-vctrs:::stop_incompatible(...) 10.     +-vctrs:::stop_vctrs(...)Execution halted
```

Only the first few lines of the output are shown. If, as for `cutadapt (WTnone)`, there were no error messages, then `-` is displayed. The information may not be readable. How to view the complete record of output and error messages for a step is described below.

As described in [Nextflow `work/` directory](./prep-riboviz-operation.md#nextflow-work-directory), when Nextflow runs, it creates a unique step-specific directory for every step in the workflow. Each step-specific directory has symbolic links to the input files for the step and a bash script with the commands to be run by Nextflow for that step. Nextflow runs this bash script within this directory which creates the step's output files. Nextflow also creates files with output and error messages and the exit code. These step-specific directories are created within a Nextflow `work/` directory located, by default, within the directory within which Nextflow is run.

The location of the `work/` subdirectory for a step can be accessed as follows, for example:

```console
$ nextflow log big_majorana -f hash,workdir -filter "name == 'collateTpms (WT3AT, WTnone)'"

38/784d89	/home/ubuntu/riboviz/work/38/784d89646ff067d5fa9bedcdd4db73
```

This directory can then be browsed:

```console
$ ls -1a /home/ubuntu/riboviz/work/38/784d89646ff067d5fa9bedcdd4db73
.
..
.command.begin
.command.err
.command.log
.command.out
.command.run
.command.sh
.exitcode
TPMs_all_CDS_all_samples.tsv
WT3AT_tpms.tsv
WTnone_tpms.tsv
```

The complete list of output and error messages captured by Nextflow for the step can be viewed in the `.command.out` and `.command.err` files. These files can be viewed:

```console
$ cat /home/ubuntu/riboviz/work/38/784d89646ff067d5fa9bedcdd4db73/.command.out 
[1] "Created by: RiboViz"
[1] "Date: 2021-06-24 05:43:14"
[1] "File: /home/ubuntu/riboviz/rscripts/collate_tpms.R"
[1] "Version: commit 144d95ab228a2da71b9a92912a24b26c37f4a64e date 2021-06-21 08:11:03 GMT"
[1] "collate_tpms.R running with parameters:"
$options
$options$tpms_file
[1] "TPMs_all_CDS_all_samples.tsv"

$options$orf_fasta
[1] NA

$options$sort_orfs
[1] FALSE

$options$digits
[1] 1

$options$help
[1] FALSE


$args
[1] "WT3AT"           "WT3AT_tpms.tsv"  "WTnone"          "WTnone_tpms.tsv"

[1] "Loading ORFs from: WT3AT_tpms.tsv"
[1] "Loading TPMs from: WT3AT_tpms.tsv"
[1] "Loading TPMs from: WTnone_tpms.tsv"
```
```console
$ cat /home/ubuntu/riboviz/work/38/784d89646ff067d5fa9bedcdd4db73/.command.err 
Parsed with column specification:
cols(
  ORF = col_character(),
  readcount = col_double(),
  rpb = col_double(),
  tpm = col_double()
)
Parsed with column specification:
cols(
  ORF = col_character(),
  readcount = col_double(),
  rpb = col_double(),
  tpm = col_double()
)
Parsed with column specification:
cols(
  ORF = col_character(),
  readcount = col_double(),
  rpb = col_double(),
  tpm = col_double()
)
Error: Can't recycle `ORF` (size 4) to match `WTnone` (size 68).
Backtrace:
     ¦
  1. +-global::CollateTpms(...)
  2. ¦ +-`%>%`(...)
  3. ¦ ¦ +-base::eval(lhs, parent, parent)
  4. ¦ ¦   +-base::eval(lhs, parent, parent)
  5. ¦ +-global::MakeTpmTable(orf_fasta, samples, sort_orfs = sort_orfs)
  6. ¦   +-dplyr::bind_cols(ORF = orfs, tpm_list[non_null_elts])
  7. ¦     +-vctrs::vec_cbind(!!!dots, .name_repair = .name_repair)
  8. +-vctrs::stop_incompatible_size(...)
  9.   +-vctrs:::stop_incompatible(...)
 10.     +-vctrs:::stop_vctrs(...)
Execution halted
```

In this example, inspecting the input files reveals the problem:

```console
$ cat /home/ubuntu/riboviz/work/38/784d89646ff067d5fa9bedcdd4db73/WT3AT_tpms.tsv 
# Created by: RiboViz
# Date: 2021-06-24 03:20:00
# File: /home/ubuntu/riboviz/rscripts/generate_stats_figs.R
# Version: commit 144d95ab228a2da71b9a92912a24b26c37f4a64e date 2021-06-21 08:11:03 GMT
ORF	readcount	rpb	tpm
YAL001C	19	0.00538243626062323	953.022776609483
YAL002W	8	0.00206611570247934	365.829752221778
YAL003W	379	0.567365269461078	100458.602437955
YAL005C	2668	1.3502024291498	239069.002530875
$ cat /home/ubuntu/riboviz/work/38/784d89646ff067d5fa9bedcdd4db73/WTnone_tpms.tsv 
# Created by: RiboViz
# Date: 2021-06-24 03:19:51
# File: /home/ubuntu/riboviz/rscripts/generate_stats_figs.R
# Version: commit 144d95ab228a2da71b9a92912a24b26c37f4a64e date 2021-06-21 08:11:03 GMT
ORF	readcount	rpb	tpm
YAL001C	4	0.00113314447592068	116.685976611765
YAL002W	8	0.00206611570247934	212.759037933641
YAL003W	1291	1.93263473053892	199013.784878156
YAL005C	4797	2.42763157894737	249986.270649977
YAL007C	42	0.060431654676259	6222.9722634749
YAL008W	6	0.0093167701863354	959.397897762878
YAL009W	4	0.00483675937122128	498.06710694018
YAL010C	5	0.00327011118378025	336.740923348209
YAL011W	16	0.00831168831168831	855.899215458761
...
```

`collate_tpms.R` expects each input file to have the same number of rows and the same ORFs. In this case they differ. Specifically, `WT3AT_tpms.tsv` only has rows for four ORFs. This may indicate a problem with the data for sample WT3AT.

As the Nextflow `work/` subdirectory includes the bash script with the command that was run and symbolic links to any input files used by the step, this bash script can rerun as-is within the subdirectory to investigate a failure in more detail. For example:

```console
$ cd /home/ubuntu/riboviz/work/38/784d89646ff067d5fa9bedcdd4db73/
$ bash .command.sh 
[1] "Created by: RiboViz"
[1] "Date: 2021-06-24 06:08:35"
...
  9.   +-vctrs:::stop_incompatible(...)
 10.     +-vctrs:::stop_vctrs(...)
Execution halted
```

For more information on the `work/` directory, and its files, see [Nextflow `work/` directory](./prep-riboviz-operation.md#nextflow-work-directory),

For more information on `nextflow log`, see Nextflow's [Tracing & visualisation](https://www.nextflow.io/docs/latest/tracing.html) or run `nextflow log -h`. To see the log fields available, run `nextflow log -l`.

---

## Generating reports

Nextflow's `-with-report`, `-with-timeline`, `with-trace` and `-with-dag` flags allow you to request that Nextflow create reports on a run and an image of the execution workflow. For example:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false \
    -with-report report.html -with-timeline timeline.html \
    -with-trace trace.tsv -with-dag workflow.svg
```

For more information on Nextflow's reports, see Nextflow [Tracing & visualisation](https://www.nextflow.io/docs/latest/tracing.html).

**Note:** If the `--validate_only` flag is provided (see [Validate configuration](#validate-configuration)) then the execution report, timeline report and DAG files will be empty and the trace report will consist of a header only, but no data.

### Troubleshooting: `WARN To render the execution DAG in the required format it is required to install Graphviz`

If you do not have Graphviz available then requesting `-with-dag workflow.svg` will cause Nextflow to show this warning on completion of the workflow.

If this happens then Nextflow will, instead, create a `nextflow.dot`, file. This can be converted to an image on a system that does have Graphviz as follows:

```console
$ dot -T svg workflow.dot -o workflow.svg
```

Alternatively you could run Nextflow using `-with-dag workflow.html`, which creates an HTML/JavaScript visualisation of the workflow.

For more information see [DAG visualisation](https://www.nextflow.io/docs/latest/tracing.html#dag-visualisation) in the Nextflow documentation.

---

## Invoking the workflow from outwith the RiboViz home directory

The workflow can be invoked from any directory, by providing the path to the workflow. For example:

```console
$ nextflow run <RIBOVIZ>/prep_riboviz.nf -params-file <CONFIG_FILE> -ansi-log false
```

For example:

```console
$ nextflow run ../riboviz/prep_riboviz.nf -params-file <CONFIG_FILE> -ansi-log false
```

Note that if `<CONFIG_FILE>` has relative paths then these will be relative to the current directory, not `<CONFIG_FILE>` or `prep_riboviz.nf`.
