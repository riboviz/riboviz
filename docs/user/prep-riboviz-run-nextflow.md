# Running the RiboViz Nextflow workflow

This page describes how to run the Nextflow workflow, `prep_riboviz.nf`. It is assumed you are familiar with [Configuring the RiboViz workflow](./prep-riboviz-config.md).

Contents:

* [Prepare input data](#prepare-input-data)
* [Set up your environment](#set-up-your-environment)
* [Configure number of processes (optional)](#configure-number-of-processes-optional)
* [Validate configuration](#validate-configuration)
  - [Skip checks for ribosome profiling data files parameter](#skip-checks-for-ribosome-profiling-data-files-parameter)
* [Run the Nextflow workflow](#run-the-nextflow-workflow)
  - [Troubleshooting: `samtools sort: couldn't allocate memory for bam_mem`](#troubleshooting-samtools-sort-couldnt-allocate-memory-for-bam_mem)
  - [Troubleshooting: deduplication and memory issues](#troubleshooting-deduplication-and-memory-issues)
* [Help](#help)
* [Incremental build](#incremental-build)
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

* If you followed [Create `setenv.sh` to configure paths](./install.md#create-setenvsh-to-configure-paths), then run:

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

## Validate configuration

The workflow supports a `--validate_only` command-line parameter which allows for the workflow configuration to be validated without running the workflow.

**Tip:** we strongly recommend validating the configuration before doing a live run on data you have not processed before.

Validate configuration, by running:

```console
$ nextflow run prep_riboviz.nf -params-file <CONFIG_FILE> --validate_only
```

where:

* `<CONFIG_FILE>`: is a YAML configuration file.

To specify values for environment variables cited as tokens in configuration parameters, (see [Environment variables and configuration tokens](./prep-riboviz-config.md#environment-variables-and-configuration-tokens)), then these should be defined in the bash shell within which the workflow is run. Alternatively, they can be provided when running the workflow:

```console
$ RIBOVIZ_SAMPLES=<SAMPLES_DIRECTORY> \
  RIBOVIZ_ORGANISMS=<ORGANISMS_DIRECTORY> \
  RIBOVIZ_DATA=<DATA_DIRECTORY> \
  nextflow run prep_riboviz.nf -params-file <CONFIG_FILE> --validate_only
```

where:

* `<SAMPLES_DIRECTORY>` is a directory with input files.
* `<ORGANISMS_DIRECTORY>` is a directory with input files.
* `<DATA_DIRECTORY>` is a directory with input files.

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

Run:

```console
$ nextflow run prep_riboviz.nf -ansi-log false -params-file <CONFIG_FILE>
```

where:

* `-ansi-log false`: requests that each invocation of a Nextflow task is displayed on a separate line.
* `<CONFIG_FILE>`: is a YAML configuration file.
* Configuration parameters can also be provided via the command-line in the form `--<PARAMETER>=<VALUE>` (for example `--make_bedgraph=FALSE`).

To specify values for environment variables cited as tokens in configuration parameters, (see [Environment variables and configuration tokens](./prep-riboviz-config.md#environment-variables-and-configuration-tokens)), then these should be defined in the bash shell within which the workflow is run. Alternatively, they can be provided when running the workflow:

```console
$ RIBOVIZ_SAMPLES=<SAMPLES_DIRECTORY> \
  RIBOVIZ_ORGANISMS=<ORGANISMS_DIRECTORY> \
  RIBOVIZ_DATA=<DATA_DIRECTORY> \
  nextflow run prep_riboviz.nf -ansi-log false -params-file <CONFIG_FILE>
```

where:

* `<SAMPLES_DIRECTORY>` is a directory with input files.
* `<ORGANISMS_DIRECTORY>` is a directory with input files.
* `<DATA_DIRECTORY>` is a directory with input files.

The workflow will then execute, displaying information on each step as it is executed:

* Indexing steps are labelled with the index prefix
* Sample-specific steps are labelled with the sample ID.
* Multiplexed file-specific steps are labelled with the file name (minus extension).
* `collateTpms` is labelled with the IDs of all the samples that are collated.

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

---

## Help

Usage and configuration information can be viewed via use of the `--help` flag:

```console
$ nextflow run prep_riboviz.nf --help
```

Note that `--help` displays RiboViz-specific workflow help, whereas `-help` display's the `nextflow run` command's in-built help.

---

## Incremental build

If processing of a sample fails then the workflow has been written to ensure that this does not prevent the processing of other samples. For example, if the file for the vignette sample `WTnone` was corrupt in some way:

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

If the workflow fails, then a `-resume` option allows it to be rerun from the point at which it failed (for example, after fixing the issue with the file for `WTnone`):

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false -resume
```

This feature also supports incremental build. For example, given a `vignette_config.yaml` which specifies only sample `WTnone`, running Nextflow gives:

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

If `WT3AT` is then added to `vignette_config.yaml` and Nextflow is run with the `-resume` option, then only the processing for `WT3AT` is done, the cached outputs to date for `WTnone` being reused, not recomputed:

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

As another example, if the configuration file had:

```yaml
make_bedgraph: FALSE
```

then bedgraphs will not be created when the workflow is run. If you decide you do want the bedgraphs you can use `-resume` in conjunction with setting the `make_bedgraph` parameter to `TRUE`, on the command-line, to create the bedgraphs. For example:

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

When Nextflow runs it prints information about the execution of each workflow task and also logs information to a log file, `nextflow.log`. If an individual workflow task fails then information on the failure will be both printed and logged. For example:

```
...
[02/fbeb79] Submitted process > collateTpms (WT3AT, WTnone)
Error executing process > 'collateTpms (WT3AT, WTnone)'

Caused by:
  Process `collateTpms (WT3AT, WTnone)` terminated with an error exit status (1)

Command executed:

  Rscript --vanilla /home/ubuntu/riboviz/rscripts/collate_tpms.R
  --sample-subdirs=False             --output-dir=.
  --tpms-file=TPMs_collated.tsv             WT3AT WTnone

Command exit status:
  1

Command output:
  [1] "Created by: RiboViz"
  [1] "Date: 2021-04-22 04:11:03"
  [1] "File: /home/ubuntu/riboviz/rscripts/collate_tpms.R"
  [1] "Version: commit 2c506e07799c17c09fc7f0e151334c0a3313a51c date 2021-04-22 10:14:44 GMT"
  [1] "collate_tpms.R running with parameters:"
  $options
  $options$output_dir
  [1] "."

  $options$tpms_file
  [1] "TPMs_collated.tsv"

  $options$sample_subdirs
  [1] FALSE

  $options$orf_fasta
  [1] NA

  $options$help
  [1] FALSE


  $args
  [1] "WT3AT"  "WTnone"

  [1] "Loading ORFs from: ./WT3AT_tpms.tsv"
  [1] "Loading TPMs from: ./WT3AT_tpms.tsv"
  [1] "Loading TPMs from: ./WTnone_tpms.tsv"

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
  Error: Can't recycle `ORF` (size 68) to match `WTnone` (size 2).
  Backtrace:

    1. ├─global::CollateTpms(...)
    2. │ ├─`%>%`(...)
    3. │ │ └─base::eval(lhs, parent, parent)
    4. │ │   └─base::eval(lhs, parent, parent)
    5. │ └─global::MakeTpmTable(output_dir, sample_subdirs, orf_fasta, samples)
    6. │   └─dplyr::bind_cols(ORF = orfs, tpm_list[non_null_elts])
    7. │     └─vctrs::vec_cbind(!!!dots, .name_repair = .name_repair)
    8. └─vctrs::stop_incompatible_size(...)
    9.   └─vctrs:::stop_incompatible(...)
   10.     └─vctrs:::stop_vctrs(...)
  Execution halted

Work dir:
  /home/ubuntu/riboviz/work/02/fbeb79153e89224250d7a8adc507fd

Tip: you can try to figure out what's wrong by changing to the process
work dir and showing the script file named `.command.sh`

Workflow finished! (failed)
```

Note the reference to the `Work dir`. Each invocation of a task has its own sub-directory within the Nextflow `work/` directory, which includes a bash script with the command that was run (`.command.sh`), its input files, its output files (if any), a file with a log of the output printed by the command (`.command.log`) and a file with its exit code (`.exit_code`). These can be used to help you understand why a task failed and also to rerun the task in isolation.

For example, for the failure above, the task's `work/` subdirectory includes:

```console
$ ls -1A /home/ubuntu/riboviz/work/02/fbeb79153e89224250d7a8adc507fd
.command.begin
.command.err
.command.log
.command.out
.command.run
.command.sh
.exitcode
TPMs_collated.tsv
WT3AT_tpms.tsv
WTnone_tpms.tsv
```

The command that was run was:

```console
$ cat work/02/fbeb79153e89224250d7a8adc507fd/.command.sh
#!/bin/bash -ue
Rscript --vanilla /home/ubuntu/riboviz/rscripts/collate_tpms.R
--sample-subdirs=False             --output-dir=.
--tpms-file=TPMs_collated.tsv             WT3AT WTnone
```

collate_tpms.R reads the `WTnone_tpms.tsv` and `WT3AT_tpms.tsv` input files and produces the output file `TPMS_collated.tsv`. In this run, `collate_tpms.R` produces an output file with no data, only a provenance, as the error arose during execution, before the data itself was output:

```console
$ head work/02/fbeb79153e89224250d7a8adc507fd/TPMs_collated.tsv
# Created by: RiboViz
# Date: 2021-04-22 04:11:03
# File: /home/ubuntu/riboviz/rscripts/collate_tpms.R
# Version: commit 2c506e07799c17c09fc7f0e151334c0a3313a51c date 2021-04-22 10:14:44 GMT
```

The command's exit code was:

```console
$ cat work/02/fbeb79153e89224250d7a8adc507fd/.exitcode
```

A non-zero exit code means that the command failed.

When the command ran it printed the following to standard output and standard error:

```console
$ cat work/02/fbeb79153e89224250d7a8adc507fd/.command.log
[1] "Created by: RiboViz"
[1] "Date: 2021-04-22 04:11:03"
[1] "File: /home/ubuntu/riboviz/rscripts/collate_tpms.R"
[1] "Version: commit 2c506e07799c17c09fc7f0e151334c0a3313a51c date 2021-04-22 10:14:44 GMT"
[1] "collate_tpms.R running with parameters:"
$options
$options$output_dir
[1] "."

$options$tpms_file
[1] "TPMs_collated.tsv"

$options$sample_subdirs
[1] FALSE

$options$orf_fasta
[1] NA

$options$help
[1] FALSE


$args
[1] "WT3AT"  "WTnone"

[1] "Loading ORFs from: ./WT3AT_tpms.tsv"
Parsed with column specification:
cols(
  ORF = col_character(),
  readcount = col_double(),
  rpb = col_double(),
  tpm = col_double()
)
[1] "Loading TPMs from: ./WT3AT_tpms.tsv"
Parsed with column specification:
cols(
  ORF = col_character(),
  readcount = col_double(),
  rpb = col_double(),
  tpm = col_double()
)
[1] "Loading TPMs from: ./WTnone_tpms.tsv"
Parsed with column specification:
cols(
  ORF = col_character(),
  readcount = col_double(),
  rpb = col_double(),
  tpm = col_double()
)
Error: Can't recycle `ORF` (size 68) to match `WTnone` (size 2).
Backtrace:
     █
  1. ├─global::CollateTpms(...)
  2. │ ├─`%>%`(...)
  3. │ │ └─base::eval(lhs, parent, parent)
  4. │ │   └─base::eval(lhs, parent, parent)
  5. │ └─global::MakeTpmTable(output_dir, sample_subdirs, orf_fasta, samples)
  6. │   └─dplyr::bind_cols(ORF = orfs, tpm_list[non_null_elts])
  7. │     └─vctrs::vec_cbind(!!!dots, .name_repair = .name_repair)
  8. └─vctrs::stop_incompatible_size(...)
  9.   └─vctrs:::stop_incompatible(...)
 10.     └─vctrs:::stop_vctrs(...)
Execution halted
```

As the Nextflow `work/` subdirectory includes both a bash script with the command that was run and symbolic links to any input files used by the task, you can rerun the task within the subdirectory to investigate a failure in more detail. For example:

```console
$ cd work/02/fbeb79153e89224250d7a8adc507fd/
$ bash .command.sh
[1] "Created by: RiboViz"
[1] "Date: 2021-04-22 04:19:24"
[1] "File: /home/ubuntu/riboviz/rscripts/collate_tpms.R"
[1] "Version: commit 2c506e07799c17c09fc7f0e151334c0a3313a51c date 2021-04-22 10:14:44 GMT"
[1] "collate_tpms.R running with parameters:"
$options
$options$output_dir
[1] "."

$options$tpms_file
[1] "TPMs_collated.tsv"

$options$sample_subdirs
[1] FALSE

$options$orf_fasta
[1] NA

$options$help
[1] FALSE


$args
[1] "WT3AT"  "WTnone"

[1] "Loading ORFs from: ./WT3AT_tpms.tsv"
Parsed with column specification:
cols(
  ORF = col_character(),
  readcount = col_double(),
  rpb = col_double(),
  tpm = col_double()
)
[1] "Loading TPMs from: ./WT3AT_tpms.tsv"
Parsed with column specification:
cols(
  ORF = col_character(),
  readcount = col_double(),
  rpb = col_double(),
  tpm = col_double()
)
[1] "Loading TPMs from: ./WTnone_tpms.tsv"
Parsed with column specification:
cols(
  ORF = col_character(),
  readcount = col_double(),
  rpb = col_double(),
  tpm = col_double()
)
Error: Can't recycle `ORF` (size 68) to match `WTnone` (size 2).
Backtrace:
     █
  1. ├─global::CollateTpms(...)
  2. │ ├─`%>%`(...)
  3. │ │ └─base::eval(lhs, parent, parent)
  4. │ │   └─base::eval(lhs, parent, parent)
  5. │ └─global::MakeTpmTable(output_dir, sample_subdirs, orf_fasta, samples)
  6. │   └─dplyr::bind_cols(ORF = orfs, tpm_list[non_null_elts])
  7. │     └─vctrs::vec_cbind(!!!dots, .name_repair = .name_repair)
  8. └─vctrs::stop_incompatible_size(...)
  9.   └─vctrs:::stop_incompatible(...)
 10.     └─vctrs:::stop_vctrs(...)
Execution halted
```

In this example, inspecting the input files reveals the problem:

```console
$ head work/02/fbeb79153e89224250d7a8adc507fd/WT3AT_tpms.tsv
...
ORF	readcount	rpb	tpm
YAL001C	19	0.00538243626062323	953.022776609483
YAL002W	8	0.00206611570247934	365.829752221778
YAL003W	379	0.567365269461078	100458.602437955
YAL005C	2668	1.3502024291498	239069.002530875
YAL007C	15	0.0215827338129496	3821.47338292102
...

$ head work/02/fbeb79153e89224250d7a8adc507fd/WTnone_tpms.tsv
...
ORF	readcount	rpb	tpm
YAL001C	4	0.00113314447592068	116.685976611765
YAL002W	8	0.00206611570247934	212.759037933641
```

`collate_tpms.R` expects each input file to have the same number of rows and the same ORF values. In this case they differ. This would indicate a problem with the upstream task in the workflow that produced these files.

For more information on the `work/` directory, and its files, see [Nextflow `work/` directory](./prep-riboviz-operation.md#nextflow-work-directory),

For more information on Nextflow's log files, see [Nextflow workflow log files](./prep-riboviz-operation.md#nextflow-workflow-log-files).

---

## Generating reports

Nextflow's `-with-report`, `-with-timeline`, `with-trace` and `-with-dag` flags allow you to request that Nextflow create reports on a run and an image of the task execution workflow. For example:

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
