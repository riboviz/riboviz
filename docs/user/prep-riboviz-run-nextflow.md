# Running the RiboViz Nextflow workflow

This page describes how to run the Nextflow workflow, `prep_riboviz.nf`. It is assumed you are familiar with [Configuring the RiboViz workflow](./prep-riboviz-config.md).

Contents:

* [Prepare input data](#prepare-input-data)
* [Set up your environment](#set-up-your-environment)
* [Configure number of processes (optional)](#configure-number-of-processes-optional)
* [Run the Nextflow workflow](#run-the-nextflow-workflow)
  - [Troubleshooting: `samtools sort: couldn't allocate memory for bam_mem`](#troubleshooting-samtools-sort-couldnt-allocate-memory-for-bam_mem)
* [Help](#help)
* [Incremental build](#incremental-build)
* [Multiplexed files](#multiplexed-files)
* [Debugging and bash scripts](#debugging-and-bash-scripts)
* [Generating reports](#generating-reports)
* [Managing your disk usage](#managing-your-disk-usage)

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

## Run the Nextflow workflow

Run:

```console
$ nextflow run prep_riboviz.nf -params-file <CONFIG_FILE> -ansi-log false
```

where:

* `<CONFIG_FILE>`: path to a YAML configuration file.
* `-ansi-log false`: requests that each invocation of a Nextflow task is displayed on a separate line.

Configuration parameters can also be provided via the command-line in the form `--<PARAMETER>=<VALUE>` (for example `--make_bedgraph=FALSE`).

Information on the key steps during processing is displayed.

Each sample-specific process is labelled with the sample ID, indexing processes are labelled with the index prefix, and multiplexed file-specific processes are labelled with the file name (minus extension).

The `collateTpms` process displays the names of all the samples that are collated.

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

Edit `prep_riboviz.nf` and change the lines:

```
	samtools view -b ${sample_sam} | samtools sort \
            -@ ${params.num_processes} -O bam -o orf_map_clean.bam -
```

to include the `samtools` flag `-m <MEMORY_DIV_PROCESSES>M` e.g.:

```
	samtools view -b ${sample_sam} | samtools sort -m 256M \
            -@ ${params.num_processes} -O bam -o orf_map_clean.bam -
```

---

## Help

Usage and configuration information can be viewed via use of the `--help` flag:

```console
$ nextflow run prep_riboviz.nf --help
```

Note that `--help` displays workflow help, whereas `-help` display's the `nextflow run` command's in-built help.

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

## Debugging and bash scripts

As described in [Nextflow `work/` directory](./prep-riboviz-operation#nextflow-work-directory), each Nextflow `work/` subdirectory has symbolic links to any input files it requires, plus a bash script with the command that was run, that specific step can be run standalone. This can be useful for debugging. For example:

```console
$ cd work/37/b11ee1d2fb315a1b72adb65c11b44/
$ cat .command.sh 
#!/bin/bash -ue
hisat2 --version
hisat2 -p 1 -N 1 -k 1 --un nonrRNA.fq -x yeast_rRNA -S rRNA_map.sam -U trim.fq
$ bash .command.sh 
/home/ubuntu/hisat2-2.1.0/hisat2-align-s version 2.1.0
64-bit
Built on login-node03
Wed Jun  7 15:53:42 EDT 2017
Compiler: gcc version 4.8.2 (GCC) 
Options: -O3 -m64 -msse2 -funroll-loops -g3 -DPOPCNT_CAPABILITY
Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
952343 reads; of these:
  952343 (100.00%) were unpaired; of these:
    467194 (49.06%) aligned 0 times
    485149 (50.94%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
50.94% overall alignment rate
```

---

## Generating reports

Nextflow's `-with-report`, `-with-timeline` and `-with-dag` flags allow you to reques that Nextflow create reports on a run and an image of the task execution workflow. For example:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false \
    -with-report report.html -with-timeline timeline.html \
    -with-dag workflow.svg
```

---

## Managing your disk usage

The workflow generates many intermediate files and some of these may be unompressed and **large**, i.e. about the same size as the input files. All these files are placed in a temporary directory (`dir_tmp`). The temporary directory's contents can be inspected for troubleshooting, if necessary.

For example, here is the volume of the outputs from a run of the vignette as documented in [Map mRNA and ribosome protected reads to transcriptome and collect data into an HDF5 file](./run-vignette.md):

| Directory |   MB |
| --------- | ---- |
| `index`   |    9 |
| `tmp`     | 1040 |
| `output`  |    3 |
| Total     | 1052 |

**Tip:** We recommend you delete temporary directories and log directories when you have completed your analyses to your satisfaction.
