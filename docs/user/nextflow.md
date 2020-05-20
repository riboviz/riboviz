# RiboViz Nextflow workflow

[Nextflow](https://www.nextflow.io/) is a workflow management system. The file `prep_riboviz.nf` is a port of the RiboViz Python workflow, `riboviz.tools.prep_riboviz`, into a Nextflow workflow. In a future release, the RiboViz Python workflow will be deprecated by this Nextflow workflow.

For more information on Nextflow, see their [Documentation](https://www.nextflow.io/docs/latest/index.html).

This page descibes how to install Nextflow and run the RiboViz Nextflow workflow.

---

## About these instructions

These instructions were written for Ubuntu 18.04 and CentOS 7.4. Other Linux flavours will require different commands to be run.

Other versions of the prerequisites, different from the versions stated, may be usable

Only minimal installation instructions are given for each prerequisite. See the documentation for each prerequisite for full instructions.

Installing some of these tools requires you to have permission to run `sudo` to install and configure software. If you don't have `sudo` access you will have to ask a local system administrator to run these commands for you.

It is assumed you have already installed RiboViz and its preprequisites as described in [Install prerequisited](./install.md).

---

## Install Nextflow

Web sites:

* [Nextflow](https://www.nextflow.io/)
* [Documentation](https://www.nextflow.io/docs/latest/index.html)
* [GitHub](https://github.com/nextflow-io/nextflow)

### Install Nextflow using conda (recommended)

Create a new conda environment from the current one with the RiboViz dependencies and activate it:

```console
$ conda create --name riboviz-nextflow --clone base
$ conda activate riboviz-nextflow
```

Install Nextflow and its dependencies (including Java):

```console
$ conda install -c bioconda nextflow
```

Check install:

```console
$ javac -version
$ java -version
$ nextflow -version

      N E X T F L O W
      version 20.01.0 build 5264
      created 12-02-2020 10:14 UTC (02:14 PDT)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

### Install Nextflow (alternative)

Install [OpenJDK](https://openjdk.java.net) 1.8:

* Ubuntu 18 users:

```console
$ sudo apt-get install -y openjdk-8-jdk-headless
```

* CentOS 7 users:

```console
$ sudo yum install -y openjdk-8-jdk-headless
```

Check install:

```console
$ javac -version
$ java -version
```

Install Nextflow:

```console
$ curl -s https://get.nextflow.io | bash
$ export PATH=$HOME/nextflow:$PATH
$ nextflow -version

      N E X T F L O W
      version 20.01.0 build 5264
      created 12-02-2020 10:14 UTC (02:14 PDT)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

Set `PATH`:

```console
$ export PATH=$HOME/nextflow:$PATH
```

Update the `setenv.sh` script (see [Create `setenv.sh` to configure Hisat2 and Bowtie paths](https://github.com/riboviz/riboviz/blob/master/docs/user/install.md#create-setenvsh-to-configure-hisat2-and-bowtie-paths)):

```
export PATH=$HOME/nextflow:$PATH
```

---

## Run Nextflow "hello" example

```console
$ nextflow run hello
N E X T F L O W  ~  version 20.01.0
Pulling nextflow-io/hello ...
downloaded from https://github.com/nextflow-io/hello.git
Launching `nextflow-io/hello` [spontaneous_magritte] - revision: 1d43afc0ec [master]
WARN: The use of `echo` method is deprecated
executor >  local (4)
[1d/bb459e] process > sayHello [100%] 4 of 4 /
Hola world!

Bonjour world!

Ciao world!

Hello world!
```

This runs Nextflow workflow [main.nf](https://github.com/nextflow-io/hello/blob/master/main.nf) from [nextflow-io/hello.git](https://github.com/nextflow-io/hello.git).

---

## Configuring the Nextflow workflow

The Nextflow workflow uses the same YAML configuration file as `riboviz.tools.prep_riboviz` (as described in [Configuring the RiboViz workflow](./prep-riboviz-config.md).

---

## Running the Nextflow workflow

Run:

```console
$ nextflow run prep_riboviz.nf -params-file <CONFIG_FILE>  -ansi-log false
```

where:

* `<CONFIG_FILE>`: path to a YAML configuration file.
* `-ansi-log false`: requests that each invocation of a Nextflow task is displayed on a separate line.

For example, to run the vignette:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.01.0
Launching `prep_riboviz.nf` [goofy_jones] - revision: 0e94b1fa62
No such file (NotHere): example_missing_file.fastq.gz
[93/e96307] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[f1/38397f] Submitted process > cutAdapters (WT3AT)
[1a/a01841] Submitted process > buildIndicesrRNA (yeast_rRNA)
[f0/5779d2] Submitted process > cutAdapters (WTnone)
[c2/497da3] Submitted process > hisat2rRNA (WTnone)
[14/479d3b] Submitted process > hisat2rRNA (WT3AT)
[cd/2e10ec] Submitted process > hisat2ORF (WTnone)
[02/1fef3f] Submitted process > trim5pMismatches (WTnone)
[9e/90df9d] Submitted process > samViewSort (WTnone)
[1e/1a3f4d] Submitted process > outputBams (WTnone)
[d8/3d85de] Submitted process > makeBedgraphs (WTnone)
[26/055ef0] Submitted process > bamToH5 (WTnone)
[d6/355a50] Submitted process > hisat2ORF (WT3AT)
[90/78ca31] Submitted process > trim5pMismatches (WT3AT)
[03/dfb679] Submitted process > samViewSort (WT3AT)
[85/ccbb35] Submitted process > outputBams (WT3AT)
[75/7e865e] Submitted process > bamToH5 (WT3AT)
[32/a55251] Submitted process > makeBedgraphs (WT3AT)
[d5/7a7e8b] Submitted process > generateStatsFigs (WTnone)
[f4/5188e3] Submitted process > generateStatsFigs (WT3AT)
Finished processing sample: WTnone
[c1/184176] Submitted process > renameTpms (WTnone)
Finished processing sample: WT3AT
[6f/716c9e] Submitted process > renameTpms (WT3AT)
[c6/bd5296] Submitted process > collateTpms (WTnone, WT3AT)
[a6/3cb086] Submitted process > countReads
Workflow finished! (OK)
```

Each sample-specific process is labelled with the sample ID, indexing processes are labelled with the index prefix, and multiplexed file-specific processes (not shown above, but an example is below) are labelled with the file name (minus extension)

The `collateTpms` process displays the names of all the samples that are collated.

---

## Incremental build

If processing of a sample fails then `prep_riboviz.nf` has been written to ensure that this does not prevent the processing of other samples. For example, if the file for the vignette sample `WTnone` was corrupt in some way:

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

If a workflow fails, then a `-resume` option allows it to be rerun from the point at which it failed (for example, after fixing the issue with the file for `WTnone`):

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

## Differences in temporary files

Please note the following differences in the files created to those of `riboviz.tools.prep_riboviz` as documented in [What the RiboViz workflow does](docs/user/prep-riboviz-operation.md).

Two additional sample-specific temporary files are produced:

* `orf_map_clean.bam`: BAM file equivalent of `orf_map_clean.sam`, ORF-mapped reads with mismatched nt trimmed.
* `orf_map_clean.bam.bai`: BAM index file for the above.

These replace:

* `pre_dedup.bam`: BAM file prior to deduplication.
* `pre_dedup.bam.bai`: BAM index file for `pre_dedup.bam`.

If deduplication is not enabled (if `dedup_umis: FALSE`) then these are copied to become the output files:

* `<SAMPLE_ID>.bam`: BAM file of reads mapped to transcripts, which can be directly used in genome browsers.
* `<SAMPLE_ID>.bam.bai`: BAM index file for `<SAMPLE_ID>.bam`.

If deduplication is enabled (if `dedup_umis: TRUE`) the following sample-specific temporary files are also produced:

* `dedup.bam`: BAM file post deduplication.
* `dedup.bam.bai`: BAM index file for the above.

These are copied to become the output files:

* `<SAMPLE_ID>.bam`: BAM file of reads mapped to transcripts, which can be directly used in genome browsers.
* `<SAMPLE_ID>.bam.bai`: BAM index file for `<SAMPLE_ID>.bam`.

---

## Nextflow's `work/` directory

Nextflow creates a `work/` directory with all the files created during execution of the workflow. Every invocation of a task - every process - has its own subdirectory within Nextflow's `work/` directory named after the process identifiers (e.g. `ad/1e7c54`). These subdirectories have:

* Input files. These are symbolic links to the input files for the task which, depending on the task, can be:
  - Output files in other `work/` subdirectories. For example, the directory for an `hisat2rRNA` proces will have input files which are symbolic links to the output files produced by a `cutAdapters` process,
  - Input files for the workflow. For example, the directory for a `cutAdapters` process will have an input file which is a symbolic link to a sample file in `vignettte/input`.
* Output files, from the invocation of the task.
* Standard output (`.command.out`) and standard error (`.command.err`) files and exit codes (`.exitcode`) for each process.
* Bash script (`.command.sh`) containing the specific command invoked by Nextflow by that process.

For example, for the process `ad/1e7c54`, an invocation of task `hisat2rRNA` for sample `WTnone`, the `work/` directory includes:

```console
$ find work/ad/1e7c54a889f21451cb07d29655e0be/ -printf '%P\t%l\n' | sort
.command.begin
.command.err
.command.log
.command.out
.command.run
.command.sh
.command.trace
.exitcode
nonrRNA.fq
rRNA_map.sam
trim.fq	/home/ubuntu/riboviz/work/5b/ff17ba3e19c5d659e54b04b08dab85/trim.fq
yeast_rRNA.1.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.1.ht2
yeast_rRNA.2.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.2.ht2
yeast_rRNA.3.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.3.ht2
yeast_rRNA.4.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.4.ht2
yeast_rRNA.5.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.5.ht2
yeast_rRNA.6.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.6.ht2
yeast_rRNA.7.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.7.ht2
yeast_rRNA.8.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.8.ht2
```

The `.ht2` files are symbolic links to the outputs of process `e5/ccf3e6`, an invocation of task `buildIndicesrRNA`.

`prep_riboviz.nf` uses Nextflow's [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) directive which allows files to be published to specific directories outwith `work/`. By default, the files in this directory are symlinked to those in `work/`. `prep_riboviz.nf` uses this to publishes files to the index (`dir_index`), temporary (`dir_tmp`), and output (`dir_out`) directories specified in the configuration file.

### `Missing` files

If an optional file for `generate_stats_figs.R` is not provided within the YAML configuration file then a `Missing_<PARAM>` file (for example `Missing_features_file`) is created within the `work/` directories for the `generateStatsFigs` process. This symbolically links to a non-existent `Missing_<PARAM>` file in your current directory. This is not an issue since the files will not be passed onto `generate_stats_figs.R` and no attempt is made to use them. They are a side-effect of using the Nextflow pattern for optional inputs, [optional inputs](https://github.com/nextflow-io/patterns/blob/master/optional-input.nf).

---

## Debugging and bash scripts

As each `work/` subdirectory has symbolic links to any input files it requires,plus a bash script with the command that was run, that specific step can be run standalone. This can be useful for debugging. For example:

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

`-with-report`, `-with-timeline` and `-with-dag` flags allow you to reques that Nextflow create reports on a run and an image of the task execution workflow. For example:

```console
$ nextflow run prep_riboviz.nf \
    -params-file vignette/vignette_config.yaml -ansi-log false \
    -with-report report.html -with-timeline timeline.html \
    -with-dag workflow.svg
```

---

## Run vignette regression test suite

There is a separate regression test suite for Nextflow, to reflect the different files produced. The syntax is is identical to that documented in [Run vignette regression test suite](docs/developer/testing.md#run-vignette-regression-test-suite), except that:

* A extra `--nextflow` flag needs to be provided.
* The `--skip-workflow` flag must be provided i.e. you must run Nextflow yourself.

Please note that *existing* regression test data, produced by `riboviz.tools.prep_riboviz` *can* be used with the outputs from Nextflow, but that tests referencing the additional files output by Nextflow, documented above, will fail.
